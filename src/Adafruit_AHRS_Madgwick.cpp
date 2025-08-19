//=============================================================================================
// Adafruit_Madgwick.cpp  (Chapter 7–style implementation aligned with X-IO /
// FusionAhrs)
//=============================================================================================
//
// Notes:
//  - Gyroscope inputs are in deg/s; internally converted to rad/s * 0.5.
//  - Accel is in g (unitless), magnetometer in arbitrary units (normalised
//  inside).
//  - Implements: initial gain ramp, accel/mag rejection with recovery triggers,
//    angular-rate recovery on gyro overflow, and helper outputs.
//
//=============================================================================================

#include "Adafruit_AHRS_Madgwick.h"

#include <float.h>
#include <math.h>

//------------------------------------------------------------------------------
// Chapter 7 constants (match FusionAhrs)
#ifndef INITIAL_GAIN
#define INITIAL_GAIN (10.0f)
#endif

#ifndef INITIALISATION_PERIOD
#define INITIALISATION_PERIOD (3.0f)
#endif

//------------------------------------------------------------------------------
// Local helpers (static)

static inline int ClampInt(const int value, const int minv, const int maxv)
{
  if (value < minv) return minv;
  if (value > maxv) return maxv;
  return value;
}

static inline float Deg2Rad(const float d)
{
  return d * 0.017453292519943295f;
}  // pi/180

//------------------------------------------------------------------------------
// Quaternion helpers (member-static in header, defined here)

AM_Quaternion Adafruit_Madgwick::qmul(const AM_Quaternion& A,
                                      const AM_Quaternion& B)
{
  AM_Quaternion C;
  C.w = A.w * B.w - A.x * B.x - A.y * B.y - A.z * B.z;
  C.x = A.w * B.x + A.x * B.w + A.y * B.z - A.z * B.y;
  C.y = A.w * B.y - A.x * B.z + A.y * B.w + A.z * B.x;
  C.z = A.w * B.z + A.x * B.y - A.y * B.x + A.z * B.w;
  return C;
}

AM_Quaternion Adafruit_Madgwick::qadd(const AM_Quaternion& A,
                                      const AM_Quaternion& B)
{
  return AM_Quaternion{A.w + B.w, A.x + B.x, A.y + B.y, A.z + B.z};
}

AM_Quaternion Adafruit_Madgwick::qnorm(const AM_Quaternion& q)
{
  const float n = sqrtf(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
  if (n <= 0.0f) return AM_Quaternion{1, 0, 0, 0};
  const float inv = 1.0f / n;
  return AM_Quaternion{q.w * inv, q.x * inv, q.y * inv, q.z * inv};
}

// q * (0, v)
AM_Quaternion Adafruit_Madgwick::qmulVec(const AM_Quaternion& q,
                                         const AM_Vector3& vHalfRad)
{
  AM_Quaternion V{0.0f, vHalfRad.x, vHalfRad.y, vHalfRad.z};
  return qmul(q, V);
}

//------------------------------------------------------------------------------
// Private utilities

void Adafruit_Madgwick::resetChapter7State(bool keepSettings)
{
  // Quaternion identity
  q0 = 1.0f;
  q1 = 0.0f;
  q2 = 0.0f;
  q3 = 0.0f;

  // Cached angles/grav
  anglesComputed = false;
  roll = pitch = yaw = 0.0f;
  grav[0] = 0.0f;
  grav[1] = 0.0f;
  grav[2] = 1.0f;

  // Settings & ramp
  if (!keepSettings)
  {
    settings = AM_Settings{};  // defaults
    settings.gain = 0.5f;
    settings.accelerationRejection = 90.0f;
    settings.magneticRejection = 90.0f;
    settings.gyroscopeRange = 0.0f;
    settings.recoveryTriggerPeriod = 0;
  }

  initialising = true;
  rampedGain = INITIAL_GAIN;
  rampedGainStep = (INITIAL_GAIN - settings.gain) / INITIALISATION_PERIOD;
  angularRateRecovery = false;

  // Accel/mag bookkeeping
  lastAccel = {0.0f, 0.0f, 0.0f};
  halfAccelFeedback = {0.0f, 0.0f, 0.0f};
  halfMagFeedback = {0.0f, 0.0f, 0.0f};

  accelerometerIgnored = false;
  accelerationRecoveryTrigger = 0;
  accelerationRecoveryTimeout = settings.recoveryTriggerPeriod;

  magnetometerIgnored = false;
  magneticRecoveryTrigger = 0;
  magneticRecoveryTimeout = settings.recoveryTriggerPeriod;

  // If gain==0 or recoveryTriggerPeriod==0, disable rejection features
  if ((settings.gain == 0.0f) || (settings.recoveryTriggerPeriod == 0))
  {
    // We mirror Fusion’s approach by pushing thresholds to FLT_MAX
    // (done in setSettings as well).
  }
}

void Adafruit_Madgwick::stepInitialisationRamp(float dt)
{
  if (!initialising) return;
  rampedGain -= rampedGainStep * dt;
  if ((rampedGain < settings.gain) || (settings.gain == 0.0f))
  {
    rampedGain = settings.gain;
    initialising = false;
    angularRateRecovery = false;
  }
}

// Half of gravity direction per convention (third column of R^T, scaled by 0.5)
AM_Vector3 Adafruit_Madgwick::halfGravity() const
{
  const float& w = q0;
  const float& x = q1;
  const float& y = q2;
  const float& z = q3;
  switch (settings.convention)
  {
    case AM_Convention_NWU:
    case AM_Convention_ENU:
      return AM_Vector3{
          x * z - w * y,
          y * z + w * x,
          w * w - 0.5f + z * z,
      };
    case AM_Convention_NED:
      return AM_Vector3{
          w * y - x * z,
          -(y * z + w * x),
          0.5f - w * w - z * z,
      };
  }
  return AM_Vector3{0, 0, 0};
}

// Half of magnetic-field direction per convention
AM_Vector3 Adafruit_Madgwick::halfMagnetic() const
{
  const float& w = q0;
  const float& x = q1;
  const float& y = q2;
  const float& z = q3;
  switch (settings.convention)
  {
    case AM_Convention_NWU:
      return AM_Vector3{
          x * y + w * z,
          w * w - 0.5f + y * y,
          y * z - w * x,
      };
    case AM_Convention_ENU:
      return AM_Vector3{
          0.5f - w * w - x * x,
          w * z - x * y,
          -(x * z + w * y),
      };
    case AM_Convention_NED:
      return AM_Vector3{
          -(x * y + w * z),
          0.5f - w * w - y * y,
          w * x - y * z,
      };
  }
  return AM_Vector3{0, 0, 0};
}

// Cross-product feedback; if dot < 0 use normalised cross (Fusion flip > 90°)
AM_Vector3 Adafruit_Madgwick::feedback(const AM_Vector3& sensor,
                                       const AM_Vector3& reference,
                                       bool flipWhenObtuse)
{
  const float dot = vdot(sensor, reference);
  AM_Vector3 c = vcross(sensor, reference);
  if (flipWhenObtuse && (dot < 0.0f))
  {
    const float m = vmag(c);
    if (m > 0.0f) return vscale(c, 1.0f / m);
  }
  return c;
}

//------------------------------------------------------------------------------
// Public: settings

void Adafruit_Madgwick::setSettings(const AM_Settings& s)
{
  settings.convention = s.convention;
  settings.gain = s.gain;

  // Gyro range: if 0 -> disabled, else use 0.98 * range (Fusion behaviour)
  settings.gyroscopeRange =
      (s.gyroscopeRange == 0.0f) ? FLT_MAX : (0.98f * s.gyroscopeRange);

  // Rejection thresholds: degrees -> internal metric pow(0.5 * sin(theta), 2)
  auto rejConv = [](float deg) -> float
  {
    if (deg == 0.0f) return FLT_MAX;
    const float rad = Deg2Rad(deg);
    const float val = 0.5f * sinf(rad);
    return val * val;
  };
  settings.accelerationRejection = rejConv(s.accelerationRejection);
  settings.magneticRejection = rejConv(s.magneticRejection);

  settings.recoveryTriggerPeriod = s.recoveryTriggerPeriod;

  accelerationRecoveryTimeout = settings.recoveryTriggerPeriod;
  magneticRecoveryTimeout = settings.recoveryTriggerPeriod;

  if ((s.gain == 0.0f) || (s.recoveryTriggerPeriod == 0))
  {
    settings.accelerationRejection = FLT_MAX;
    settings.magneticRejection = FLT_MAX;
  }

  if (!initialising)
  {
    rampedGain = settings.gain;
  }
  rampedGainStep = (INITIAL_GAIN - settings.gain) / INITIALISATION_PERIOD;
}

//------------------------------------------------------------------------------
// Update (with magnetometer) — mirrors FusionAhrsUpdate

void Adafruit_Madgwick::update(float gx, float gy, float gz, float ax, float ay,
                               float az, float mx, float my, float mz, float dt)
{
  // Store accelerometer
  lastAccel = {ax, ay, az};

  // Reinitialise if gyroscope range exceeded (deg/s)
  if ((fabsf(gx) > settings.gyroscopeRange) ||
      (fabsf(gy) > settings.gyroscopeRange) ||
      (fabsf(gz) > settings.gyroscopeRange))
  {
    // Soft reset while maintaining settings; preserve quaternion like Fusion
    AM_Quaternion qPrev{q0, q1, q2, q3};
    resetChapter7State(/*keepSettings=*/true);
    q0 = qPrev.w;
    q1 = qPrev.x;
    q2 = qPrev.y;
    q3 = qPrev.z;
    angularRateRecovery = true;
  }

  // Ramp down gain during initialisation
  stepInitialisationRamp(dt);

  // Calculate direction of gravity indicated by algorithm (scaled by 0.5)
  const AM_Vector3 halfG = halfGravity();

  // ---------------------------
  // Accelerometer feedback
  AM_Vector3 halfAccelFB{0, 0, 0};
  accelerometerIgnored = true;
  const bool accelIsZero = (ax == 0.0f) && (ay == 0.0f) && (az == 0.0f);
  if (!accelIsZero)
  {
    // Normalise accelerometer; form feedback vs halfGravity
    const AM_Vector3 aN = vnorm(AM_Vector3{ax, ay, az});
    halfAccelFeedback = feedback(aN, halfG, /*flipWhenObtuse=*/true);

    // Gate accelerometer if error above threshold (compare squared magnitudes)
    if (initialising ||
        (vmag2(halfAccelFeedback) <= settings.accelerationRejection))
    {
      accelerometerIgnored = false;
      accelerationRecoveryTrigger -= 9;
    }
    else
    {
      accelerationRecoveryTrigger += 1;
    }

    // Recovery window logic
    if (accelerationRecoveryTrigger > accelerationRecoveryTimeout)
    {
      accelerationRecoveryTimeout = 0;
      accelerometerIgnored = false;
    }
    else
    {
      accelerationRecoveryTimeout = settings.recoveryTriggerPeriod;
    }
    accelerationRecoveryTrigger = ClampInt(accelerationRecoveryTrigger, 0,
                                           settings.recoveryTriggerPeriod);

    // Apply accelerometer feedback if not ignored
    if (!accelerometerIgnored)
    {
      halfAccelFB = halfAccelFeedback;
    }
  }

  // ---------------------------
  // Magnetometer feedback
  AM_Vector3 halfMagFB{0, 0, 0};
  magnetometerIgnored = true;
  const bool magIsZero = (mx == 0.0f) && (my == 0.0f) && (mz == 0.0f);
  if (!magIsZero)
  {
    // Direction of magnetic field indicated by algorithm (scaled by 0.5)
    const AM_Vector3 halfM = halfMagnetic();

    // Cross(halfGravity, magnetometer), normalised, then feedback vs halfM
    const AM_Vector3 hxc = vcross(halfG, AM_Vector3{mx, my, mz});
    const AM_Vector3 hxcN = vnorm(hxc);
    halfMagFeedback = feedback(hxcN, halfM, /*flipWhenObtuse=*/true);

    // Gate magnetometer based on threshold
    if (initialising || (vmag2(halfMagFeedback) <= settings.magneticRejection))
    {
      magnetometerIgnored = false;
      magneticRecoveryTrigger -= 9;
    }
    else
    {
      magneticRecoveryTrigger += 1;
    }

    // Recovery window logic
    if (magneticRecoveryTrigger > magneticRecoveryTimeout)
    {
      magneticRecoveryTimeout = 0;
      magnetometerIgnored = false;
    }
    else
    {
      magneticRecoveryTimeout = settings.recoveryTriggerPeriod;
    }
    magneticRecoveryTrigger =
        ClampInt(magneticRecoveryTrigger, 0, settings.recoveryTriggerPeriod);

    // Apply magnetometer feedback if not ignored
    if (!magnetometerIgnored)
    {
      halfMagFB = halfMagFeedback;
    }
  }

  // ---------------------------
  // Convert gyroscope to radians/s scaled by 0.5
  const AM_Vector3 halfGyro = vscale(AM_Vector3{gx, gy, gz}, Deg2Rad(0.5f));

  // Apply feedback to gyroscope
  const AM_Vector3 fbSum = vadd(halfAccelFB, halfMagFB);
  const AM_Vector3 adjustedHalfGyro = vadd(halfGyro, vscale(fbSum, rampedGain));

  // Integrate quaternion rate: q += q ⊗ (0, adjustedHalfGyro) * dt
  AM_Quaternion q{q0, q1, q2, q3};
  q = qadd(q, qmulVec(q, vscale(adjustedHalfGyro, dt)));

  // Normalise quaternion
  q = qnorm(q);
  q0 = q.w;
  q1 = q.x;
  q2 = q.y;
  q3 = q.z;
  anglesComputed = false;
}

//------------------------------------------------------------------------------
// Update (IMU only) — simply calls update with zero magnetometer and zero
// heading during init

void Adafruit_Madgwick::updateIMU(float gx, float gy, float gz, float ax,
                                  float ay, float az, float dt)
{
  update(gx, gy, gz, ax, ay, az, 0.0f, 0.0f, 0.0f, dt);

  // Zero heading during initialisation (Fusion behaviour)
  if (initialising)
  {
    setHeadingDegrees(0.0f);
  }
}

//------------------------------------------------------------------------------
// Heading reset (useful when no magnetometer)

void Adafruit_Madgwick::setHeadingDegrees(float headingDeg)
{
  const float w = q0, x = q1, y = q2, z = q3;

  // yaw = atan2(w*z + x*y, 0.5 - y^2 - z^2)  (Fusion’s expression)
  const float yaw = atan2f(w * z + x * y, 0.5f - y * y - z * z);

  const float halfYawMinusHeading = 0.5f * (yaw - Deg2Rad(headingDeg));
  const AM_Quaternion rot{cosf(halfYawMinusHeading), 0.0f, 0.0f,
                          -sinf(halfYawMinusHeading)};

  AM_Quaternion q{q0, q1, q2, q3};
  q = qmul(rot, q);
  q = qnorm(q);
  q0 = q.w;
  q1 = q.x;
  q2 = q.y;
  q3 = q.z;
  anglesComputed = false;
}

//------------------------------------------------------------------------------
// Gravity & acceleration helpers (Chapter 7 outputs)

AM_Vector3 Adafruit_Madgwick::gravityBodyFromQuat() const
{
  // third column of the transposed rotation matrix (not scaled)
  const float& w = q0;
  const float& x = q1;
  const float& y = q2;
  const float& z = q3;
  return AM_Vector3{
      2.0f * (x * z - w * y),
      2.0f * (y * z + w * x),
      2.0f * (w * w - 0.5f + z * z),
  };
}

AM_Vector3 Adafruit_Madgwick::accelEarthFromBody(const AM_Vector3& aBody) const
{
  // R * aBody (explicit expanded form using quaternion)
  const float& w = q0;
  const float& x = q1;
  const float& y = q2;
  const float& z = q3;

  const float qwqw = w * w;
  const float qwqx = w * x;
  const float qwqy = w * y;
  const float qwqz = w * z;
  const float qxqy = x * y;
  const float qxqz = x * z;
  const float qyqz = y * z;

  AM_Vector3 e;
  e.x = 2.0f * ((qwqw - 0.5f + x * x) * aBody.x + (qxqy - qwqz) * aBody.y +
                (qxqz + qwqy) * aBody.z);
  e.y = 2.0f * ((qxqy + qwqz) * aBody.x + (qwqw - 0.5f + y * y) * aBody.y +
                (qyqz - qwqx) * aBody.z);
  e.z = 2.0f * ((qxqz - qwqy) * aBody.x + (qyqz + qwqx) * aBody.y +
                (qwqw - 0.5f + z * z) * aBody.z);

  // Remove gravity depending on convention
  switch (settings.convention)
  {
    case AM_Convention_NWU:
    case AM_Convention_ENU:
      e.z -= 1.0f;
      break;
    case AM_Convention_NED:
      e.z += 1.0f;
      break;
  }
  return e;
}

AM_Vector3 Adafruit_Madgwick::getLinearAccelerationBody() const
{
  const AM_Vector3 g = gravityBodyFromQuat();
  switch (settings.convention)
  {
    case AM_Convention_NWU:
    case AM_Convention_ENU:
      return vsub(lastAccel, g);
    case AM_Convention_NED:
      return vadd(lastAccel, g);
  }
  return AM_Vector3{0, 0, 0};
}

AM_Vector3 Adafruit_Madgwick::getEarthAcceleration() const
{
  return accelEarthFromBody(lastAccel);
}

//------------------------------------------------------------------------------
// Diagnostics (match Fusion getters)

AM_InternalStates Adafruit_Madgwick::getInternalStates() const
{
  // accelerationError = asin(2 * |halfAccelFeedback|) in degrees
  const float accelErrRad =
      asinf(fminf(1.0f, fmaxf(-1.0f, 2.0f * vmag(halfAccelFeedback))));
  const float magErrRad =
      asinf(fminf(1.0f, fmaxf(-1.0f, 2.0f * vmag(halfMagFeedback))));

  AM_InternalStates s;
  s.accelerationErrorDegrees = accelErrRad * 57.29577951308232f;
  s.accelerometerIgnored = accelerometerIgnored;
  s.accelerationRecoveryRatio = (settings.recoveryTriggerPeriod == 0)
                                    ? 0.0f
                                    : ((float)accelerationRecoveryTrigger /
                                       (float)settings.recoveryTriggerPeriod);

  s.magneticErrorDegrees = magErrRad * 57.29577951308232f;
  s.magnetometerIgnored = magnetometerIgnored;
  s.magneticRecoveryRatio = (settings.recoveryTriggerPeriod == 0)
                                ? 0.0f
                                : ((float)magneticRecoveryTrigger /
                                   (float)settings.recoveryTriggerPeriod);
  return s;
}

AM_Flags Adafruit_Madgwick::getFlags() const
{
  AM_Flags f;
  f.initialising = initialising;
  f.angularRateRecovery = angularRateRecovery;
  f.accelerationRecovery =
      (accelerationRecoveryTrigger > accelerationRecoveryTimeout);
  f.magneticRecovery = (magneticRecoveryTrigger > magneticRecoveryTimeout);
  return f;
}

//------------------------------------------------------------------------------
// Classic utilities retained from original Adafruit code

float Adafruit_Madgwick::invSqrt(float x)
{
  float halfx = 0.5f * x;
  union
  {
    float f;
    long i;
  } conv = {x};
  conv.i = 0x5f3759df - (conv.i >> 1);
  conv.f *= 1.5f - (halfx * conv.f * conv.f);
  conv.f *= 1.5f - (halfx * conv.f * conv.f);
  return conv.f;
}

void Adafruit_Madgwick::computeAngles()
{
  // Same Euler extraction as original (roll/pitch/yaw) and gravity
  roll = atan2f(q0 * q1 + q2 * q3, 0.5f - q1 * q1 - q2 * q2);
  pitch = asinf(-2.0f * (q1 * q3 - q0 * q2));
  yaw = atan2f(q1 * q2 + q0 * q3, 0.5f - q2 * q2 - q3 * q3);

  grav[0] = 2.0f * (q1 * q3 - q0 * q2);
  grav[1] = 2.0f * (q0 * q1 + q2 * q3);
  grav[2] = 2.0f * (q0 * q0 - 0.5f + q3 * q3);

  anglesComputed = true;
}