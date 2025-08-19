//=============================================================================================
// Adafruit_Madgwick.h  (Chapter 7–style header, aligned with X-IO / FusionAhrs)
//=============================================================================================
//
// Orientation estimation using Madgwick's AHRS algorithm with the practical
// improvements described in Chapter 7 of S.O.H. Madgwick’s thesis:
//
//  • Ramped gain during initialisation
//  • Optional angular-rate recovery when gyro range is exceeded
//  • Accelerometer / magnetometer validity gating and recovery triggers
//  • Internal state & flags reporting
//  • Gravity / linear / Earth-frame acceleration helpers
//
// This header keeps the familiar Adafruit API (Euler getters, quaternion I/O,
// update/IMU overloads) but adds a settings struct and diagnostics mirroring
// the X-IO FusionAhrs interface.
//
//=============================================================================================

#ifndef __Adafruit_Madgwick_h__
#define __Adafruit_Madgwick_h__

#include <float.h>
#include <math.h>

#include "Adafruit_AHRS_FusionInterface.h"

// ------------------------------------
// Chapter 7–style constants (X-IO)
#define INITIAL_GAIN (10.0f)          // aggressive startup gain
#define INITIALISATION_PERIOD (3.0f)  // seconds to ramp down to nominal gain

// ------------------------------------
// Small helper types (X-IO–like, but minimal)
struct AM_Quaternion
{
  float w, x, y, z;
};

struct AM_Vector3
{
  float x, y, z;
};

// Coordinate-frame convention (match FusionConvention*)
enum AM_Convention : uint8_t
{
  AM_Convention_NWU = 0,
  AM_Convention_ENU = 1,
  AM_Convention_NED = 2
};

// Settings exposed in Chapter 7 implementation (match FusionAhrsSettings)
struct AM_Settings
{
  AM_Convention convention = AM_Convention_NWU;

  // “Nominal” gain for steady-state (Chapter 7 uses ramp from INITIAL_GAIN)
  // If gain == 0, accel/mag rejection features are disabled.
  float gain = 0.5f;

  // Gyro range in deg/s; if 0 -> disabled (no angular-rate recovery).
  // If non-zero, exceeding ~0.98 * range triggers a soft reset + recovery.
  float gyroscopeRange = 0.0f;

  // Rejection thresholds in DEGREES of error (0 disables -> treated as
  // FLT_MAX). Internally stored as squared(0.5 * sin(threshold_rad)).
  float accelerationRejection = 90.0f;
  float magneticRejection = 90.0f;

  // Number of update periods that must persist for “recovery” un-gating.
  // 0 disables recovery tracking.
  int recoveryTriggerPeriod = 0;
};

// Diagnostic “internal states” (match FusionAhrsInternalStates)
struct AM_InternalStates
{
  float accelerationErrorDegrees = 0.0f;  // asin(2|halfAccelFeedback|) in deg
  bool accelerometerIgnored = false;
  float accelerationRecoveryRatio = 0.0f;  // 0..1 of trigger / period

  float magneticErrorDegrees = 0.0f;  // asin(2|halfMagFeedback|) in deg
  bool magnetometerIgnored = false;
  float magneticRecoveryRatio = 0.0f;  // 0..1 of trigger / period
};

// Status flags (match FusionAhrsFlags)
struct AM_Flags
{
  bool initialising = false;
  bool angularRateRecovery = false;
  bool accelerationRecovery = false;
  bool magneticRecovery = false;
};

//=============================================================================================
// Class
//=============================================================================================

class Adafruit_Madgwick : public Adafruit_AHRS_FusionInterface
{
 public:
  // -----------------------------
  // Construction & configuration
  Adafruit_Madgwick() : Adafruit_Madgwick(0.1f) {}
  explicit Adafruit_Madgwick(float beta_nominal)
  {
    // Back-compat: allow constructing with a classic “beta”
    settings.gain = beta_nominal;
    resetChapter7State(/*keepSettings=*/false);
  }

  // Set sampling frequency; used by overloads that don’t pass dt explicitly
  void begin(float sampleFrequency)
  {
    invSampleFreq = (sampleFrequency > 0.0f) ? (1.0f / sampleFrequency) : 0.0f;
  }

  // Chapter 7 API: settings
  void setSettings(const AM_Settings& s);
  AM_Settings getSettings() const { return settings; }

  // Chapter 7 API: reset (keeps settings)
  void reset() { resetChapter7State(/*keepSettings=*/true); }

  // --------------------------------
  // Core update APIs (back compatible)
  void update(float gx, float gy, float gz, float ax, float ay, float az,
              float mx, float my, float mz)
  {
    update(gx, gy, gz, ax, ay, az, mx, my, mz, invSampleFreq);
  }

  void updateIMU(float gx, float gy, float gz, float ax, float ay, float az)
  {
    updateIMU(gx, gy, gz, ax, ay, az, invSampleFreq);
  }

  // Chapter 7–style update (explicit dt)
  void update(float gx, float gy, float gz, float ax, float ay, float az,
              float mx, float my, float mz, float dt);

  void updateIMU(float gx, float gy, float gz, float ax, float ay, float az,
                 float dt);

  // Optional: heading “reset” (useful when running 6DoF only)
  void setHeadingDegrees(float headingDeg);

  // --------------------------------
  // Quaternion access (back compatible)
  void getQuaternion(float* w, float* x, float* y, float* z) const
  {
    *w = q0;
    *x = q1;
    *y = q2;
    *z = q3;
  }
  void setQuaternion(float w, float x, float y, float z)
  {
    q0 = w;
    q1 = x;
    q2 = y;
    q3 = z;
    anglesComputed = false;
  }

  // Euler outputs (back compatible)
  float getRoll()
  {
    ensureAngles();
    return roll * 57.29578f;
  }
  float getPitch()
  {
    ensureAngles();
    return pitch * 57.29578f;
  }
  float getYaw()
  {
    ensureAngles();
    return yaw * 57.29578f + 180.0f;
  }
  float getRollRadians()
  {
    ensureAngles();
    return roll;
  }
  float getPitchRadians()
  {
    ensureAngles();
    return pitch;
  }
  float getYawRadians()
  {
    ensureAngles();
    return yaw;
  }

  // Gravity vector in body frame (3rd column of R^T)
  void getGravityVector(float* x, float* y, float* z)
  {
    ensureAngles();
    *x = grav[0];
    *y = grav[1];
    *z = grav[2];
  }

  // Chapter 7 helpers: linear acceleration (body) = accel - gravity (g)
  AM_Vector3 getLinearAccelerationBody() const;

  // Chapter 7 helpers: Earth-frame acceleration with gravity removed (g)
  AM_Vector3 getEarthAcceleration() const;

  // Chapter 7 diagnostics
  AM_InternalStates getInternalStates() const;
  AM_Flags getFlags() const;

  // --------------------------------
  // Legacy beta access (maps to settings.gain)
  float getBeta() const { return settings.gain; }
  void setBeta(float beta)
  {
    settings.gain = beta;
    if (!initialising) rampedGain = settings.gain;
  }

 private:
  // ================================
  // Core state (classic)
  static float invSqrt(float x);
  float invSampleFreq = 0.0f;

  // Quaternion of sensor (body) relative to Earth
  float q0 = 1.0f, q1 = 0.0f, q2 = 0.0f, q3 = 0.0f;

  // Cached Euler & gravity
  mutable bool anglesComputed = false;
  mutable float roll = 0.0f, pitch = 0.0f, yaw = 0.0f;
  mutable float grav[3] = {0.0f, 0.0f, 1.0f};
  void computeAngles();
  void ensureAngles()
  {
    if (!anglesComputed) computeAngles();
  }

  // ================================
  // Chapter 7 settings & state
  AM_Settings settings;  // user settings (Chapter 7)

  // Initialisation ramp
  bool initialising = true;
  float rampedGain = INITIAL_GAIN;  // current gain during ramp
  float rampedGainStep =
      (INITIAL_GAIN - 0.5f) / INITIALISATION_PERIOD;  // updated in setSettings

  // Angular-rate recovery (reinit if gyro exceeds range)
  bool angularRateRecovery = false;

  // Latest accelerometer (for linear/Earth-frame accel outputs)
  AM_Vector3 lastAccel = {0.0f, 0.0f, 0.0f};

  // Feedback (scaled by 0.5) for diagnostics and gating
  AM_Vector3 halfAccelFeedback = {0.0f, 0.0f, 0.0f};
  AM_Vector3 halfMagFeedback = {0.0f, 0.0f, 0.0f};

  // Accel gating & recovery
  bool accelerometerIgnored = false;
  int accelerationRecoveryTrigger = 0;
  int accelerationRecoveryTimeout = 0;  // copied from settings

  // Mag gating & recovery
  bool magnetometerIgnored = false;
  int magneticRecoveryTrigger = 0;
  int magneticRecoveryTimeout = 0;  // copied from settings

  // ----------------
  // Internal helpers
  void resetChapter7State(bool keepSettings);

  // Gain ramp update per step (dt seconds)
  void stepInitialisationRamp(float dt);

  // Compute 0.5 * gravity direction in the chosen convention (R^T third column
  // * 0.5)
  AM_Vector3 halfGravity() const;

  // Compute 0.5 * magnetic-field direction in the chosen convention
  AM_Vector3 halfMagnetic() const;

  // Cross-product feedback (optionally “flipped” if dot < 0 as in Fusion)
  static AM_Vector3 feedback(const AM_Vector3& sensor,
                             const AM_Vector3& reference, bool flipWhenObtuse);

  // Basic vector ops (header-only declarations; defined in .cpp)
  static inline AM_Vector3 vadd(const AM_Vector3& a, const AM_Vector3& b)
  {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
  }
  static inline AM_Vector3 vsub(const AM_Vector3& a, const AM_Vector3& b)
  {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
  }
  static inline AM_Vector3 vscale(const AM_Vector3& a, float s)
  {
    return {a.x * s, a.y * s, a.z * s};
  }
  static inline float vdot(const AM_Vector3& a, const AM_Vector3& b)
  {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }
  static inline AM_Vector3 vcross(const AM_Vector3& a, const AM_Vector3& b)
  {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
  }
  static inline float vmag2(const AM_Vector3& a) { return vdot(a, a); }
  static inline float vmag(const AM_Vector3& a) { return sqrtf(vmag2(a)); }
  static inline AM_Vector3 vnorm(const AM_Vector3& a)
  {
    float n = vmag(a);
    return (n > 0.0f) ? vscale(a, 1.0f / n) : AM_Vector3{0, 0, 0};
  }

  // Quaternion ops used internally (defined in .cpp)
  static AM_Quaternion qmul(const AM_Quaternion& A, const AM_Quaternion& B);
  static AM_Quaternion qadd(const AM_Quaternion& A, const AM_Quaternion& B);
  static AM_Quaternion qnorm(const AM_Quaternion& q);
  static AM_Quaternion qmulVec(const AM_Quaternion& q,
                               const AM_Vector3& vHalfRad);  // (0, v) * q, etc.

  // Frame helpers
  AM_Vector3 gravityBodyFromQuat() const;  // 3rd column of R^T
  AM_Vector3 accelEarthFromBody(const AM_Vector3& aBody) const;  // R * aBody
};

#endif  // __Adafruit_Madgwick_h__