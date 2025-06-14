 #include <Arduino.h> #include <BasicLinearAlgebra.h>

// --- Configuration Constants --- namespace C { constexpr int PHASE_PINS[3] = {A0, A1, A2}; constexpr int COMMON_PIN = 9; constexpr float ELECTRICAL_STEPS_PER_REV = 6.0f; constexpr float TWO_PI = 6.28318531f; constexpr float SQRT3 = 1.73205081f; constexpr unsigned int SETTLE_US = 10;

// Calibration
constexpr int CAL_POINTS = 16;
constexpr int CAL_PULSE_US[CAL_POINTS] = {50,75,100,150,200,250,300,400,500,600,700,800,900,1000,1200,1500};
constexpr int MIN_PULSE = CAL_PULSE_US[0];
constexpr int MAX_PULSE = CAL_PULSE_US[CAL_POINTS-1];
constexpr int CAL_SAMPLES = 256;
constexpr float MAX_VEL = 15.0f;

// Filter parameters
constexpr float Q_KF = 0.5f, R_KF = 200.0f, EPS_KF = 1e-6f;
constexpr float Q_ANG = 0.0005f, Q_VEL = 0.005f, R_ANG = 0.2f, EPS_EKF = 1e-6f;
constexpr float OUTLIER_THRESH = 3.14159f;

// Motion detection
constexpr int VEL_HIST = 15, ANG_HIST = 7;
constexpr float VEL_STEP = 0.04f, MIN_CONSIST = 0.05f, NOISE_THRESH = 0.015f;
constexpr unsigned long LOCK_TIMEOUT = 250000UL;
constexpr float REV_THRESH = 2.356f;
constexpr int DIR_COUNT = 5;
constexpr float TREND_THRESH = 0.008f;

}

// Lookup table struct LUT { int pulse; float offset; }; static LUT lut[C::CAL_POINTS];

// Kalman filter matrices using Mat3 = BLA::Matrix<3,3>; using Vec3 = BLA::Matrix<3,1>; using Mat2 = BLA::Matrix<2,2>; using Vec2 = BLA::Matrix<2,1>;

static Vec3 x_kf; static Mat3 P_kf, Q_kf, R_kf, I3; static Vec2 x_ekf; static Mat2 P_ekf, Q_ekf, R_ekf, I2, P_init;

// History struct Hist { float ang, vel; unsigned long t; }; static Hist vel_hist[C::VEL_HIST], ang_hist[C::ANG_HIST];

// State static unsigned long last_nz_vel, last_ekf_us; static long rev_count, mech_steps; static float unwrap_ang, last_valid_ang, smooth_vel, vel_trend, last_stable_vel; static int vel_idx, ang_idx, consist_dir; static bool vel_full, ang_full, locked, stopping;

// Forward declarations void calibrate(); void sense(); float getAngle(); void updateVel(float, unsigned long); void updateAng(float, unsigned long); float getSmooth(); float getTrend(); bool isStopping(); bool isConsistent(float); void updateRev(float); inline float normAngle(float a) noexcept { while (a > C::TWO_PI) a -= C::TWO_PI; while (a < 0) a += C::TWO_PI; return a; }

void setup() { Serial.begin(115200); while (!Serial) delay(10); Serial.println(F("Init encoder..."));

// Initialize identity matrices
I3 = Mat3::Identity();
I2 = Mat2::Identity();

// Process noise & measurement noise
Q_kf = Mat3::Identity() * C::Q_KF;
R_kf = Mat3::Identity() * C::R_KF;
Q_ekf = Mat2::Zero(); Q_ekf(0,0) = C::Q_ANG; Q_ekf(1,1) = C::Q_VEL;
R_ekf = Mat2::Zero(); R_ekf(0,0) = C::R_ANG;

// ADC & pins
analogReference(INTERNAL);
ADCSRA = (ADCSRA & ~7) | 5;
pinMode(C::COMMON_PIN, OUTPUT);
digitalWrite(C::COMMON_PIN, LOW);

// Calibrate and init state
calibrate();
x_kf.Fill(0); P_kf = Mat3::Identity() * 100;
for (auto &h : vel_hist) h = {};
for (auto &h : ang_hist) h = {};
sense();
float init_ang = getAngle();
x_ekf(0) = init_ang; x_ekf(1) = 0;
P_ekf = Mat2::Identity();
P_init = P_ekf;
last_ekf_us = micros();
unwrap_ang = last_valid_ang = init_ang;

Serial.println(F("Ready"));

}

void ekf(float meas_ang) { // Time update unsigned long now = micros(); float dt = (now - last_ekf_us) * 1e-6f; last_ekf_us = now;

// State prediction
Vec2 pred{x_ekf(0) + x_ekf(1)*dt, x_ekf(1)};
Mat2 A{{1, dt},{0, 1}};
P_ekf = A * P_ekf * A.transpose() + Q_ekf;

// Measurement update
float y = normAngle(meas_ang - pred(0));
if (fabs(y) > C::OUTLIER_THRESH) return;

float S = P_ekf(0,0) + R_ekf(0,0) + C::EPS_EKF;
Vec2 K{P_ekf(0,0)/S, P_ekf(1,0)/S};

x_ekf(0) = normAngle(pred(0) + K(0)*y);
x_ekf(1) = pred(1) + K(1)*y;

Mat2 I = Mat2::Identity();
P_ekf = (I - BLA::Matrix<2,2>{{K(0),0},{K(1),0}}) * P_ekf;

// Update motion stats
updateVel(x_ekf(1), now);
updateAng(x_ekf(0), now);
smooth_vel = getSmooth();
vel_trend = getTrend();
stopping = isStopping();

// Lock logic
bool moving = fabs(x_ekf(1)) > C::VEL_STEP;
if (moving && isConsistent(x_ekf(1))) {
    int dir = x_ekf(1) > 0;
    int last_dir = last_stable_vel > 0;
    consist_dir = (fabs(last_stable_vel) > C::VEL_STEP && dir == last_dir) ? consist_dir+1 : 1;
    last_stable_vel = x_ekf(1);
} else {
    consist_dir = 0;
}

if (!locked && moving && consist_dir >= C::DIR_COUNT) {
    locked = true;
    Serial.println(F("Locked"));
} else if (locked) {
    float diff = normAngle(x_ekf(0) - unwrap_ang);
    unwrap_ang += diff;
    updateRev(diff);
    if (moving) last_nz_vel = now;
    if (stopping || (!moving && now - last_nz_vel > C::LOCK_TIMEOUT)) {
        locked = false;
        Serial.println(F("Unlocked"));
    }
}

}

// ... rest of functions unchanged (calibrate, sense, getAngle, etc.)

