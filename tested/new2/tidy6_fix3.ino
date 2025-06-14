#include <Arduino.h>
#include <BasicLinearAlgebra.h>

// Configuration constants
namespace C {
    constexpr int PHASE_PINS[3] = {A0, A1, A2};
    constexpr int COMMON_PIN = 9;
    constexpr float ELECTRICAL_STEPS_PER_REV = 6.0f;
    constexpr float TWO_PI = 6.28318531f;
    constexpr float SQRT3 = 1.73205081f;
    constexpr int SETTLE_US = 10;
    
    // Calibration
    constexpr int CAL_POINTS = 16;
    constexpr int CAL_PULSE_US[CAL_POINTS] = {50,75,100,150,200,250,300,400,500,600,700,800,900,1000,1200,1500};
    constexpr int MIN_PULSE = 50, MAX_PULSE = 1500;
    constexpr int CAL_SAMPLES = 256;
    constexpr float MAX_VEL = 15.0f;
    constexpr int MIN_SAMPLES = 4, MAX_SAMPLES = 64;
    
    // Filter parameters
    constexpr float Q_KF = 0.5f, R_KF = 200.0f, EPS_KF = 1e-6f;
    constexpr float Q_ANG = 0.0005f, Q_VEL = 0.005f, R_ANG = 0.2f, EPS_EKF = 1e-6f;
    constexpr float OUTLIER_THRESH = 3.14159f;
    
    // Motion detection
    constexpr int VEL_HIST = 15, ANG_HIST = 7;
    constexpr float VEL_STEP = 0.04f, MIN_CONSIST = 0.05f, NOISE_THRESH = 0.015f;
    constexpr unsigned long LOCK_TIMEOUT = 250000UL; // microseconds
    constexpr float REV_THRESH = 2.356f; // 0.75 * PI
    constexpr int DIR_COUNT = 5;
    constexpr float TREND_THRESH = 0.008f;
}

// Lookup table for pulse offset calibration
struct LUT { int pulse; float offset; };
LUT lut[C::CAL_POINTS];

// Kalman filter matrices
#define N3 3
#define N2 2
BLA::Matrix<N3,1> x_kf;
BLA::Matrix<N3,N3> P_kf, Q_kf, R_kf, I3;
BLA::Matrix<N2,1> x_ekf;
BLA::Matrix<N2,N2> P_ekf, Q_ekf, I2, P_init;
BLA::Matrix<1,1> R_ekf;

// History and state
struct Hist { float ang, vel; unsigned long t; };
Hist vel_hist[C::VEL_HIST], ang_hist[C::ANG_HIST];
int vel_idx, ang_idx;
bool vel_full, ang_full, locked;

unsigned long last_nz_vel, last_ekf_us;
long rev_count, mech_steps;
float unwrap_ang, last_valid_ang, smooth_vel, vel_trend, last_stable_vel;
int consist_dir;
bool stopping;

void setup() {
    Serial.begin(115200);
    while (!Serial) delay(10);
    Serial.println(F("Init encoder..."));
    
    // Initialize matrices
    I3.Fill(0); I2.Fill(0);
    for(int i=0; i<N3; i++) I3(i,i) = 1;
    for(int i=0; i<N2; i++) I2(i,i) = 1;
    
    Q_kf.Fill(0); R_kf.Fill(0);
    for(int i=0; i<N3; i++) { Q_kf(i,i) = C::Q_KF; R_kf(i,i) = C::R_KF; }
    
    Q_ekf.Fill(0);
    Q_ekf(0,0) = C::Q_ANG; Q_ekf(1,1) = C::Q_VEL;
    R_ekf(0,0) = C::R_ANG;
    
    // Setup ADC and pins
    analogReference(INTERNAL);
    ADCSRA = (ADCSRA & ~7) | 5; // Set prescaler to 32
    pinMode(C::COMMON_PIN, OUTPUT);
    digitalWrite(C::COMMON_PIN, LOW);
    
    // Calibrate lookup table
    calibrate();
    
    // Initialize state
    x_kf.Fill(0);
    P_kf.Fill(0);
    for(int i=0; i<N3; i++) P_kf(i,i) = 100;
    
    for(int i=0; i<C::VEL_HIST; i++) vel_hist[i] = {0,0,0};
    for(int i=0; i<C::ANG_HIST; i++) ang_hist[i] = {0,0,0};
    
    sense();
    x_ekf(0) = getAngle();
    x_ekf(1) = 0;
    last_valid_ang = x_ekf(0);
    P_ekf.Fill(0);
    P_ekf(0,0) = P_ekf(1,1) = 1;
    P_init = P_ekf;
    
    last_ekf_us = micros();
    unwrap_ang = x_ekf(0);
    
    Serial.println(F("Ready"));
}

void loop() {
    sense();
    float ang = getAngle();
    
    if(!isnan(ang) && !isinf(ang)) {
        last_valid_ang = ang;
        ekf(ang);
    }
    
    mech_steps = round(unwrap_ang / C::TWO_PI * C::ELECTRICAL_STEPS_PER_REV);
    log();
    
    if(Serial.available()) {
        if(toupper(Serial.read()) == 'T') {
            Serial.println(F("Recalibrating..."));
            calibrate();
            // Reset state
            sense();
            x_ekf(0) = getAngle();
            x_ekf(1) = 0;
            last_valid_ang = x_ekf(0);
            P_ekf = P_init;
            last_ekf_us = micros();
            rev_count = 0;
            unwrap_ang = x_ekf(0);
            locked = false;
            vel_idx = ang_idx = 0;
            vel_full = ang_full = false;
            smooth_vel = vel_trend = 0;
            consist_dir = 0;
            stopping = false;
            Serial.println(F("Reset complete"));
        }
    }
}

void calibrate() {
    Serial.println(F("Calibrating LUT..."));
    for(int i=0; i<C::CAL_POINTS; i++) {
        int pulse = C::CAL_PULSE_US[i];
        long sum = 0;
        
        for(int s=0; s<C::CAL_SAMPLES; s++) {
            digitalWrite(C::COMMON_PIN, HIGH);
            delayMicroseconds(pulse);
            digitalWrite(C::COMMON_PIN, LOW);
            delayMicroseconds(C::SETTLE_US);
            sum += analogRead(C::PHASE_PINS[0]);
        }
        
        lut[i] = {pulse, (float)sum / C::CAL_SAMPLES};
    }
    
    Serial.println(F("Pulse\tOffset"));
    for(int i=0; i<C::CAL_POINTS; i++) {
        Serial.print(lut[i].pulse); Serial.print('\t');
        Serial.println(lut[i].offset, 2);
    }
}

float getOffset(int pulse) {
    if(pulse <= lut[0].pulse) return lut[0].offset;
    if(pulse >= lut[C::CAL_POINTS-1].pulse) return lut[C::CAL_POINTS-1].offset;
    
    int i = 0;
    while(i < C::CAL_POINTS-1 && lut[i+1].pulse < pulse) i++;
    
    float t = (float)(pulse - lut[i].pulse) / (lut[i+1].pulse - lut[i].pulse);
    return lut[i].offset + t * (lut[i+1].offset - lut[i].offset);
}

int getPulse() {
    float vel = fabs(x_ekf(1));
    float norm = constrain(vel / C::MAX_VEL, 0.0f, 1.0f);
    return C::MIN_PULSE + norm * (C::MAX_PULSE - C::MIN_PULSE);
}

int getSamples() {
    float vel = fabs(x_ekf(1));
    float norm = constrain(vel / C::MAX_VEL, 0.0f, 1.0f);
    return C::MAX_SAMPLES - norm * (C::MAX_SAMPLES - C::MIN_SAMPLES);
}

void sense() {
    long sums[3] = {0};
    int pulse = getPulse();
    int samples = getSamples();
    float offset = getOffset(pulse);
    
    for(int i=0; i<samples; i++) {
        digitalWrite(C::COMMON_PIN, HIGH);
        delayMicroseconds(pulse);
        digitalWrite(C::COMMON_PIN, LOW);
        delayMicroseconds(C::SETTLE_US);
        
        for(int j=0; j<3; j++) sums[j] += analogRead(C::PHASE_PINS[j]);
    }
    
    float readings[3];
    for(int i=0; i<3; i++) readings[i] = (float)sums[i] / samples - offset;
    
    kf(readings);
}

void kf(const float z[3]) {
    BLA::Matrix<N3,1> y;
    for(int i=0; i<3; i++) y(i) = z[i] - x_kf(i);
    
    BLA::Matrix<N3,N3> S = P_kf + R_kf;
    for(int i=0; i<N3; i++) S(i,i) += C::EPS_KF;
    
    BLA::Matrix<N3,N3> K = P_kf * BLA::Inverse(S);
    x_kf += K * y;
    P_kf = (I3 - K) * P_kf + Q_kf;
}

float getAngle() {
    float a = x_kf(0), b = x_kf(1), c = x_kf(2);
    float alpha = (2.0f/3.0f) * (a - 0.5f * (b + c));
    float beta = (b - c) / C::SQRT3;
    float ang = atan2(beta, alpha);
    return (isnan(ang) || isinf(ang)) ? last_valid_ang : ang;
}

float norm(float a) { return atan2(sin(a), cos(a)); }

void ekf(float meas_ang) {
    unsigned long now = micros();
    float dt = (now - last_ekf_us) * 1e-6f;
    last_ekf_us = now;
    
    // Predict
    float pred_ang = x_ekf(0) + x_ekf(1) * dt;
    BLA::Matrix<N2,N2> A = {{1, dt}, {0, 1}};
    P_ekf = A * P_ekf * BLA::Transpose(A) + Q_ekf;
    
    // Update
    float y = norm(meas_ang - pred_ang);
    if(fabs(y) > C::OUTLIER_THRESH) return;
    
    float S = R_ekf(0,0) + P_ekf(0,0) + C::EPS_EKF;
    float K0 = P_ekf(0,0) / S, K1 = P_ekf(1,0) / S;
    
    x_ekf(0) = norm(pred_ang + K0 * y);
    x_ekf(1) += K1 * y;
    
    P_ekf(0,0) *= (1 - K0);
    P_ekf(0,1) *= (1 - K0);
    P_ekf(1,0) -= K1 * P_ekf(0,0);
    P_ekf(1,1) -= K1 * P_ekf(0,1);
    
    // Update history
    updateVel(x_ekf(1), now);
    updateAng(x_ekf(0), now);
    smooth_vel = getSmooth();
    vel_trend = getTrend();
    stopping = isStopping();
    
    // Locking logic
    bool moving = fabs(x_ekf(1)) > C::VEL_STEP;
    if(moving && isConsistent(x_ekf(1))) {
        float dir = x_ekf(1) > 0 ? 1 : -1;
        float last_dir = last_stable_vel > 0 ? 1 : -1;
        if(fabs(last_stable_vel) > C::VEL_STEP && dir == last_dir) {
            consist_dir++;
        } else consist_dir = 0;
        last_stable_vel = x_ekf(1);
    } else consist_dir = 0;
    
    if(!locked && moving && consist_dir >= C::DIR_COUNT) {
        locked = true;
        Serial.println(F("Locked"));
    } else if(locked) {
        float diff = norm(x_ekf(0) - unwrap_ang);
        unwrap_ang += diff;
        updateRev(diff);
        if(moving) last_nz_vel = now;
        
        if(stopping || (!moving && (now - last_nz_vel > C::LOCK_TIMEOUT))) {
            locked = false;
            Serial.println(F("Unlocked"));
        }
    }
}

void updateRev(float diff) {
    if(diff > C::REV_THRESH) rev_count--;
    else if(diff < -C::REV_THRESH) rev_count++;
}

void updateVel(float v, unsigned long t) {
    vel_hist[vel_idx] = {0, v, t};
    vel_idx = (vel_idx + 1) % C::VEL_HIST;
    if(vel_idx == 0) vel_full = true;
}

void updateAng(float a, unsigned long t) {
    ang_hist[ang_idx] = {a, 0, t};
    ang_idx = (ang_idx + 1) % C::ANG_HIST;
    if(ang_idx == 0) ang_full = true;
}

float getSmooth() {
    if(!vel_full && vel_idx < 3) return x_ekf(1);
    float sum = 0, w_sum = 0;
    int cnt = vel_full ? C::VEL_HIST : vel_idx;
    for(int i=0; i<cnt; i++) {
        int idx = (vel_idx - 1 - i + C::VEL_HIST) % C::VEL_HIST;
        float w = exp(-i * 0.2f);
        sum += vel_hist[idx].vel * w;
        w_sum += w;
    }
    return w_sum > 0 ? sum / w_sum : x_ekf(1);
}

float getTrend() {
    if(!vel_full && vel_idx < 5) return 0;
    int cnt = min(vel_full ? C::VEL_HIST : vel_idx, 5);
    float sx = 0, sy = 0, sxy = 0, sx2 = 0;
    for(int i=0; i<cnt; i++) {
        int idx = (vel_idx - 1 - i + C::VEL_HIST) % C::VEL_HIST;
        float x = i, y = vel_hist[idx].vel;
        sx += x; sy += y; sxy += x*y; sx2 += x*x;
    }
    float denom = cnt * sx2 - sx * sx;
    return fabs(denom) < 1e-6f ? 0 : (cnt * sxy - sx * sy) / denom;
}

bool isConsistent(float v) {
    if(!vel_full && vel_idx < 3) return false;
    int cnt = min(vel_full ? C::VEL_HIST : vel_idx, 5);
    float var = 0;
    for(int i=0; i<cnt; i++) {
        int idx = (vel_idx - 1 - i + C::VEL_HIST) % C::VEL_HIST;
        float diff = vel_hist[idx].vel - v;
        var += diff * diff;
    }
    return sqrt(var / cnt) < C::NOISE_THRESH;
}

bool isStopping() {
    return fabs(smooth_vel) < C::MIN_CONSIST && vel_trend < 0 && 
           fabs(vel_trend) > C::TREND_THRESH && isConsistent(smooth_vel);
}

void log() {
    static unsigned long last = 0;
    if(millis() - last < 100) return;
    last = millis();
    
    Serial.print(F("A:")); Serial.print(degrees(x_ekf(0)), 1);
    Serial.print(F(" V:")); Serial.print(x_ekf(1), 2);
    Serial.print(F(" R:")); Serial.print(rev_count);
    Serial.print(F(" S:")); Serial.print(mech_steps);
    Serial.print(F(" L:")); Serial.println(locked);
}
