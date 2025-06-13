#include <Arduino.h>
#include <BasicLinearAlgebra.h>

namespace Config {
    constexpr int PHASE_A_PIN = A0, PHASE_B_PIN = A1, PHASE_C_PIN = A2, COMMON_PIN = 9;
    constexpr float ELECTRICAL_STEPS_PER_REV = 6.0f, TWO_PI_F = 2.0f * PI, SQRT_3 = 1.7320508f;
    constexpr int NUM_SAMPLES = 16;
    constexpr float ANGLE_UPDATE_DEADBAND = 0.02f, VELOCITY_STEP_THRESHOLD = 0.2f;
    constexpr unsigned long ZERO_VEL_LOCK_TIMEOUT = 500000UL;
    int EXCITATION_PWM_VALUE = 1;
    constexpr int EXCITATION_PULSE_WIDTH_US = 1000, FLOATING_SETTLING_TIME_US = 10;
    constexpr int TARGET_QUIESCENT_AVG_ADC = 80, PWM_ADJUST_STEP = 1, PWM_CALIBRATION_ITERATIONS = 2000, PWM_CALIBRATION_TOLERANCE = 1;
    constexpr float Q_DIAGONAL_VALUE_KF1 = 1.0f, R_DIAGONAL_VALUE_KF1 = 100.0f;
    constexpr float Q_ANGLE_EKF = 0.001f, Q_VEL_EKF = 0.1f, R_ANGLE_EKF = 0.5f;
}

#define N_STATES_KF1 3
BLA::Matrix<N_STATES_KF1, 1> x_hat_kf1;
BLA::Matrix<N_STATES_KF1, N_STATES_KF1> P_kf1, Q_kf1, R_kf1, I_matrix_kf1;

#define N_STATES_EKF 2
BLA::Matrix<N_STATES_EKF, 1> x_ekf;
BLA::Matrix<N_STATES_EKF, N_STATES_EKF> P_ekf, Q_ekf, I_matrix_ekf;
BLA::Matrix<1, 1> R_ekf;

bool ekfLocked = false;
unsigned long lastNonZeroVelMicros = 0, prev_ekf_micros = 0;
long revolution_count = 0, estimated_mechanical_steps = 0;
float unwrapped_electrical_angle = 0.0f;

void initializeMatrices(), autoTuneExcitationPwm(), executeSensingCycle(), applyMultivariateKalmanFilter(const int[]);
float getInferredElectricalAngle(), normalizeAngle(float);
void applyEKF(float), logData();
int checkedAnalogRead(int);

void setup() {
    Serial.begin(115200); while (!Serial) delay(10);
    Serial.println(F("Initializing BLDC Sensorless Encoder..."));
    initializeMatrices();
    analogReference(INTERNAL);
    ADCSRA = (ADCSRA & ~0x07) | 0x05;
    TCCR1B = (TCCR1B & 0b11111000) | 0x01;
    pinMode(Config::COMMON_PIN, OUTPUT); digitalWrite(Config::COMMON_PIN, LOW);
    autoTuneExcitationPwm();

    x_hat_kf1.Fill((float)Config::TARGET_QUIESCENT_AVG_ADC);
    P_kf1.Fill(0.0f); for (int i = 0; i < N_STATES_KF1; i++) P_kf1(i, i) = 100.0f;

    executeSensingCycle();
    x_ekf(0) = getInferredElectricalAngle(); x_ekf(1) = 0.0f;
    P_ekf.Fill(0.0f); P_ekf(0, 0) = 1.0f; P_ekf(1, 1) = 1.0f;
    prev_ekf_micros = micros(); unwrapped_electrical_angle = x_ekf(0);
    Serial.println(F("System ready. Rotate motor to begin estimation."));
}

void loop() {
    executeSensingCycle();
    float inferred_angle = getInferredElectricalAngle();
    applyEKF(inferred_angle);
    estimated_mechanical_steps = round(unwrapped_electrical_angle / Config::TWO_PI_F * Config::ELECTRICAL_STEPS_PER_REV);
    logData();

    if (Serial.available() > 0) {
        char command = toupper(Serial.read());
        while (Serial.available() > 0) Serial.read();
        if (command == 'T') {
            Serial.println(F("\nRe-tuning requested..."));
            autoTuneExcitationPwm();
            executeSensingCycle();
            x_ekf(0) = getInferredElectricalAngle(); x_ekf(1) = 0.0f;
            P_ekf(0, 0) = 1.0f; P_ekf(1, 1) = 1.0f;
            prev_ekf_micros = micros(); revolution_count = 0;
            unwrapped_electrical_angle = x_ekf(0); ekfLocked = false;
            Serial.println(F("Re-tuning complete. System ready."));
        }
    }
}

void initializeMatrices() {
    I_matrix_kf1.Fill(0.0f); I_matrix_ekf.Fill(0.0f);
    for (int i = 0; i < N_STATES_KF1; i++) I_matrix_kf1(i, i) = 1.0f;
    for (int i = 0; i < N_STATES_EKF; i++) I_matrix_ekf(i, i) = 1.0f;

    Q_kf1.Fill(0.0f); R_kf1.Fill(0.0f);
    for (int i = 0; i < N_STATES_KF1; i++) {
        Q_kf1(i, i) = Config::Q_DIAGONAL_VALUE_KF1;
        R_kf1(i, i) = Config::R_DIAGONAL_VALUE_KF1;
    }

    Q_ekf.Fill(0.0f); Q_ekf(0, 0) = Config::Q_ANGLE_EKF; Q_ekf(1, 1) = Config::Q_VEL_EKF;
    R_ekf(0, 0) = Config::R_ANGLE_EKF;
} 

void autoTuneExcitationPwm() {
    int sum = 0, prev_avg = 0;
    for (int iter = 0; iter < Config::PWM_CALIBRATION_ITERATIONS; iter++) {
        sum = 0;
        for (int i = 0; i < Config::NUM_SAMPLES; i++) {
            digitalWrite(Config::COMMON_PIN, LOW);
            analogWrite(Config::COMMON_PIN, Config::EXCITATION_PWM_VALUE);
            delayMicroseconds(Config::EXCITATION_PULSE_WIDTH_US);
            pinMode(Config::COMMON_PIN, INPUT);
            delayMicroseconds(Config::FLOATING_SETTLING_TIME_US);
            int reading = analogRead(Config::PHASE_A_PIN);
            sum += reading;
        }
        int avg = sum / Config::NUM_SAMPLES;
        if (abs(avg - Config::TARGET_QUIESCENT_AVG_ADC) < Config::PWM_CALIBRATION_TOLERANCE) break;
        Config::EXCITATION_PWM_VALUE += (avg < Config::TARGET_QUIESCENT_AVG_ADC) ? Config::PWM_ADJUST_STEP : -Config::PWM_ADJUST_STEP;
        Config::EXCITATION_PWM_VALUE = constrain(Config::EXCITATION_PWM_VALUE, 0, 255);
        prev_avg = avg;
    }
    Serial.print(F("Calibrated PWM value: "));
    Serial.println(Config::EXCITATION_PWM_VALUE);
}

void executeSensingCycle() {
    int readings[3] = {0};
    for (int i = 0; i < Config::NUM_SAMPLES; i++) {
        digitalWrite(Config::COMMON_PIN, LOW);
        analogWrite(Config::COMMON_PIN, Config::EXCITATION_PWM_VALUE);
        delayMicroseconds(Config::EXCITATION_PULSE_WIDTH_US);
        pinMode(Config::COMMON_PIN, INPUT);
        delayMicroseconds(Config::FLOATING_SETTLING_TIME_US);
        readings[0] += checkedAnalogRead(Config::PHASE_A_PIN);
        readings[1] += checkedAnalogRead(Config::PHASE_B_PIN);
        readings[2] += checkedAnalogRead(Config::PHASE_C_PIN);
    }
    for (int i = 0; i < 3; i++) readings[i] /= Config::NUM_SAMPLES;
    applyMultivariateKalmanFilter(readings);
}

void applyMultivariateKalmanFilter(const int readings[3]) {
    BLA::Matrix<N_STATES_KF1, 1> z, y, S;
    for (int i = 0; i < 3; i++) z(i) = readings[i];
    y = z - x_hat_kf1;
    S = P_kf1 + R_kf1;
    auto K = P_kf1 * BLA::Inverse(S);
    x_hat_kf1 += K * y;
    P_kf1 = (I_matrix_kf1 - K) * P_kf1 + Q_kf1;
}

float getInferredElectricalAngle() {
    float a = x_hat_kf1(0), b = x_hat_kf1(1), c = x_hat_kf1(2);
    float alpha = (2.0f / 3.0f) * (a - 0.5f * (b + c));
    float beta = (2.0f / 3.0f) * ((b - c) / Config::SQRT_3);
    return atan2f(beta, alpha);
}

void applyEKF(float measured_angle) {
    unsigned long now = micros();
    float dt = (now - prev_ekf_micros) / 1e6f;
    prev_ekf_micros = now;

    float angle = x_ekf(0), velocity = x_ekf(1);
    float predicted_angle = normalizeAngle(angle + velocity * dt);
    float predicted_velocity = velocity;

    BLA::Matrix<N_STATES_EKF, 1> x_pred = {predicted_angle, predicted_velocity};
    BLA::Matrix<N_STATES_EKF, N_STATES_EKF> A = {{1, dt}, {0, 1}};
    P_ekf = A * P_ekf * BLA::Transpose(A) + Q_ekf;

    float y = normalizeAngle(measured_angle - predicted_angle);
    BLA::Matrix<1, 1> S = R_ekf + P_ekf(0, 0);
    BLA::Matrix<N_STATES_EKF, 1> K = {P_ekf(0, 0) / S(0, 0), P_ekf(1, 0) / S(0, 0)};
    x_ekf = x_pred + K * y;

    BLA::Matrix<N_STATES_EKF, N_STATES_EKF> I_KH;
    I_KH = I_matrix_ekf;
    I_KH(0, 0) -= K(0, 0);
    I_KH(1, 0) -= K(1, 0);
    P_ekf = I_KH * P_ekf;

    float vel = x_ekf(1);
    if (fabsf(vel) > Config::VELOCITY_STEP_THRESHOLD) {
        if (!ekfLocked) {
            ekfLocked = true;
            Serial.println(F("Angle locked."));
        }
        lastNonZeroVelMicros = now;
        float angle_diff = normalizeAngle(x_ekf(0) - unwrapped_electrical_angle);
        unwrapped_electrical_angle += angle_diff;
        if (angle_diff > PI) revolution_count--;
        else if (angle_diff < -PI) revolution_count++;
    } else if (ekfLocked && (now - lastNonZeroVelMicros > Config::ZERO_VEL_LOCK_TIMEOUT)) {
        ekfLocked = false;
        Serial.println(F("Angle unlocked due to zero velocity."));
    }
}

float normalizeAngle(float a) {
    while (a < -PI) a += Config::TWO_PI_F;
    while (a >  PI) a -= Config::TWO_PI_F;
    return a;
}

int checkedAnalogRead(int pin) {
    int val = analogRead(pin);
    if (val < 0) val = 0;
    if (val > 1023) val = 1023;
    return val;
}

void logData() {
    Serial.print("Enc angle: "); Serial.print(x_ekf(0));
    Serial.print("\tUnwrapped: "); Serial.print(unwrapped_electrical_angle);
    Serial.print("\tSteps: "); Serial.print(estimated_mechanical_steps);
    Serial.print("\tVel: "); Serial.println(x_ekf(1));
}
