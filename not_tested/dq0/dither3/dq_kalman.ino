#include <math.h>

// --- Pin and System Constants ---
const int PHASE_A_PIN = A0;
const int PHASE_B_PIN = A1;
const int PHASE_C_PIN = A2;
const int DITHER_PIN = 9;

const float VREF = 5.0;
const int ADC_MAX = 1023;
const float OFFSET_VOLTAGE = VREF / 2.0;
const float DEG_PER_RAD = 180.0 / PI;
const float SQRT3_DIV_2 = 0.8660254; // Pre-calculate for efficiency

// --- Dither and Oversampling Parameters ---
const float MAX_DITHER_AMP = 12.0;
const float DITHER_FREQ = 100.0;

// --- Kalman Filter Tuning Parameters ---
// Process Noise (Q): How much we trust our model (lower = more trust, smoother)
// Represents the expected variance of change per second.
const float Q_VOLTAGE = 0.01;   // Variance for Vd, Vq
const float Q_OMEGA = 1.0;      // Variance for omega (velocity)

// Measurement Noise (R): How much we trust the ADC measurement (higher = more noise)
const float R_BASE_NOISE = 0.05; // Base ADC noise variance

// --- Data Structures ---
struct Complex {
  float re, im; // Represents alpha, beta
};

struct DQ {
  float d, q;
};

// --- Kalman Filter State Variables ---
// State vector x = [Vd, Vq, omega]'
float x_state[3] = {0, 0, 0}; // Vd, Vq, omega

// State covariance matrix P
float P[3][3] = {
  {1, 0, 0},
  {0, 1, 0},
  {0, 0, 1}
};

// --- Global Kinematic Variables ---
float master_theta = 0; // The master filtered angle from the PLL
float master_omega = 0; // The master filtered velocity from the PLL
unsigned long prev_time = 0;

void setup() {
  Serial.begin(115200);
  pinMode(DITHER_PIN, OUTPUT);
  analogWrite(DITHER_PIN, 127);
  prev_time = micros(); // Use micros for higher precision dt
}

void loop() {
  // --- High-Precision Non-Blocking Timer ---
  unsigned long current_time = micros();
  float dt = (current_time - prev_time) * 1e-6;
  if (dt <= 0.0) return;
  prev_time = current_time;

  // --- 1. Sensor Reading (Adaptive Oversampling) ---
  float v_mag_est = sqrt(x_state[0] * x_state[0] + x_state[1] * x_state[1]);
  int oversample_count = 1;
  if (v_mag_est < 0.05) oversample_count = 8;
  else if (v_mag_est < 0.2) oversample_count = 4;
  else if (v_mag_est < 0.5) oversample_count = 2;
  
  float Va = 0, Vb = 0, Vc = 0;
  for (int i = 0; i < oversample_count; ++i) {
    Va += analogRead(PHASE_A_PIN);
    Vb += analogRead(PHASE_B_PIN);
    Vc += analogRead(PHASE_C_PIN);
  }
  Va = (Va / oversample_count) * VREF / ADC_MAX - OFFSET_VOLTAGE;
  Vb = (Vb / oversample_count) * VREF / ADC_MAX - OFFSET_VOLTAGE;
  Vc = (Vc / oversample_count) * VREF / ADC_MAX - OFFSET_VOLTAGE;

  // --- 2. Coordinate Transforms (abc -> alpha-beta -> dq) ---
  // Clarke Transform (abc -> alpha-beta)
  Complex z_alpha_beta;
  z_alpha_beta.re = (2.0/3.0) * (Va - 0.5 * Vb - 0.5 * Vc);
  z_alpha_beta.im = (2.0/3.0) * (SQRT3_DIV_2 * (Vb - Vc));

  // Park Transform (alpha-beta -> dq)
  // This uses the *previous* step's angle to project the new measurement.
  float cos_theta = cos(master_theta);
  float sin_theta = sin(master_theta);
  DQ z_dq_measured;
  z_dq_measured.d = z_alpha_beta.re * cos_theta + z_alpha_beta.im * sin_theta;
  z_dq_measured.q = -z_alpha_beta.re * sin_theta + z_alpha_beta.im * cos_theta;
  
  // --- 3. Kalman Filter ---
  // A. Prediction Step
  // State prediction is simple (x_hat = F*x), with F=Identity
  // Covariance prediction: P = F*P*F' + Q => P = P + Q*dt
  P[0][0] += Q_VOLTAGE * dt;
  P[1][1] += Q_VOLTAGE * dt;
  P[2][2] += Q_OMEGA * dt;

  // B. Update Step
  // Measurement residual (innovation)
  float y[2];
  y[0] = z_dq_measured.d - x_state[0]; // Vd error
  y[1] = z_dq_measured.q - x_state[1]; // Vq error
  
  // Measurement Noise Covariance (Adaptive)
  float R = R_BASE_NOISE / (1.0 + 50.0 * v_mag_est); // Trust measurements less when signal is small

  // Innovation covariance: S = H*P*H' + R
  // H = [[1,0,0],[0,1,0]], so H*P*H' is the top-left 2x2 of P
  float S[2][2];
  S[0][0] = P[0][0] + R;
  S[0][1] = P[0][1];
  S[1][0] = P[1][0];
  S[1][1] = P[1][1] + R;
  
  // Inverse of S (innovation covariance)
  float S_inv[2][2];
  float det_S = S[0][0] * S[1][1] - S[0][1] * S[1][0];
  det_S = (fabs(det_S) < 1e-9) ? 1e-9 : det_S; // Avoid division by zero
  S_inv[0][0] =  S[1][1] / det_S;
  S_inv[0][1] = -S[0][1] / det_S;
  S_inv[1][0] = -S[1][0] / det_S;
  S_inv[1][1] =  S[0][0] / det_S;

  // Kalman Gain: K = P*H' * S_inv
  float K[3][2];
  for(int i=0; i<3; ++i) {
    K[i][0] = P[i][0] * S_inv[0][0] + P[i][1] * S_inv[1][0];
    K[i][1] = P[i][0] * S_inv[0][1] + P[i][1] * S_inv[1][1];
  }

  // Update state: x = x + K*y
  for(int i=0; i<3; ++i) {
    x_state[i] += K[i][0] * y[0] + K[i][1] * y[1];
  }
  
  // Update covariance: P = (I - K*H)*P
  float I_KH[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
  for(int i=0; i<3; ++i) {
    I_KH[i][0] -= K[i][0];
    I_KH[i][1] -= K[i][1];
  }
  
  float P_new[3][3] = {{0}};
  for(int i=0; i<3; ++i) {
    for(int j=0; j<3; ++j) {
      for(int k=0; k<3; ++k) {
        P_new[i][j] += I_KH[i][k] * P[k][j];
      }
    }
  }
  for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) P[i][j] = P_new[i][j];

  // --- 4. State Integration and PLL Feedback ---
  master_omega = x_state[2];
  master_theta += master_omega * dt;
  
  // Normalize theta to keep it within -PI to PI
  master_theta = fmod(master_theta + PI, 2 * PI) - PI;

  // --- 5. Adaptive Dithering ---
  float dither_amp = MAX_DITHER_AMP / (1.0 + 100.0 * v_mag_est);
  float triangle_wave = 2.0 * fabs(fmod(current_time * 1e-6 * DITHER_FREQ, 1.0) - 0.5);
  int pwm_value = 127 + int((triangle_wave - 0.5) * 2.0 * dither_amp);
  analogWrite(DITHER_PIN, pwm_value);

  // --- Output ---
  Serial.print("θ(deg):"); Serial.print(master_theta * DEG_PER_RAD, 2);
  Serial.print("\tω(rad/s):"); Serial.print(master_omega, 2);
  Serial.print("\tVd:"); Serial.print(x_state[0], 2);
  Serial.print("\tVq:"); Serial.print(x_state[1], 2);
  Serial.print("\tOS:"); Serial.println(oversample_count);
}
