 const int PHASE_A_PIN = A0;
const int PHASE_B_PIN = A1;
const int PHASE_C_PIN = A2;
const int DITHER_PIN = 9;

const float VREF = 5.0;
const int ADC_MAX = 1023;
const float OFFSET_VOLTAGE = VREF / 2.0;
const float DEG_PER_RAD = 180.0 / PI;
const float DT = 0.01; // seconds

const float LAMBDA = 0.15; // momentum smoothing
const float MAX_DITHER_AMP = 12; // max PWM units (8-bit)
const float DITHER_FREQ = 100.0; // Hz
const int MAX_OVERSAMPLE = 8;

struct Complex {
  float re, im;

  Complex operator+(const Complex& b) const { return {re + b.re, im + b.im}; }
  Complex operator*(float s) const { return {re * s, im * s}; }
  float magnitude() const { return sqrt(re * re + im * im); }
  float angle() const { return atan2(im, re); }
};

Complex z_curr = {0, 0};
Complex momentum = {0, 0};

float prev_theta = 0;
float velocity = 0;

unsigned long t_start;

void setup() {
  Serial.begin(115200);
  pinMode(DITHER_PIN, OUTPUT);
  analogWrite(DITHER_PIN, 127); // midpoint
  t_start = millis();
}

void loop() {
  float time = (millis() - t_start) / 1000.0;

  // --- Adaptive Oversampling Setup ---
  int oversample_count = 1;
  float Va = 0, Vb = 0, Vc = 0;

  // Initial rough sample to estimate magnitude
  float V1 = analogRead(PHASE_A_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;
  float V2 = analogRead(PHASE_B_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;
  float V3 = analogRead(PHASE_C_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;

  float Valpha_est = (2.0 / 3.0) * (V1 - 0.5 * V2 - 0.5 * V3);
  float Vbeta_est = (2.0 / 3.0) * ((sqrt(3.0) / 2.0) * (V2 - V3));
  float z_mag = sqrt(Valpha_est * Valpha_est + Vbeta_est * Vbeta_est);

  // Oversampling factor (inversely proportional to signal magnitude)
  if (z_mag < 0.05)      oversample_count = 8;
  else if (z_mag < 0.1)  oversample_count = 6;
  else if (z_mag < 0.2)  oversample_count = 4;
  else if (z_mag < 0.5)  oversample_count = 2;

  // --- Oversample Analog Inputs ---
  for (int i = 0; i < oversample_count; ++i) {
    Va += analogRead(PHASE_A_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;
    Vb += analogRead(PHASE_B_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;
    Vc += analogRead(PHASE_C_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;
  }

  Va /= oversample_count;
  Vb /= oversample_count;
  Vc /= oversample_count;

  // --- Clarke Transform ---
  float Valpha = (2.0 / 3.0) * (Va - 0.5 * Vb - 0.5 * Vc);
  float Vbeta  = (2.0 / 3.0) * ((sqrt(3.0) / 2.0) * (Vb - Vc));
  z_curr = {Valpha, Vbeta};
  z_mag = z_curr.magnitude();

  // --- Adaptive Dither Amplitude ---
  float dither_amp = MAX_DITHER_AMP / (1.0 + 20.0 * z_mag); // smaller amplitude at high signal
  float triangle = 2.0 * fabs(2.0 * (time * DITHER_FREQ - floor(time * DITHER_FREQ + 0.5))) - 1.0;
  int pwm_value = 127 + int(triangle * dither_amp);
  analogWrite(DITHER_PIN, pwm_value);

  // --- Complex Filtered Momentum ---
  momentum = momentum * (1.0 - LAMBDA) + z_curr * LAMBDA;

  float theta = momentum.angle();
  if (theta < 0) theta += 2 * PI;

  // --- Velocity Estimate ---
  float dtheta = theta - prev_theta;
  if (dtheta > PI) dtheta -= 2 * PI;
  if (dtheta < -PI) dtheta += 2 * PI;

  velocity = dtheta / DT;
  prev_theta = theta;

  // --- Correlation (cosine similarity) ---
  float m_mag = momentum.magnitude();
  float corr = 0;
  if (z_mag > 0 && m_mag > 0) {
    corr = (z_curr.re * momentum.re + z_curr.im * momentum.im) / (z_mag * m_mag);
  }

  // --- Output ---
  Serial.print("θ (deg): "); Serial.print(theta * DEG_PER_RAD, 2);
  Serial.print("\tω (rad/s): "); Serial.print(velocity, 2);
  Serial.print("\t|z|: "); Serial.print(z_mag, 3);
  Serial.print("\tCorr: "); Serial.print(corr, 3);
  Serial.print("\tDither: "); Serial.print(dither_amp, 2);
  Serial.print("\tOS: "); Serial.println(oversample_count);

  delay(int(DT * 1000));
}
