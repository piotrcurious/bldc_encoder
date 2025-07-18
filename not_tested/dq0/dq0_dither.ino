const int PHASE_A_PIN = A0;
const int PHASE_B_PIN = A1;
const int PHASE_C_PIN = A2;
const int DITHER_PIN = 9;

const float VREF = 5.0;
const int ADC_MAX = 1023;
const float DEG_PER_RAD = 180.0 / PI;

const float LAMBDA = 0.15;        // Filter coefficient
const float DT = 0.01;            // Time step (s)

struct Complex {
  float re, im;

  Complex operator+(const Complex& b) const { return {re + b.re, im + b.im}; }
  Complex operator*(float s) const { return {re * s, im * s}; }
  float magnitude() const { return sqrt(re * re + im * im); }
  float angle() const { return atan2(im, re); }
};

Complex z_curr = {0, 0};    // Current stator vector
Complex momentum = {0, 0};  // Filtered complex momentum

float prev_theta = 0;
float velocity = 0;

unsigned long t_start;

void setup() {
  Serial.begin(115200);
  pinMode(DITHER_PIN, OUTPUT);
  analogWrite(DITHER_PIN, 127);  // Start with 50% duty cycle (VREF/2)

  t_start = millis();
}

void loop() {
  float time = (millis() - t_start) / 1000.0;

  // --- Dither Injection ---
  // Triangle wave dither oscillating around 127 ± 8 (approx. ±0.04V after filter)
  // 100 Hz triangular wave mapped to 8-bit PWM value
  const int dither_amplitude = 8;
  const float dither_freq = 100.0; // Hz

  float triangle = 2.0 * fabs(2.0 * (time * dither_freq - floor(time * dither_freq + 0.5))) - 1.0;
  int pwm_value = 127 + int(triangle * dither_amplitude);
  analogWrite(DITHER_PIN, pwm_value);

  // --- ADC and Signal Conditioning ---
  float OFFSET_VOLTAGE = VREF / 2.0;

  float Va = analogRead(PHASE_A_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;
  float Vb = analogRead(PHASE_B_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;
  float Vc = analogRead(PHASE_C_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;

  // Clarke Transform (abc → αβ)
  float Valpha = (2.0 / 3.0) * (Va - 0.5 * Vb - 0.5 * Vc);
  float Vbeta  = (2.0 / 3.0) * ((sqrt(3.0) / 2.0) * (Vb - Vc));
  z_curr = {Valpha, Vbeta};

  // Complex Momentum Filtering
  momentum = momentum * (1.0 - LAMBDA) + z_curr * LAMBDA;

  float theta = momentum.angle();
  if (theta < 0) theta += 2 * PI;

  // Phase unwrapping for velocity
  float dtheta = theta - prev_theta;
  if (dtheta > PI) dtheta -= 2 * PI;
  if (dtheta < -PI) dtheta += 2 * PI;

  velocity = dtheta / DT;
  prev_theta = theta;

  // Optional: correlation
  float z_mag = z_curr.magnitude();
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
  Serial.print("\tPWM: "); Serial.println(pwm_value);

  delay(int(DT * 1000));
} 
