 const int PHASE_A_PIN = A0;
const int PHASE_B_PIN = A1;
const int PHASE_C_PIN = A2;

const float VREF = 5.0;
const int ADC_MAX = 1023;
const float OFFSET_VOLTAGE = VREF / 2.0;

const float DEG_PER_RAD = 180.0 / PI;
const float LAMBDA = 0.15;        // Momentum filter coefficient
const float DT = 0.01;            // Time step in seconds (10ms)

struct Complex {
  float re;
  float im;

  Complex operator+(const Complex& b) {
    return {re + b.re, im + b.im};
  }

  Complex operator*(float s) {
    return {re * s, im * s};
  }

  float magnitude() {
    return sqrt(re * re + im * im);
  }

  float angle() {
    return atan2(im, re);
  }
};

Complex z_curr = {0, 0};     // Current stator vector
Complex momentum = {0, 0};   // Filtered momentum

float prev_theta = 0;
float velocity = 0;

void setup() {
  Serial.begin(115200);
}

void loop() {
  // Read voltages and remove common-mode offset
  float Va = analogRead(PHASE_A_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;
  float Vb = analogRead(PHASE_B_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;
  float Vc = analogRead(PHASE_C_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;

  // Clarke Transform (abc → αβ)
  float Valpha = (2.0 / 3.0) * (Va - 0.5 * Vb - 0.5 * Vc);
  float Vbeta  = (2.0 / 3.0) * ((sqrt(3.0) / 2.0) * (Vb - Vc));

  // Current complex vector (z_n)
  z_curr = {Valpha, Vbeta};

  // Filter: momentum = (1 - λ) * momentum + λ * z_curr
  momentum = momentum * (1.0 - LAMBDA) + z_curr * LAMBDA;

  // Angle estimation
  float theta = atan2(momentum.im, momentum.re);
  if (theta < 0) theta += 2 * PI;

  // Velocity estimation via phase difference
  float dtheta = theta - prev_theta;
  if (dtheta > PI) dtheta -= 2 * PI;
  if (dtheta < -PI) dtheta += 2 * PI;

  velocity = dtheta / DT;
  prev_theta = theta;

  // Optional: magnitude correlation (normalized dot product)
  float z_mag = z_curr.magnitude();
  float m_mag = momentum.magnitude();
  float corr = 0;
  if (z_mag > 0 && m_mag > 0) {
    corr = (z_curr.re * momentum.re + z_curr.im * momentum.im) / (z_mag * m_mag);
  }

  // Output
  Serial.print("θ (deg): "); Serial.print(theta * DEG_PER_RAD, 2);
  Serial.print("\tω (rad/s): "); Serial.print(velocity, 2);
  Serial.print("\t|z|: "); Serial.print(z_mag, 3);
  Serial.print("\tCorr: "); Serial.println(corr, 3);

  delay(int(DT * 1000));
}
