 const int PHASE_A_PIN = A0;
const int PHASE_B_PIN = A1;
const int PHASE_C_PIN = A2;

const float VREF = 5.0;           // Reference voltage
const int ADC_MAX = 1023;         // 10-bit ADC
const float OFFSET_VOLTAGE = VREF / 2.0;

const float DEG_PER_RAD = 180.0 / PI;

void setup() {
  Serial.begin(115200);
}

void loop() {
  // Step 1: Read and normalize analog inputs
  float Va = analogRead(PHASE_A_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;
  float Vb = analogRead(PHASE_B_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;
  float Vc = analogRead(PHASE_C_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;

  // Step 2: Clarke Transformation (abc → αβ0)
  float Valpha = (2.0 / 3.0) * (Va - 0.5 * Vb - 0.5 * Vc);
  float Vbeta  = (2.0 / 3.0) * ((sqrt(3.0) / 2.0) * (Vb - Vc));
  float Vzero  = (1.0 / 3.0) * (Va + Vb + Vc);

  // Step 3: Estimate electrical angle θ from alpha-beta (e.g., sensorless, observer, or simplified from atan2)
  float theta = atan2(Vbeta, Valpha);  // In radians
  if (theta < 0) theta += 2 * PI;

  // Step 4: Park Transformation (αβ0 → DQ0)
  float cosTheta = cos(theta);
  float sinTheta = sin(theta);

  float Vd =  Valpha * cosTheta + Vbeta * sinTheta;
  float Vq = -Valpha * sinTheta + Vbeta * cosTheta;
  float V0 = Vzero;

  // Serial Output
  Serial.print("Theta (deg): "); Serial.print(theta * DEG_PER_RAD, 2);
  Serial.print("\tVd: "); Serial.print(Vd, 3);
  Serial.print("\tVq: "); Serial.print(Vq, 3);
  Serial.print("\tV0: "); Serial.println(V0, 3);

  delay(10);
}
