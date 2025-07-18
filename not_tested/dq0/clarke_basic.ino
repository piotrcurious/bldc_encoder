 // Pin configuration
const int PHASE_A_PIN = A0;
const int PHASE_B_PIN = A1;
const int PHASE_C_PIN = A2;

// ADC characteristics
const float VREF = 5.0;          // Supply/reference voltage
const int ADC_MAX = 1023;        // 10-bit ADC

// Center offset (e.g., voltage divider gives mid-point at VREF / 2)
const float OFFSET_VOLTAGE = VREF / 2.0;

// Helper constants
const float TWO_THIRD = 2.0 / 3.0;
const float ONE_OVER_SQRT3 = 1.0 / sqrt(3.0);
const float DEG_PER_RAD = 180.0 / PI;

void setup() {
  Serial.begin(115200);
}

void loop() {
  // Read and convert ADC values to voltage (assume common-mode already handled)
  float Va = analogRead(PHASE_A_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;
  float Vb = analogRead(PHASE_B_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;
  float Vc = analogRead(PHASE_C_PIN) * VREF / ADC_MAX - OFFSET_VOLTAGE;

  // Clarke Transformation (to αβ frame)
  float Valpha = Va;
  float Vbeta = (Vb - Vc) * ONE_OVER_SQRT3;

  // Rotor electrical angle (atan2 gives angle in radians)
  float theta_rad = atan2(Vbeta, Valpha);
  if (theta_rad < 0) theta_rad += 2 * PI; // Normalize to [0, 2π)

  // Optionally convert to degrees
  float theta_deg = theta_rad * DEG_PER_RAD;

  // Output
  Serial.print("Electrical Angle (rad): ");
  Serial.print(theta_rad, 4);
  Serial.print("\t(deg): ");
  Serial.println(theta_deg, 2);

  delay(10);
}
