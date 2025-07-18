 /*
 * BLDC Rotor Position Estimation using DQ0 Transformation
 *
 * This sketch estimates the electrical angle of a BLDC motor's rotor
 * by reading the analog voltages from its three phases (A, B, C).
 *
 * It assumes the analog signals have been properly conditioned (amplified/filtered)
 * and are centered around a DC offset. For example, if the signal swings
 * between 0V and 5V, the offset would be at 2.5V.
 *
 * The core logic involves two steps:
 * 1. Clarke Transformation: Converts the three-phase system (abc) into a
 * two-phase orthogonal stationary system (alpha-beta).
 * 2. Park Transformation (implied): The angle of the resulting alpha-beta
 * vector is calculated. This angle IS the rotor's electrical position (theta),
 * which is the input required for a full Park transformation.
 *
 * Hardware Setup:
 * - Connect the conditioned analog output for Phase A to Arduino pin A0.
 * - Connect the conditioned analog output for Phase B to Arduino pin A1.
 * - Connect the conditioned analog output for Phase C to Arduino pin A2.
 * - Ensure the Arduino and the signal conditioning circuit share a common ground.
 */

// --- Pin Definitions ---
const int PHASE_A_PIN = A0; // Analog input for Phase A
const int PHASE_B_PIN = A1; // Analog input for Phase B
const int PHASE_C_PIN = A2; // Analog input for Phase C

// --- Constants ---
// The ADC value corresponding to the center voltage (DC offset).
// For a 10-bit ADC (0-1023) with a 5V reference and a 2.5V offset,
// the value is 1023 * (2.5 / 5.0) = 511.5, which we round to 512.
// **ADJUST THIS VALUE** to match your specific hardware's offset.
const int ADC_OFFSET = 512;

// Pre-calculated mathematical constants for efficiency
const float SQRT_3 = 1.73205081;
const float ONE_OVER_SQRT_3 = 1.0 / SQRT_3;
const float TWO_THIRDS = 2.0 / 3.0;
const float ONE_THIRD = 1.0 / 3.0;
const float RAD_TO_DEG = 180.0 / PI;

void setup() {
  // Initialize Serial communication for debugging and outputting the angle.
  // A baud rate of 115200 is fast and reliable.
  Serial.begin(115200);
  Serial.println("BLDC Rotor Position Estimation Initialized...");

  // pinMode for analog inputs is not strictly required on most Arduino boards
  // as they default to inputs, but it's good practice for clarity.
  pinMode(PHASE_A_PIN, INPUT);
  pinMode(PHASE_B_PIN, INPUT);
  pinMode(PHASE_C_PIN, INPUT);
}

void loop() {
  // 1. Read Raw ADC Values
  // Read the instantaneous voltage from each of the three phases.
  int rawA = analogRead(PHASE_A_PIN);
  int rawB = analogRead(PHASE_B_PIN);
  int rawC = analogRead(PHASE_C_PIN);

  // 2. Center the Signals
  // Subtract the DC offset to get signed values representing the AC component.
  // This is crucial for the transformations to work correctly.
  float phaseA = (float)(rawA - ADC_OFFSET);
  float phaseB = (float)(rawB - ADC_OFFSET);
  float phaseC = (float)(rawC - ADC_OFFSET);

  // 3. Clarke Transformation (abc -> alpha-beta)
  // This transforms the three-phase stationary reference frame into a
  // two-phase orthogonal stationary reference frame (alpha, beta).
  // We use the power-invariant form of the transformation.
  float v_alpha = TWO_THIRDS * phaseA - ONE_THIRD * phaseB - ONE_THIRD * phaseC;
  float v_beta = ONE_OVER_SQRT_3 * (phaseB - phaseC);

  // 4. Calculate Rotor Angle (theta)
  // The angle of the vector formed by v_alpha and v_beta in the complex plane
  // corresponds to the electrical angle of the rotor.
  // atan2 is used because it correctly handles all four quadrants and avoids
  // division by zero, giving a result from -PI to +PI.
  float rotor_angle_rad = atan2(v_beta, v_alpha);

  // 5. Convert and Normalize Angle
  // Convert the angle from radians to degrees for easier interpretation.
  float rotor_angle_deg = rotor_angle_rad * RAD_TO_DEG;

  // Optional: Map the angle from [-180, 180] to [0, 360] degrees.
  if (rotor_angle_deg < 0) {
    rotor_angle_deg += 360.0;
  }

  // 6. Output the Result
  // Print the calculated electrical angle to the Serial Monitor.
  // You can use the Arduino IDE's Serial Plotter to visualize the angle changing.
  Serial.print("Alpha: ");
  Serial.print(v_alpha);
  Serial.print(" | Beta: ");
  Serial.print(v_beta);
  Serial.print(" | Angle (deg): ");
  Serial.println(rotor_angle_deg);

  // A small delay to prevent flooding the serial port and to allow the ADC
  // to settle. Adjust as needed for your application's speed requirements.
  delay(2);
}
