/*
 * Novel Mini BLDC 3-Phase Controller for Arduino Uno
 *
 * Description:
 * This code implements a sensorless BLDC motor controller using a unique
 * push-pull topology. It's designed for small motors where the MCU's
 * internal clamping diodes can handle the inductive spikes.
 *
 * Connections:
 * - PWM Pin 9: Connect to the common (center-tap) of the BLDC motor.
 * - Pin A0 (14): Phase A
 * - Pin A1 (15): Phase B
 * - Pin A2 (16): Phase C
 *
 * How it Works:
 * 1. Push-Pull Driver: The PWM pin acts as the high-side or low-side driver.
 * - When PWM is HIGH, an active phase is set to LOW, current flows PWM -> Phase.
 * - When PWM is LOW, an active phase is set to HIGH, current flows Phase -> PWM.
 * 2. Sensorless BEMF Detection: The position of the rotor is detected by measuring
 * the Back EMF on the currently unused (floating) phase.
 * 3. PWM Synchronization: BEMF is read during the PWM-off cycle to avoid noise.
 * This is achieved by synchronizing with the Timer1 overflow interrupt.
 * 4. 6-Step Commutation: The motor is driven using a standard 6-step sequence.
 * The timing of the steps is determined by BEMF zero-crossings.
 */

// --- Pin Definitions ---
const int PHASE_A_PIN = A0;
const int PHASE_B_PIN = A1;
const int PHASE_C_PIN = A2;
const int PWM_COMMON_PIN = 9; // Must be a Timer1 PWM pin (9 or 10)

// --- Motor Control Parameters ---
volatile int commutation_step = 0;
volatile bool commutate_now = false;
volatile unsigned long last_commutation_time = 0;
volatile unsigned long current_commutation_delay = 5000; // Initial delay in microseconds

// PWM duty cycle (0-255). Start slow.
// This value sets the ICR1 register for Timer1, effectively setting the PWM frequency.
// A lower value gives a higher frequency. Let's set it for ~16 kHz.
// F_CPU / (2 * N * TOP) => 16,000,000 / (2 * 1 * 500) = 16 kHz
const int PWM_TOP_VALUE = 500;
// This controls the duty cycle (speed). Value should be less than PWM_TOP_VALUE.
volatile int motor_speed_pwm = 120; // Start at a low speed

// --- BEMF Detection Parameters ---
const int BEMF_THRESHOLD_HIGH = 400; // ADC value threshold for rising BEMF
const int BEMF_THRESHOLD_LOW = 200;  // ADC value threshold for falling BEMF
volatile bool bemf_high = false;     // Flag to track BEMF state

// --- Startup Parameters ---
const int STARTUP_STEPS = 200; // Number of forced steps to start the motor
int startup_step_counter = 0;

void setup() {
  Serial.begin(9600);
  Serial.println("BLDC Controller Initializing...");

  // --- Initialize Timer1 for Fast PWM on Pin 9 ---
  // We will control the timer registers directly for precise control.
  pinMode(PWM_COMMON_PIN, OUTPUT);

  // TCCR1A: Clear on compare match for OC1A (Pin 9), Fast PWM mode with TOP=ICR1
  TCCR1A = (1 << COM1A1) | (1 << WGM11);
  // TCCR1B: Set prescaler to 1 (no prescaling), Fast PWM mode
  TCCR1B = (1 << WGM13) | (1 << WGM12) | (1 << CS10);

  // Set the PWM frequency.
  ICR1 = PWM_TOP_VALUE;

  // Set the initial motor speed (duty cycle).
  OCR1A = motor_speed_pwm;

  // Enable the Timer1 Overflow Interrupt. This is our synchronization point.
  TIMSK1 |= (1 << TOIE1);

  // Enable global interrupts
  sei();

  Serial.println("Controller Ready. Starting motor...");
}

// --- Commutation Logic ---
// This function sets the pin modes and states according to the 6-step sequence.
void commutate() {
  // The logic is based on which phase is driven HIGH, which is driven LOW,
  // and which is left floating (INPUT) to read BEMF.

  switch (commutation_step) {
    case 0: // A->B, C is floating
      // Set phase C to input to read BEMF
      pinMode(PHASE_C_PIN, INPUT);
      // Set phase A to be driven when PWM is low
      digitalWrite(PHASE_A_PIN, HIGH);
      pinMode(PHASE_A_PIN, OUTPUT);
      // Set phase B to be driven when PWM is high
      digitalWrite(PHASE_B_PIN, LOW);
      pinMode(PHASE_B_PIN, OUTPUT);
      break;

    case 1: // C->B, A is floating
      pinMode(PHASE_A_PIN, INPUT);
      digitalWrite(PHASE_C_PIN, HIGH);
      pinMode(PHASE_C_PIN, OUTPUT);
      digitalWrite(PHASE_B_PIN, LOW);
      pinMode(PHASE_B_PIN, OUTPUT);
      break;

    case 2: // C->A, B is floating
      pinMode(PHASE_B_PIN, INPUT);
      digitalWrite(PHASE_C_PIN, HIGH);
      pinMode(PHASE_C_PIN, OUTPUT);
      digitalWrite(PHASE_A_PIN, LOW);
      pinMode(PHASE_A_PIN, OUTPUT);
      break;

    case 3: // B->A, C is floating
      pinMode(PHASE_C_PIN, INPUT);
      digitalWrite(PHASE_B_PIN, HIGH);
      pinMode(PHASE_B_PIN, OUTPUT);
      digitalWrite(PHASE_A_PIN, LOW);
      pinMode(PHASE_A_PIN, OUTPUT);
      break;

    case 4: // B->C, A is floating
      pinMode(PHASE_A_PIN, INPUT);
      digitalWrite(PHASE_B_PIN, HIGH);
      pinMode(PHASE_B_PIN, OUTPUT);
      digitalWrite(PHASE_C_PIN, LOW);
      pinMode(PHASE_C_PIN, OUTPUT);
      break;

    case 5: // A->C, B is floating
      pinMode(PHASE_B_PIN, INPUT);
      digitalWrite(PHASE_A_PIN, HIGH);
      pinMode(PHASE_A_PIN, OUTPUT);
      digitalWrite(PHASE_C_PIN, LOW);
      pinMode(PHASE_C_PIN, OUTPUT);
      break;
  }
}

// --- Timer1 Overflow Interrupt Service Routine ---
// This function is called automatically at the start of every PWM cycle.
ISR(TIMER1_OVF_vect) {
  // Check if we are still in the forced startup phase
  if (startup_step_counter < STARTUP_STEPS) {
    // During startup, we commutate at a fixed rate to get the motor spinning.
    unsigned long current_time = micros();
    if (current_time - last_commutation_time >= current_commutation_delay) {
      last_commutation_time = current_time;
      commutation_step = (commutation_step + 1) % 6;
      commutate();
      startup_step_counter++;

      // Gradually decrease the delay to speed up
      if (current_commutation_delay > 1000) {
        current_commutation_delay -= 20;
      }
    }
    return; // Exit ISR until startup is complete
  }


  // --- BEMF Detection Logic (after startup) ---
  int bemf_reading = 0;
  bool zero_crossed = false;

  // Read BEMF from the floating pin for the current step
  // We look for a rising or falling edge depending on the step
  switch (commutation_step) {
    case 0: // A->B, C floating (expect falling BEMF)
      bemf_reading = analogRead(PHASE_C_PIN);
      if (bemf_high && bemf_reading < BEMF_THRESHOLD_LOW) {
        zero_crossed = true;
        bemf_high = false;
      }
      break;
    case 1: // C->B, A floating (expect rising BEMF)
      bemf_reading = analogRead(PHASE_A_PIN);
      if (!bemf_high && bemf_reading > BEMF_THRESHOLD_HIGH) {
        zero_crossed = true;
        bemf_high = true;
      }
      break;
    case 2: // C->A, B floating (expect falling BEMF)
      bemf_reading = analogRead(PHASE_B_PIN);
      if (bemf_high && bemf_reading < BEMF_THRESHOLD_LOW) {
        zero_crossed = true;
        bemf_high = false;
      }
      break;
    case 3: // B->A, C floating (expect rising BEMF)
      bemf_reading = analogRead(PHASE_C_PIN);
      if (!bemf_high && bemf_reading > BEMF_THRESHOLD_HIGH) {
        zero_crossed = true;
        bemf_high = true;
      }
      break;
    case 4: // B->C, A floating (expect falling BEMF)
      bemf_reading = analogRead(PHASE_A_PIN);
      if (bemf_high && bemf_reading < BEMF_THRESHOLD_LOW) {
        zero_crossed = true;
        bemf_high = false;
      }
      break;
    case 5: // A->C, B floating (expect rising BEMF)
      bemf_reading = analogRead(PHASE_B_PIN);
      if (!bemf_high && bemf_reading > BEMF_THRESHOLD_HIGH) {
        zero_crossed = true;
        bemf_high = true;
      }
      break;
  }

  if (zero_crossed) {
    // Zero crossing detected! Now we wait for the 30-degree phase delay
    // before the next commutation.
    // The delay is half of the time since the last commutation.
    unsigned long time_since_last_comm = micros() - last_commutation_time;
    current_commutation_delay = time_since_last_comm / 2;
    
    // We can't use delay() in an ISR, so we set a flag to commutate in the main loop.
    // However, a better approach for timing critical code is to handle it all
    // within the interrupt context using timers. For simplicity, we will set a flag.
    // A more advanced method would use another timer to trigger the commutation.
    commutate_now = true; 
  }
}

void loop() {
  // The main loop's primary job is to trigger commutation when the ISR signals it.
  // This avoids doing complex logic inside the ISR.
  if (commutate_now) {
    // Wait for the calculated delay before commutating
    delayMicroseconds(current_commutation_delay);

    // Advance to the next step and commutate
    commutation_step = (commutation_step + 1) % 6;
    commutate();
    
    // Reset flag and record time
    commutate_now = false;
    last_commutation_time = micros();
  }

  // You can add code here to change the motor speed via Serial input
  if (Serial.available() > 0) {
    int new_speed = Serial.parseInt();
    if (new_speed > 50 && new_speed < PWM_TOP_VALUE) {
      motor_speed_pwm = new_speed;
      OCR1A = motor_speed_pwm; // Update PWM duty cycle
      Serial.print("Speed set to: ");
      Serial.println(new_speed);
    }
  }
}
