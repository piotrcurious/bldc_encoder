 // Enhanced Three-Phase Resolver/Encoder with Advanced Signal Processing
// Optimized for Arduino Uno/Nano with improved performance and reliability

#include <math.h>

// === HARDWARE CONFIGURATION ===
const int PHASE_A_PIN = A0;
const int PHASE_B_PIN = A1;
const int PHASE_C_PIN = A2;
const int DITHER_PIN = 9;

// === SYSTEM CONSTANTS ===
const float VREF = 5.0f;
const int ADC_MAX = 1023;
const float OFFSET_VOLTAGE = VREF * 0.5f;
const float DEG_PER_RAD = 57.2957795f; // More precise than 180/PI
const float TWO_PI = 6.28318531f;
const float PI_FLOAT = 3.14159265f;

// === TIMING CONFIGURATION ===
const float DT = 0.01f; // 100Hz update rate
const unsigned long DT_MICROS = 10000UL; // 10ms in microseconds

// === SIGNAL PROCESSING PARAMETERS ===
const float LAMBDA = 0.15f; // Momentum smoothing factor
const float MAX_DITHER_AMP = 12.0f; // Max PWM dither amplitude
const float DITHER_FREQ = 100.0f; // Dither frequency in Hz
const int MAX_OVERSAMPLE = 8;

// === ADAPTIVE THRESHOLDS ===
const float NOISE_THRESHOLD = 0.01f; // Minimum signal level
const float HIGH_SIGNAL_THRESHOLD = 0.5f;
const float MED_SIGNAL_THRESHOLD = 0.2f;
const float LOW_SIGNAL_THRESHOLD = 0.05f;

// === COMPLEX NUMBER CLASS ===
struct Complex {
  float re, im;

  // Constructor
  Complex(float r = 0.0f, float i = 0.0f) : re(r), im(i) {}

  // Operators
  Complex operator+(const Complex& other) const {
    return Complex(re + other.re, im + other.im);
  }

  Complex operator*(float scalar) const {
    return Complex(re * scalar, im * scalar);
  }

  // Fast magnitude calculation using approximation for low-precision needs
  float magnitude() const {
    return sqrtf(re * re + im * im);
  }

  float magnitudeSquared() const {
    return re * re + im * im;
  }

  // Angle calculation with proper quadrant handling
  float angle() const {
    return atan2f(im, re);
  }

  // Normalization
  Complex normalize() const {
    float mag = magnitude();
    if (mag > NOISE_THRESHOLD) {
      return Complex(re / mag, im / mag);
    }
    return Complex(0, 0);
  }
};

// === GLOBAL VARIABLES ===
Complex z_current(0, 0);
Complex momentum(0, 0);
float previous_theta = 0.0f;
float angular_velocity = 0.0f;
unsigned long last_update_time = 0;
unsigned long program_start_time = 0;

// Performance monitoring
unsigned long loop_count = 0;
float max_loop_time = 0.0f;
float avg_loop_time = 0.0f;

// === SETUP FUNCTION ===
void setup() {
  Serial.begin(115200);
  while (!Serial) { ; } // Wait for serial port to connect
  
  // Initialize dither pin
  pinMode(DITHER_PIN, OUTPUT);
  analogWrite(DITHER_PIN, 127); // Set to midpoint
  
  // Initialize timing
  program_start_time = micros();
  last_update_time = program_start_time;
  
  // ADC optimization for faster conversions
  // Set ADC prescaler to 16 (default is 128)
  // This increases ADC speed at slight cost of precision
  ADCSRA &= ~PS_128;
  ADCSRA |= PS_16;
  
  Serial.println(F("Enhanced Three-Phase Resolver Initialized"));
  Serial.println(F("Theta(deg)\tOmega(rad/s)\tMagnitude\tCorrelation\tDither\tOS\tLoop(us)"));
}

// === FAST ADC READ FUNCTION ===
inline int fastAnalogRead(int pin) {
  // Direct register manipulation for faster ADC reads
  ADMUX = (1 << REFS0) | (pin & 0x07);
  ADCSRA |= (1 << ADSC);
  while (ADCSRA & (1 << ADSC));
  return ADC;
}

// === VOLTAGE CONVERSION FUNCTION ===
inline float adcToVoltage(int adc_value) {
  return (adc_value * VREF / ADC_MAX) - OFFSET_VOLTAGE;
}

// === ADAPTIVE OVERSAMPLING FUNCTION ===
int calculateOversampleCount(float signal_magnitude) {
  if (signal_magnitude < LOW_SIGNAL_THRESHOLD) return 8;
  if (signal_magnitude < MED_SIGNAL_THRESHOLD) return 4;
  if (signal_magnitude < HIGH_SIGNAL_THRESHOLD) return 2;
  return 1;
}

// === CLARKE TRANSFORM FUNCTION ===
Complex clarkeTransform(float Va, float Vb, float Vc) {
  // Optimized Clarke transform coefficients
  const float SQRT3_2 = 0.866025404f; // sqrt(3)/2
  const float TWO_THIRDS = 0.666666667f;
  
  float Valpha = TWO_THIRDS * (Va - 0.5f * Vb - 0.5f * Vc);
  float Vbeta = TWO_THIRDS * SQRT3_2 * (Vb - Vc);
  
  return Complex(Valpha, Vbeta);
}

// === DITHER GENERATION FUNCTION ===
void generateDither(float time_seconds, float signal_magnitude) {
  // Adaptive dither amplitude based on signal strength
  float dither_amplitude = MAX_DITHER_AMP / (1.0f + 20.0f * signal_magnitude);
  
  // Generate triangle wave using optimized calculation
  float phase = time_seconds * DITHER_FREQ;
  float triangle = 2.0f * fabsf(2.0f * (phase - floorf(phase + 0.5f))) - 1.0f;
  
  // Apply dither to PWM
  int pwm_value = 127 + (int)(triangle * dither_amplitude);
  pwm_value = constrain(pwm_value, 0, 255);
  analogWrite(DITHER_PIN, pwm_value);
}

// === ANGLE WRAPPING FUNCTION ===
float wrapAngle(float angle) {
  while (angle > PI_FLOAT) angle -= TWO_PI;
  while (angle < -PI_FLOAT) angle += TWO_PI;
  return angle;
}

// === MAIN LOOP ===
void loop() {
  unsigned long loop_start = micros();
  unsigned long current_time = micros();
  
  // Check if it's time for the next update
  if (current_time - last_update_time < DT_MICROS) {
    return; // Skip this iteration to maintain consistent timing
  }
  
  float elapsed_time = (current_time - program_start_time) * 1e-6f; // Convert to seconds
  
  // === PHASE 1: QUICK SIGNAL ASSESSMENT ===
  float Va_est = adcToVoltage(fastAnalogRead(PHASE_A_PIN));
  float Vb_est = adcToVoltage(fastAnalogRead(PHASE_B_PIN));
  float Vc_est = adcToVoltage(fastAnalogRead(PHASE_C_PIN));
  
  Complex z_estimate = clarkeTransform(Va_est, Vb_est, Vc_est);
  float estimated_magnitude = z_estimate.magnitude();
  
  // === PHASE 2: ADAPTIVE OVERSAMPLING ===
  int oversample_count = calculateOversampleCount(estimated_magnitude);
  
  float Va_sum = 0, Vb_sum = 0, Vc_sum = 0;
  
  for (int i = 0; i < oversample_count; ++i) {
    Va_sum += adcToVoltage(fastAnalogRead(PHASE_A_PIN));
    Vb_sum += adcToVoltage(fastAnalogRead(PHASE_B_PIN));
    Vc_sum += adcToVoltage(fastAnalogRead(PHASE_C_PIN));
  }
  
  // Average the samples
  float Va_avg = Va_sum / oversample_count;
  float Vb_avg = Vb_sum / oversample_count;
  float Vc_avg = Vc_sum / oversample_count;
  
  // === PHASE 3: CLARKE TRANSFORM ===
  z_current = clarkeTransform(Va_avg, Vb_avg, Vc_avg);
  float current_magnitude = z_current.magnitude();
  
  // === PHASE 4: MOMENTUM FILTERING ===
  momentum = momentum * (1.0f - LAMBDA) + z_current * LAMBDA;
  
  // === PHASE 5: ANGLE CALCULATION ===
  float theta = momentum.angle();
  if (theta < 0) theta += TWO_PI;
  
  // === PHASE 6: VELOCITY ESTIMATION ===
  float delta_theta = wrapAngle(theta - previous_theta);
  float actual_dt = (current_time - last_update_time) * 1e-6f;
  angular_velocity = delta_theta / actual_dt;
  previous_theta = theta;
  
  // === PHASE 7: CORRELATION CALCULATION ===
  float momentum_magnitude = momentum.magnitude();
  float correlation = 0.0f;
  if (current_magnitude > NOISE_THRESHOLD && momentum_magnitude > NOISE_THRESHOLD) {
    correlation = (z_current.re * momentum.re + z_current.im * momentum.im) / 
                  (current_magnitude * momentum_magnitude);
  }
  
  // === PHASE 8: DITHER GENERATION ===
  generateDither(elapsed_time, current_magnitude);
  
  // === PHASE 9: PERFORMANCE MONITORING ===
  unsigned long loop_time = micros() - loop_start;
  if (loop_time > max_loop_time) max_loop_time = loop_time;
  avg_loop_time = (avg_loop_time * loop_count + loop_time) / (loop_count + 1);
  loop_count++;
  
  // === PHASE 10: OUTPUT ===
  Serial.print(theta * DEG_PER_RAD, 2);
  Serial.print(F("\t\t"));
  Serial.print(angular_velocity, 2);
  Serial.print(F("\t\t"));
  Serial.print(current_magnitude, 3);
  Serial.print(F("\t\t"));
  Serial.print(correlation, 3);
  Serial.print(F("\t\t"));
  Serial.print(MAX_DITHER_AMP / (1.0f + 20.0f * current_magnitude), 2);
  Serial.print(F("\t"));
  Serial.print(oversample_count);
  Serial.print(F("\t"));
  Serial.println(loop_time);
  
  // Performance report every 1000 loops
  if (loop_count % 1000 == 0) {
    Serial.print(F("Performance - Max: "));
    Serial.print(max_loop_time);
    Serial.print(F("us, Avg: "));
    Serial.print(avg_loop_time, 1);
    Serial.println(F("us"));
  }
  
  last_update_time = current_time;
}

// === ADC PRESCALER DEFINITIONS ===
#define PS_16  (1 << ADPS2)
#define PS_128 (1 << ADPS2) | (1 << ADPS1) | (1 << ADPS0)
