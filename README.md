# bldc_encoder
junkbox allowing use of 3 phase stepper motor from hdd as encoder 
Tune EKF Parameters (Q_ANGLE_EKF, Q_VEL_EKF, R_ANGLE_EKF): These are critical for the performance of the second layer.
Q_VEL_EKF (Process Noise for Velocity): This is very important. If your motor changes speed frequently or is subject to external disturbances, increase this. If it runs at a fairly constant speed, keep it smaller. A larger Q_VEL_EKF allows the estimated velocity to change more quickly, making the system more responsive but potentially more susceptible to noise.
R_ANGLE_EKF (Measurement Noise for Inferred Angle): This reflects how much you trust the atan2 derived angle. If the g_valA/B/C values are very noisy even after the first KF, or if the atan2 model isn't perfect for your motor, a higher R_ANGLE_EKF makes the EKF rely more on its internal velocity prediction and less on the instantaneous angle measurement.
