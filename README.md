[//]: # (Image References)
[ekf]: ./ekf_end.png

![ekf]

# Extended Kalman Filter Project

This project implement an extended kalman filter to estimate the state of a moving object of interest with noisy lidar and radar measurements. 

# Current solution
The current solution have been developed using xcode 8.3.3 and currently delivers and RMSE about 0.94, 0.84, 0.34, 0.41.  
The code is fairely organized and it appears to run smoothly, basic optimization have been put in place.
For simplicity 2 simple switches useLaser, useRadar have been added to main.cpp to enable or disable the measurement coming from a specific sensor.
