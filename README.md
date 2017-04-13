#Extended Kalman Filter Project

##Purpose

This project implements an Extended Kalman Filter for tracking a pedestrian or other target with both Lidar and Radar.  This probabilistic algorithm first predicts position and process covariance, then updates these values using a sensor measurement.  The Kalman Filter state tracks position and velocity in two dimensions (four total dimensions).

##Handling Lidar _and_ Radar

For the purposes of this project, prediction is a linear process for both sensors.  The measurement step is linear for Lidar, but nonlinear for Radar.  Both the standard measurement step and the Extended Kalman Filter measurement step are used, depending on the sensor supplying each measurement, and a single 4-dimensional state vector is maintained.

##Screenshot from Udacity Kalman Tracker
![tracker] (images/Screenshot_from_2017-04-11.png "Udacity Tracking Screenshot")
