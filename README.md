# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program
   
### Simulator.
You can download the Simulator which contains the Path Planning Project from [here](https://github.com/udacity/self-driving-car-sim/releases).

### Goals
In this project our goal is to safely navigate around a virtual highway with other traffic that is driving +-10 MPH of the 50 MPH speed limit. We are provided with the car's localization and sensor fusion data, there is also a sparse map list of waypoints around the highway. The car should try to go as close as possible to the 50 MPH speed limit, which means passing slower traffic when possible, note that other cars will try to change lanes too. The car should avoid hitting other cars at all cost as well as driving inside of the marked road lanes at all times, unless going from one lane to another. The car should be able to make one complete loop around the 6946m highway. Since the car is trying to go 50 MPH, it should take a little over 5 minutes to complete 1 loop. Also the car should not experience total acceleration over 10 m/s^2 and jerk that is greater than 50 m/s^3.

#### The map of the highway is in data/highway_map.txt
Each waypoint in the list contains  [x,y,s,dx,dy] values. x and y are the waypoint's map coordinate position, the s value is the distance along the road to get to that waypoint in meters, the dx and dy values define the unit normal vector pointing outward of the highway loop.

The highway's waypoints loop around so the frenet s value, distance along the road, goes from 0 to 6945.554.

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./path_planning`.

## Trajectory Generation
To be able to change lanes without exceeding the max jerk, I use the spline interpolation library to generate a smooth set of points between given start and intermediate points. For each iteration, two start points are obtained from car's previous coordinates. I convert the start points to frenet coordinates to obtain s and d values. This helps me derive a set of three end points of the trajectory (more on this later). I convert all the five points of the trajectory to car's coordinates before generating spline. This is done because the library has a limitation that x should be monotonically increasing. Once we have the spline, we can generate a set of 50 coordinates to pass it to the simulator.



