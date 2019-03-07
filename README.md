# Bathymetric Graph SLAM

SLAM framework for underwater vehicles.
The algorithm gets a set of bathymetric submaps as input and corrects the global map constructed while refining the vehicle trajectory through a map-to-map registration followed by a pose-graph optimization. 

## Dependencies (Ubuntu 16.04)
* AUVLIB [here](https://github.com/nilsbore/auvlib) 
* PCL  http://pointclouds.org/
* G2O https://github.com/RainerKuemmerle/g2o


## Building

Clone this repository and create a `build` folder under the root, then execute
```
cd build
cmake ..
make -j4
```

### Available apps
There are two applications avalible under the `bin` folder.
In order to test the framework with simulated data in the form of .pdc files, use the toy dataset `map_small` under `sim_data`.
You can visualize both the ground truth map and vehicle trajectory in the visualizer. To start the optimization process, hit "q".
```
./test_slam_simulation /path/to/repository/sim_data/map_small/
```

To run the SLAM solution with real data from a bathymetric survey, currently the input is in the form of a cereal file containing all the necessary information from your bathymetric files.
Check the AUVLIB [data project](https://github.com/nilsbore/auvlib/tree/master/example_projects/data_project) for tools to parse your files. 
```
./test_slam_real --folder /path/to/datasets/"your_data".cereal
```

### Generating your own data from the SMARC simulator
You can generate bigger and more complex bathymetric surveys from user-defined environments within the SMARC simulator.
To install it, follow the instructions [here](https://github.com/smarc-project/rosinstall).
To create and save the surveys with the default AUV and environment, run
```
roslaunch smarc_bringup auv_scenarios.launch
roslaunch smarc_bringup auv_system.launch
```
Through the pop-up window, navigate the AUV over the area to map. You can press UP, DOWN,LEFT, RIGHT to steer, w to start the thruster, and s to stop it.
When you have finished with the survey, move the `.pdc` files with the submaps generated from `~/.ros/` to `sim_data` and run the framework as explained above.