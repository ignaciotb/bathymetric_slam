# Bathymetric Graph SLAM

Baseline SLAM framework for underwater vehicles.
The algorithm gets a set of bathymetric submaps as input and corrects the global map constructed while refining the vehicle trajectory through a map-to-map registration followed by a pose graph optimization. 

## Dependencies (Ubuntu 16.04)
* AUVLIB [here](https://github.com/nilsbore/auvlib) 
* PCL  http://pointclouds.org/
* G2O https://github.com/RainerKuemmerle/g2o
* Ceres http://ceres-solver.org/


## Building

Clone this repository and create a `build` folder under the root, then execute
```
cd build
cmake ..
make -j4
```

### Available apps
There are two applications avalible under the `bin` folder.
If you compile the apps with the macro "INTERACTIVE" = 1 [here](https://github.com/ignaciotb/bathymetric_slam/tree/master/src/apps/src), you'll have to hit "q" for every step to be executed. This will allow you to visualize the before/after of the GICP registration per submap and the global graph optimization.

##### SLAM with simulated data
In order to test the framework with simulated data in the form of .pdc files, use the toy dataset `map_small` under `sim_data`. Specify the name of the output file with the graph in G2O format as a second argument. The current script optimizes the graph with Ceres, but with this file you can run your optimization in G2O if preferred. 
You can visualize both the ground truth map and vehicle trajectory in the visualizer. To start the optimization process, hit "q".
```
./test_slam_simulation --folder /path/to/folder/
```
The simulation outputs a measure of the error contained in the map, as well as the height maps and error plots as .png files.
To increase the complexity of the sim dataset, add more Gaussian noise to the vehicle's position estimate.
In order to adapt the performance of the algorithm to the dataset, tune the weights of the edges of the pose-graph, the GICP registration parameters and the Ceres solver.

##### SLAM with real data
To run the SLAM solution with real data from a bathymetric survey, currently the input is in the form of a cereal file containing all the necessary information from your data files.
You can find a real survey carried out with an ROV [here](https://strands.pdc.kth.se/public/IROS-2019-Bathymetry/). Download it, adjust the framework values, and test it.
```
./test_slam_real --file /path/to/datasets/your_data.cereal
```

### Generating your own data from the SMARC simulator
You can generate bigger and more complex bathymetric surveys from user-defined environments within the SMARC simulator.
To install it, follow the instructions [here](https://github.com/smarc-project/rosinstall).
To create and save the surveys with the default AUV and environment, run
```
roslaunch smarc_bringup auv_scenarios.launch
roslaunch smarc_bringup auv_system.launch
rosrun rviz rviz
```
Through the pop-up window, navigate the AUV over the area to map. You can press UP, DOWN,LEFT, RIGHT to steer, w to start the thruster, and s to stop it.
When you have finished with the survey, move the `.pdc` files with the submaps generated from `~/.ros/` to a new folder and run the framework as explained above.

### Generating your own cereal files from real surveys
Take a look at the [AUVLIB](https://github.com/nilsbore/auvlib) toolbox in order to parse MBES, SSS, navigation, etc data from the most common formats into .cereal files.
