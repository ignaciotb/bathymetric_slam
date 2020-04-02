# Bathymetric Graph SLAM

Baseline SLAM framework for underwater vehicles.
The algorithm gets a set of bathymetric submaps as input and corrects the global map constructed while refining the vehicle trajectory through a map-to-map registration followed by a pose graph optimization. 


![real_data_example](https://github.com/ignaciotb/bathymetric_slam/blob/master/img/graph_borno.png)


## Paper introducing and applying the method
The method implemented is described [in this paper](https://ieeexplore.ieee.org/abstract/document/8968241) and used [in this one](https://arxiv.org/abs/2003.10931)
```
@inproceedings{torroba2019towards,
  title={Towards Autonomous Industrial-Scale Bathymetric Surveying},
  author={Torroba, Ignacio and Bore, Nils and Folkesson, John},
  booktitle={2019 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)},
  pages={6377--6382},
  year={2019},
  organization={IEEE}
}

@article{torroba2020pointnetkl,
  title={PointNetKL: Deep Inference for GICP Covariance Estimation in Bathymetric SLAM},
  author={Torroba, Ignacio and Sprague, Christopher Iliffe and Bore, Nils and Folkesson, John},
  journal={arXiv preprint arXiv:2003.10931},
  year={2020}
}
```

## Dependencies (tested on Ubuntu 16.04 and 18.04)
* AUVLIB [here](https://github.com/ignaciotb/auvlib) 
* PCL  http://pointclouds.org/
* G2O https://github.com/RainerKuemmerle/g2o
* Ceres http://ceres-solver.org/

Note that for G2O to be used by this repo you need to install it at a system level.
From the G2O build folder, run  
```
sudo make install
```

## Building

Clone this repository and create a `build` folder under the root, then execute
```
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ..
make -j4
make install
```

Finally, add the following line to your ~/.bashrc file adapted to your own installation
```
export PATH=$PATH:/path/to/folder/bathymetric_slam/install/share
```
### Available apps
Under `bin` folder.
The process outputs .png images with the maps of bathymetry and consistency error.
The current script optimizes the graph with Ceres, but the app outputs a "graph.g2o" file which you can solve with G2O if preferred. 

##### SLAM with simulated data
In order to test the framework with data from the [SMARC simulator](https://github.com/smarc-project), use the toy dataset `map_small` under `sim_data`. 
You can visualize both the ground truth map and vehicle trajectory in the visualizer. To start the optimization process, hit "q".
```
./bathy_slam_real --simulation yes --bathy_survey ../sim_data/map_small/
```
The simulation outputs a measure of the error contained in the map, as well as the height maps and error plots as .png files.
To increase the complexity of the sim dataset, increase the Gaussian noise to the vehicle's position estimate.
In order to adapt the performance of the algorithm to the dataset, adjust the weights of the edges of the pose-graph accordingly and tune the GICP and the Ceres solver parameters.
The algorithm is by default tuned for the toy example `map_small`.
You can find a bigger and more challenging dataset from the simulator [here](https://strands.pdc.kth.se/public/IROS-2019-Bathymetry/).

##### SLAM with real data
To run the SLAM solution with real data from a bathymetric survey, currently the input is in the form of a cereal file containing all the necessary information from your data files.
You can find a real survey carried out with an ROV [here](https://strands.pdc.kth.se/public/IROS-2019-Bathymetry/). Download it, adjust the framework values, and test it.
```
./bathy_slam_real --simulation no --bathy_survey /path/to/datasets/mbes_pings.cereal 
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
Take a look at the [AUVLIB](https://github.com/nilsbore/auvlib) toolbox in order to parse real MBES, SSS, navigation, etc data from the most common formats into .cereal files.
