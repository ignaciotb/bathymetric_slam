#ifndef UTILS_VISUALIZATION_HPP
#define UTILS_VISUALIZATION_HPP

#include <fstream>
#include <iostream>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/visualization/pcl_visualizer.h>

#include <eigen3/Eigen/Core>

#include "submaps_tools/submaps.hpp"
#include "graph_optimization/graph_construction.hpp"
#include "graph_optimization/ceres_optimizer.hpp"

typedef pcl::PointCloud<pcl::PointXYZ> PointCloudT;
typedef pcl::PointXYZ PointT;
using pcl::visualization::PCLVisualizer;

class SubmapsVisualizer{

private:

    int vp1_;
    int vp2_;
    PCLVisualizer viewer_;

public:

    SubmapsVisualizer(PCLVisualizer& viewer);

    void setVisualizer(SubmapsVec &submaps_set, int num);

    void updateVisualizer(const SubmapsVec& submaps_set);

    void plotPoseGraphG2O(const GraphConstructor& graph);

    void plotPoseGraphCeres(const ceres::optimizer::MapOfPoses& poses, SubmapsVec &submaps_set);

};

#endif // UTILS_VISUALIZATION_HPP
