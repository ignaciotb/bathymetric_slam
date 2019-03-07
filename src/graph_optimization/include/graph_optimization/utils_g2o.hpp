#ifndef UTILS_G2O_HPP
#define UTILS_G2O_HPP

#include <fstream>
#include <iostream>

#include "g2o/stuff/sampler.h"

#include <eigen3/Eigen/Core>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/transforms.h>

#include "submaps_tools/submaps.hpp"

using namespace std;
using namespace g2o;

typedef pcl::PointCloud<pcl::PointXYZ> PointCloudT;
typedef pcl::PointXYZ PointT;
typedef g2o::GaussianSampler<Eigen::Vector3d, Eigen::Matrix3d> GaussianGen;

Matrix<double, 6, 6> generateGaussianNoise(GaussianGen& transSampler,
                                           GaussianGen& rotSampler);

void addNoiseToSubmap(GaussianGen& transSampler,
                      GaussianGen& rotSampler,
                      const Matrix<double, 6, 6> &information,
                      SubmapObj& submap, Isometry3f &poseDRt);


#endif // UTILS_G2O_HPP
