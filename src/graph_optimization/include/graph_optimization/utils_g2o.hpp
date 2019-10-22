/* Copyright 2019 Ignacio Torroba (torroba@kth.se)
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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

//#include "graph_optimization/graph_construction.hpp"

using namespace std;
using namespace g2o;

typedef pcl::PointCloud<pcl::PointXYZ> PointCloudT;
typedef pcl::PointXYZ PointT;
typedef g2o::GaussianSampler<Eigen::Vector3d, Eigen::Matrix3d> GaussianGen;

Matrix<double, 6, 6> generateGaussianNoise(GaussianGen& transSampler,
                                           GaussianGen& rotSampler);

void addNoiseToSubmap(GaussianGen& transSampler,
                      GaussianGen& rotSampler,
                      SubmapObj& submap);

void addNoiseToMap(GaussianGen& transSampler,
                   GaussianGen& rotSampler,
                   SubmapsVec& submap_set);

#endif // UTILS_G2O_HPP
