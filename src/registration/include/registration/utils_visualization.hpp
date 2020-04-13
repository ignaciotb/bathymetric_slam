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

#ifndef UTILS_VISUALIZATION_HPP
#define UTILS_VISUALIZATION_HPP

#include <fstream>
#include <iostream>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/visualization/pcl_visualizer.h>

#include <Eigen/Core>

#include "submaps_tools/submaps.hpp"
#include "graph_optimization/graph_construction.hpp"
#include "graph_optimization/ceres_optimizer.hpp"


typedef pcl::PointCloud<pcl::PointXYZ> PointCloudT;
typedef pcl::PointCloud<pcl::PointXYZRGB> PointCloudRGB;
typedef pcl::PointXYZ PointT;
using pcl::visualization::PCLVisualizer;
typedef std::tuple<Eigen::Isometry3f, Eigen::Array3f> grid_point;
typedef std::vector<grid_point, Eigen::aligned_allocator<grid_point> > grid;

class SubmapsVisualizer{

private:

    int vp1_;
    int vp2_;
    PCLVisualizer viewer_;
    int sample_cnt;
    int traj_cnt_;

public:

    SubmapsVisualizer(PCLVisualizer& viewer);

    void setVisualizer(SubmapsVec &submaps_set, int num);

    void updateVisualizer(const SubmapsVec& submaps_set);

//    void plotPoseGraphG2O(const GraphConstructor& graph);

    void plotPoseGraphG2O(const GraphConstructor& graph, const SubmapsVec &submaps_set);

    void plotPoseGraphCeres(SubmapsVec &submaps_set);

    void addCoordinateSystem(const Eigen::Isometry3f& tf);

    void removeCoordinateSystem();

    void visualizeGrid(const grid& grid);

    void addSubmap(const SubmapObj& submap_i, int vp);

    void plotMBESPing(const SubmapObj& submap_i, float spam, float n_beams, int vp);
};

#endif // UTILS_VISUALIZATION_HPP
