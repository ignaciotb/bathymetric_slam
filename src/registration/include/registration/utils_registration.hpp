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

#ifndef UTILS_HPP
#define UTILS_HPP

#include <fstream>
#include <iostream>

#include <pcl/common/transforms.h>
#include <pcl/point_types.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/registration/icp.h>

#include <Eigen/Geometry>
#include <Eigen/Core>



using namespace std;
typedef pcl::PointCloud<pcl::PointXYZ> PointCloudT;
typedef pcl::PointXYZ PointT;


struct submap {
public:   
    int submap_id;
    PointCloudT submap_pcl;
    Eigen::Matrix4f tf_submap_map;
    std::vector<Eigen::VectorXf, Eigen::aligned_allocator<Eigen::VectorXf>> auv_pose_;   // AUV pose for each ping
    std::vector<Eigen::Matrix3f, Eigen::aligned_allocator<Eigen::Matrix3f>> pcl_covariances;
    Eigen::MatrixXf cov_submap_frame;   // 6x6 cov of submap frame
    unsigned int pings_num_;
    unsigned int beams_per_ping_;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


struct icp_match{
public:

    icp_match(PointT trg_point, PointT src_point, Eigen::Vector3f normal, Eigen::Vector3f error_mean, Eigen::Matrix3f error_sigma){
        trg_point_ = trg_point;
        src_point_ = src_point;
        normal_ = normal;
        error_mean_ = error_mean;
        error_sigma_ = error_sigma;
    }

    PointT trg_point_;
    PointT src_point_;
    Eigen::Vector3f normal_;
    // Components of error pdf
    Eigen::Vector3f error_mean_;
    Eigen::Matrix3f error_sigma_;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

void subsampleMap(submap& submap_input);

void outlierFilter(submap& submap_input);

Eigen::Vector3f computePCAPcl(PointCloudT& set_Ai);

Eigen::Affine3d create_rotation_matrix(double ax, double ay, double az);

Eigen::Matrix4f inverseTfMatrix(Eigen::Matrix4f tf_mat);

void plotSubmapsSet(const std::vector<submap, Eigen::aligned_allocator<submap> > &submaps_set);

void print4x4Matrix (const Eigen::Matrix4f & matrix);

#endif // UTILS_HPP
