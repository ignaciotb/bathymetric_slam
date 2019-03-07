#ifndef SUBMAPS_HPP
#define SUBMAPS_HPP

#include <fstream>
#include <iostream>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/common/transforms.h>
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

#include <pcl/common/projection_matrix.h>
#include <pcl/filters/uniform_sampling.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <eigen3/Eigen/Core>

#include "data_tools/std_data.h"

using namespace std;
using namespace Eigen;
typedef std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd>> PointsT;
typedef pcl::PointCloud<pcl::PointXYZ> PointCloudT;
typedef pcl::PointXYZ PointT;
typedef std::vector<Vector3d, aligned_allocator<Vector3d>> corners;

class SubmapObj{

private:

    std::pair<int, corners> getSubmapCorners(const SubmapObj& submap);

    bool checkSubmapsOverlap(const corners submap_i_corners, const corners submap_k_corners);

    bool pointToLine(const Vector3d seg_a, const Vector3d seg_b, const Vector3d point_c);

public:

    int submap_id_;
    unsigned int swath_id_;
    PointCloudT submap_pcl_;
    std::vector<int> overlaps_idx_;
    Eigen::Vector3d colors_;
    Eigen::Isometry3f submap_tf_;
    Eigen::Matrix<double,6,6> submap_info_;
    Eigen::MatrixXd auv_tracks_;

    SubmapObj();

    SubmapObj(const unsigned int& submap_id, const unsigned int& swath_id, PointCloudT& submap_pcl);

    void findOverlaps(std::vector<SubmapObj, Eigen::aligned_allocator<SubmapObj>> &submaps_set);

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

typedef std::vector<SubmapObj, Eigen::aligned_allocator<SubmapObj>> SubmapsVec;

void readSubmapFile(const string submap_str, PointCloudT::Ptr submap_pcl);

std::vector<std::string> checkFilesInDir(DIR *dir);

std::vector<SubmapObj, Eigen::aligned_allocator<SubmapObj> > readSubmapsInDir(const string& dir_path);

PointsT pclToMatrixSubmap(const SubmapsVec& submaps_set);

PointsT trackToMatrixSubmap(const SubmapsVec& submaps_set);

double computeInfoInSubmap(const SubmapObj& submap);

SubmapsVec parseSubmapsAUVlib(std_data::pt_submaps& ss);

PointsT trackofSubmap(const SubmapsVec& submaps_set);

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


#endif // SUBMAPS_HPP
