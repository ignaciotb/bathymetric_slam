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

#ifndef GICP_REG_HPP
#define GICP_REG_HPP

#include <math.h>
#include <Eigen/Core>

#include <pcl/registration/gicp.h>
#include <pcl/registration/warp_point_rigid.h>
#include <pcl/registration/warp_point_rigid_3d.h>
#include <pcl/registration/transformation_estimation_lm.h>
#include <pcl/features/normal_3d.h>
#include <pcl/search/impl/search.hpp>
#include <pcl/Vertices.h>
#include <pcl/features/from_meshes.h>

#include "submaps_tools/submaps.hpp"

#include "data_tools/benchmark.h"

typedef std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d> > CovsVec;
typedef boost::shared_ptr <CovsVec > CovsVecPtr;

class SubmapRegistration {

private:
    pcl::GeneralizedIterativeClosestPoint<PointT, PointT> gicp_;
    Eigen::Matrix4f ret_tf_;
    benchmark::track_error_benchmark benchmark_;

public:

    SubmapRegistration();

    bool gicpSubmapRegistration(SubmapObj &trg_submap, SubmapObj &src_submap);

    void transformSubmap(SubmapObj& submap);

    SubmapObj constructTrgSubmap(const SubmapsVec& submaps_set, std::vector<int> &overlaps);

    double consistencyErrorOverlap(const SubmapObj& trg_submap, const SubmapObj& src_submap);

    ~SubmapRegistration();

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

};


#endif // GICP_REG_HPP
