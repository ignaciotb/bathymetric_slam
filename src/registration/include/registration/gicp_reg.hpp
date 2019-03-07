#ifndef GICP_REG_HPP
#define GICP_REG_HPP

#include <eigen3/Eigen/Core>
#include <pcl/registration/gicp.h>
#include <pcl/registration/warp_point_rigid.h>
#include <pcl/registration/warp_point_rigid_3d.h>
#include <pcl/registration/transformation_estimation_lm.h>

#include "submaps_tools/submaps.hpp"

#include "data_tools/benchmark.h"

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
