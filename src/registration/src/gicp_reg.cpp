#include "registration/gicp_reg.hpp"

using namespace std;
using namespace Eigen;
using PointsT = std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd>>;



SubmapRegistration::SubmapRegistration(){

    ret_tf_ = Eigen::Matrix4f::Identity();
    benchmark_ = benchmark::track_error_benchmark();

    // Constrain GICP to x,y, yaw
    pcl::registration::WarpPointRigid3D<PointT, PointT>::Ptr warp_fcn
      (new pcl::registration::WarpPointRigid3D<PointT, PointT>);

    pcl::registration::TransformationEstimationLM<PointT, PointT>::Ptr te
            (new pcl::registration::TransformationEstimationLM<PointT, PointT>);
    te->setWarpFunction (warp_fcn);
    gicp_.setTransformationEstimation(te);

    // GICP parameters
    gicp_.setMaxCorrespondenceDistance(2);
    gicp_.setTransformationEpsilon (1e-3);
    gicp_.setMaximumIterations(100);
}

SubmapRegistration::~SubmapRegistration(){

}

SubmapObj SubmapRegistration::constructTrgSubmap(const SubmapsVec& submaps_set, std::vector<int>& overlaps){

    // Merge submaps in overlaps into submap_trg
    SubmapObj submap_trg;
    for(SubmapObj submap_j: submaps_set){
        if(std::find(overlaps.begin(), overlaps.end(), submap_j.submap_id_) != overlaps.end()){
            submap_trg.submap_pcl_ += submap_j.submap_pcl_;
        }
    }

    return submap_trg;
}

void SubmapRegistration::transformSubmap(SubmapObj& submap){

    // Apply tranform to submap frame
    Isometry3f rel_tf = Isometry3f (Isometry3f(Translation3f(ret_tf_.block<3,1>(0,3)))*
                                    Isometry3f(Quaternionf(ret_tf_.block<3,3>(0,0)).normalized()));

    submap.submap_tf_ = rel_tf * submap.submap_tf_;
}

double SubmapRegistration::consistencyErrorOverlap(const SubmapObj& trg_submap,
                                                   const SubmapObj& src_submap){

    // Compute consistency error in overlapped area
    PointsT submaps;
    submaps.push_back(trg_submap.submap_pcl_.getMatrixXfMap(3,4,0).transpose().cast<double>());
    submaps.push_back(src_submap.submap_pcl_.getMatrixXfMap(3,4,0).transpose().cast<double>());

    Eigen::MatrixXd error_vals;
    double consistency_rms_error;
    std::vector<std::vector<std::vector<MatrixXd>>> grid_maps = benchmark_.create_grids_from_matrices(submaps);
    tie(consistency_rms_error, error_vals) = benchmark_.compute_consistency_error(grid_maps);

    return (consistency_rms_error);

}


bool SubmapRegistration::gicpSubmapRegistration(SubmapObj& trg_submap, SubmapObj& src_submap){

    // Copy the originals to work over them
    PointCloudT::Ptr src_pcl_ptr (new PointCloudT(src_submap.submap_pcl_));
    PointCloudT::Ptr trg_pcl_ptr (new PointCloudT(trg_submap.submap_pcl_));

    // The Iterative Closest Point algorithm
    gicp_.setInputSource(src_pcl_ptr);
    gicp_.setInputTarget(trg_pcl_ptr);
    gicp_.align (src_submap.submap_pcl_);

    // Apply transform to submap
    ret_tf_ =  gicp_.getFinalTransformation();
    this->transformSubmap(src_submap);

    // TODO: compute this value properly
    bool convergence = (gicp_.hasConverged())? true: false;
    return convergence;
}

