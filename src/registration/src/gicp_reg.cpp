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
    gicp_.setMaxCorrespondenceDistance(10);
    gicp_.setTransformationEpsilon (1e-4);
//    gicp_.setEuclideanFitnessEpsilon(1e-5);
    gicp_.setMaximumIterations(200);
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

    // Compute GICP covs externally
    pcl::NormalEstimation<PointT, pcl::Normal> ne;
    pcl::search::KdTree<PointT>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ> ());
    ne.setSearchMethod(tree);
    ne.setRadiusSearch(20);

    pcl::PointCloud<pcl::Normal>::Ptr normals_src(new pcl::PointCloud<pcl::Normal>);
    ne.setInputCloud(src_pcl_ptr);
    ne.compute(*normals_src);

    pcl::PointCloud<pcl::Normal>::Ptr normals_trg(new pcl::PointCloud<pcl::Normal>);
    ne.setInputCloud(trg_pcl_ptr);
    ne.compute(*normals_trg);

    CovsVec covs_src;
    CovsVecPtr covs_src_ptr;
    pcl::features::computeApproximateCovariances(*src_pcl_ptr, *normals_src, covs_src);
    covs_src_ptr.reset(new CovsVec(covs_src));

    CovsVec covs_trg;
    CovsVecPtr covs_trg_ptr;
    pcl::features::computeApproximateCovariances(*trg_pcl_ptr, *normals_trg, covs_trg);
    covs_trg_ptr.reset(new CovsVec(covs_trg));

    // Compute GICP vanilla information matrix
    Eigen::Matrix3d avg_cov_src;
    avg_cov_src.setZero();
    for(unsigned int n=0; n<covs_src_ptr->size(); n++){
        avg_cov_src += covs_src_ptr->at(n);
    }
    avg_cov_src /= (covs_src_ptr->size());

    Eigen::Matrix3d avg_cov_trg;
    avg_cov_trg.setZero();
    for(unsigned int n=0; n<covs_trg_ptr->size(); n++){
        avg_cov_trg += covs_trg_ptr->at(n);
    }
    avg_cov_trg /= (covs_trg_ptr->size());

    // The Iterative Closest Point algorithm
    gicp_.setInputSource(src_pcl_ptr);
    gicp_.setInputTarget(trg_pcl_ptr);
    gicp_.setSourceCovariances(covs_src_ptr);
    gicp_.setTargetCovariances(covs_trg_ptr);
    gicp_.align (src_submap.submap_pcl_);

    // Apply transform to submap
    ret_tf_ =  gicp_.getFinalTransformation();
    this->transformSubmap(src_submap);

    // Set GICP information matrix
    src_submap.submap_lc_info_.setZero();
    Eigen::VectorXd info_diag(4);
    info_diag << 10000.0, 10000.0, 10000.0, 1000.0;
    src_submap.submap_lc_info_.bottomRightCorner(4,4) = info_diag.asDiagonal();

    Eigen::Matrix3d rot = ret_tf_.topLeftCorner(3,3).cast<double>();
    Eigen::Matrix3d gicp_cov = avg_cov_trg + rot*avg_cov_src*rot.transpose();
    src_submap.submap_lc_info_.topLeftCorner(2,2) = gicp_cov.topLeftCorner(2,2).inverse();
//    std::cout << gicp_cov.topLeftCorner(3,3) << std::endl;

    for(int x=0; x<src_submap.submap_lc_info_.array().size(); x++){
        if(isnan(src_submap.submap_lc_info_.array()(x))){
            throw std::runtime_error("Nan components in the matrix");
            std::exit(0);
        }
    }

    // TODO: compute this value properly
    bool convergence = (gicp_.hasConverged())? true: false;
    return convergence;
}

