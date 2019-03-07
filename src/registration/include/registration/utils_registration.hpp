#ifndef UTILS_HPP
#define UTILS_HPP

#include <fstream>
#include <iostream>

#include <pcl/common/transforms.h>
#include <pcl/point_types.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/registration/icp.h>

#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Core>



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
