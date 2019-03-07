#include "registration/utils_visualization.hpp"
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */


using namespace std;
using pcl::visualization::PointCloudColorHandlerCustom;

SubmapsVisualizer::SubmapsVisualizer(pcl::visualization::PCLVisualizer &viewer):
                                    viewer_(viewer){

}


void SubmapsVisualizer::setVisualizer(SubmapsVec& submaps_set, int num){

    // The color we will be using
    float black = 0.0;  // Black
    float white = 1.0 - black;

    // Create one/two vertically separated viewports
    if (num == 2){
        // Second viewport
        viewer_.createViewPort (0.0, 0.0, 0.5, 1.0, vp1_);
        viewer_.createViewPort (0.5, 0.0, 1.0, 1.0, vp2_);
        viewer_.addText ("Results", 10, 15, 16, white, white, white, "info_2", vp2_);
        viewer_.setBackgroundColor (white, white, white, vp2_);
        viewer_.addCoordinateSystem(5.0, "reference_frame", vp2_);
    }
    else{
        viewer_.createViewPort (0.0, 0.0, 1.0, 1.0, vp1_);
    }


    // Adding text descriptions in each viewport
    viewer_.addText ("Ground truth", 10, 15, 16, white, white, white, "info_1", vp1_);

    // Set GT pointclouds
    unsigned int i = 0;
    PointCloudT::Ptr submap_ptr (new PointCloudT);
    for(SubmapObj& submap: submaps_set){
        submap_ptr.reset(new PointCloudT(submap.submap_pcl_));
        submap.colors_ = Vector3d(rand() % 256, rand() % 256, rand() % 256);
        PointCloudColorHandlerCustom<PointT> cloud_color(submap_ptr, submap.colors_[0], submap.colors_[1], submap.colors_[2]);
        viewer_.addPointCloud(submap_ptr, cloud_color, "gt_cloud_" + std::to_string(i), vp1_);
        viewer_.addCoordinateSystem(3.0, submap.submap_tf_, "gt_cloud_" + std::to_string(i), vp1_);
        i++;
    }

    // Set GT trajectory (linearized)
    Vector3i dr_color = Vector3i(rand() % 256, rand() % 256, rand() % 256);
    for(unsigned int j=0; j<submaps_set.size()-1; j++){
        SubmapObj submap_fr = submaps_set.at(j+1);
        Eigen::Vector3d from_ps = submap_fr.submap_tf_.translation().cast<double>();
        SubmapObj submap_to = submaps_set.at(j);
        Eigen::Vector3d to_ps = submap_to.submap_tf_.translation().cast<double>();
        viewer_.addArrow(PointT(from_ps[0],from_ps[1],from_ps[2]), PointT(to_ps[0],to_ps[1],to_ps[2]),
                dr_color[0], dr_color[1], dr_color[2], false, "gt_dr_edge_" + std::to_string(j), vp1_);
    }

    // Set background color
    viewer_.setBackgroundColor (white, white, white, vp1_);

    // Set camera position and orientation
    viewer_.setSize (1920/2, 1080/2);
    viewer_.addCoordinateSystem(5.0, "gt_reference_frame", vp1_);
}


void SubmapsVisualizer::updateVisualizer(const SubmapsVec& submaps_set){

    viewer_.removeAllPointClouds(vp2_);
    viewer_.removeAllCoordinateSystems(vp2_);
    viewer_.addCoordinateSystem(5.0, "reference_frame", vp2_);

    // Update pointclouds
    unsigned int i = 0;
    PointCloudT::Ptr submap_ptr (new PointCloudT);
    for(const SubmapObj& submap: submaps_set){
        submap_ptr.reset(new PointCloudT(submap.submap_pcl_));
        PointCloudColorHandlerCustom<PointT> cloud_color(submap_ptr, submap.colors_[0], submap.colors_[1], submap.colors_[2]);
        viewer_.addPointCloud(submap_ptr, cloud_color, "cloud_" + std::to_string(i), vp2_);
        viewer_.addCoordinateSystem(3.0, submap.submap_tf_, "cloud_" + std::to_string(i), vp2_);
        i++;
    }
    viewer_.spinOnce();
}


void SubmapsVisualizer::plotPoseGraphG2O(const GraphConstructor& graph){

    // TODO: check to and from vertices
    Vector3i dr_color = Vector3i(rand() % 256, rand() % 256, rand() % 256);
    unsigned int i = 0;
    for(EdgeSE3* edgeDR: graph.drEdges_){
        VertexSE3* from = static_cast<VertexSE3*>(edgeDR->vertex(1));
        VertexSE3* to = static_cast<VertexSE3*>(edgeDR->vertex(0));
        Eigen::Vector3d from_ps = from->estimate().translation();
        Eigen::Vector3d to_ps = to->estimate().translation();
        viewer_.addArrow(PointT(from_ps[0],from_ps[1],from_ps[2]), PointT(to_ps[0],to_ps[1],to_ps[2]),
                dr_color[0], dr_color[1], dr_color[2], false, "dr_edge_" + std::to_string(i), vp2_);
        i++;
    }

    i = 0;
    Vector3i lc_color = Vector3i(rand() % 256, rand() % 256, rand() % 256);
    for(EdgeSE3* edgeLC: graph.edges_){
        VertexSE3* from = static_cast<VertexSE3*>(edgeLC->vertex(0));
        VertexSE3* to = static_cast<VertexSE3*>(edgeLC->vertex(1));
        Eigen::Vector3d from_ps = from->estimate().translation();
        Eigen::Vector3d to_ps = to->estimate().translation();
        viewer_.addArrow(PointT(from_ps[0],from_ps[1],from_ps[2]), PointT(to_ps[0],to_ps[1],to_ps[2]),
                lc_color[0], lc_color[1], lc_color[2], false, "lc_edge_" + std::to_string(i), vp2_);
        i++;
    }
    viewer_.spinOnce();
}

void SubmapsVisualizer::plotPoseGraphCeres(const ceres::optimizer::MapOfPoses& poses, SubmapsVec& submaps_set){

    // Clean initial graph
    viewer_.removeAllPointClouds(vp2_);
    viewer_.removeAllCoordinateSystems(vp2_);
    viewer_.removeAllShapes(vp2_);
    viewer_.addCoordinateSystem(5.0, "reference_frame", vp2_);

    // Update pointclouds
    unsigned int i = 0;
    PointCloudT::Ptr submap_ptr (new PointCloudT);
    for(SubmapObj& submap: submaps_set){
        // Final pose of submap_i
        ceres::optimizer::Pose3d pose_i = poses.at(i);
        Eigen::Quaterniond q = Eigen::AngleAxisd(pose_i.q.x(), Eigen::Vector3d::UnitX())
                               * Eigen::AngleAxisd(pose_i.q.y(), Eigen::Vector3d::UnitY())
                               * Eigen::AngleAxisd(pose_i.q.z(), Eigen::Vector3d::UnitZ());

        Isometry3f final_tf = (Isometry3f) q.cast<float>();
        final_tf.translation() = pose_i.p.cast<float>();

        // Transform submap_i pcl and tf
        pcl::transformPointCloud(submap.submap_pcl_, submap.submap_pcl_, final_tf.matrix() * submap.submap_tf_.inverse().matrix());
        submap.submap_tf_ = final_tf;

        // Plot pcl and tf
        submap_ptr.reset(new PointCloudT(submap.submap_pcl_));
        PointCloudColorHandlerCustom<PointT> cloud_color(submap_ptr, submap.colors_[0], submap.colors_[1], submap.colors_[2]);
        viewer_.addPointCloud(submap_ptr, cloud_color, "cloud_" + std::to_string(i), vp2_);
        viewer_.addCoordinateSystem(3.0, submap.submap_tf_, "cloud_" + std::to_string(i), vp2_);
        i++;
    }

    // Plot final trajectory estimate
    Vector3i dr_color = Vector3i(rand() % 256, rand() % 256, rand() % 256);
    for(unsigned int j=0; j<submaps_set.size()-1; j++){
        SubmapObj submap_fr = submaps_set.at(j+1);
        Eigen::Vector3d from_ps = submap_fr.submap_tf_.translation().cast<double>();
        SubmapObj submap_to = submaps_set.at(j);
        Eigen::Vector3d to_ps = submap_to.submap_tf_.translation().cast<double>();
        viewer_.addArrow(PointT(from_ps[0],from_ps[1],from_ps[2]), PointT(to_ps[0],to_ps[1],to_ps[2]),
                dr_color[0], dr_color[1], dr_color[2], false, "final_dr_edge_" + std::to_string(j), vp2_);
    }

    viewer_.spinOnce();
}

