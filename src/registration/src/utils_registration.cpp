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

#include "registration/utils_registration.hpp"

Eigen::Affine3d create_rotation_matrix(double ax, double ay, double az) {
  Eigen::Affine3d rx = Eigen::Affine3d(Eigen::AngleAxisd(ax, Eigen::Vector3d(1, 0, 0)));
  Eigen::Affine3d ry = Eigen::Affine3d(Eigen::AngleAxisd(ay, Eigen::Vector3d(0, 1, 0)));
  Eigen::Affine3d rz = Eigen::Affine3d(Eigen::AngleAxisd(az, Eigen::Vector3d(0, 0, 1)));
  return rz * ry * rx;
}


Eigen::Matrix4f inverseTfMatrix(Eigen::Matrix4f tf_mat){

    Eigen::Matrix3f R_inv = tf_mat.topLeftCorner(3,3).transpose();

    Eigen::Matrix4f tf_mat_inv = Eigen::Matrix4f::Identity();
    tf_mat_inv.topLeftCorner(3,3) = R_inv;
    tf_mat_inv.topRightCorner(3,1) = R_inv * (-tf_mat.topRightCorner(3,1));

    return tf_mat_inv;

}


void plotSubmapsSet(const std::vector<submap, Eigen::aligned_allocator<submap>>& submaps_set){
    // Plot submaps
    pcl::visualization::PCLVisualizer viewer ("ICP demo");
    viewer.setSize (1920, 1080);  // Visualiser window size

    PointCloudT::Ptr submap_i_ptr;
    int cnt = 0;
    for(const submap& submap_i: submaps_set){
//        if(submap_i.submap_id == 0 || submap_i.submap_id == 24){
            submap_i_ptr.reset(new PointCloudT(submap_i.submap_pcl));
            pcl::transformPointCloud(*submap_i_ptr, *submap_i_ptr, inverseTfMatrix(submap_i.tf_submap_map));
            viewer.addPointCloud (submap_i_ptr, "cloud_" + std::to_string(cnt), 0);
//        }
        ++cnt;
    }

    while (!viewer.wasStopped()) {
        viewer.spinOnce ();
    }
}


void print4x4Matrix (const Eigen::Matrix4f & matrix) {
    printf ("Rotation matrix :\n");
    printf ("    | %6.3f %6.3f %6.3f | \n", matrix (0, 0), matrix (0, 1), matrix (0, 2));
    printf ("R = | %6.3f %6.3f %6.3f | \n", matrix (1, 0), matrix (1, 1), matrix (1, 2));
    printf ("    | %6.3f %6.3f %6.3f | \n", matrix (2, 0), matrix (2, 1), matrix (2, 2));
    printf ("Translation vector :\n");
    printf ("t = < %6.3f, %6.3f, %6.3f >\n\n", matrix (0, 3), matrix (1, 3), matrix (2, 3));
}

Eigen::Vector3f computePCAPcl(PointCloudT& set_Ai){

    // Compute the mean of the PCL
    Eigen::Vector4f com_Ai;
    pcl::compute3DCentroid(set_Ai, com_Ai);

    // Compute covariance matrix
    Eigen::Matrix3f cov_mat;
    pcl::computeCovarianceMatrixNormalized (set_Ai, com_Ai, cov_mat);
//    pcl::computeCovarianceMatrix(set_Ai_demean, cov_mat);

    // Extract eigenvalues and eigenvector from cov matrix
    Eigen::EigenSolver<Eigen::Matrix3f> eigenSolver;
    eigenSolver.compute(cov_mat, true);

    Eigen::MatrixXf eigenVectors = eigenSolver.eigenvectors().real();
    Eigen::VectorXf eigenvalues = eigenSolver.eigenvalues().real();

    // Return eigenvector with smallest eigenvalue
    std::vector<std::pair<double, int>> pidx;
    for (unsigned int i = 0 ; i<cov_mat.cols(); i++){
        pidx.push_back(std::make_pair(eigenvalues(i), i));
    }
    sort(pidx.begin(), pidx.end());

    return eigenVectors.col(pidx[0].second).transpose();
}


void subsampleMap(submap& submap_input){

    // Compute PCA of input pointcloud
    Eigen::Vector3f plane_normal = computePCAPcl(submap_input.submap_pcl);

    // Create best fitting plane
    PointT sum_pointsj(0,0,0);
    double sum_traces = 0;
    double prj_trace;
    for(unsigned int i=0; i<submap_input.submap_pcl.points.size(); i++){
        prj_trace = submap_input.pcl_covariances.at(i).trace();
        sum_traces += std::pow(prj_trace, -2);
        sum_pointsj.getArray3fMap() += submap_input.submap_pcl.points.at(i).getArray3fMap() * (float)std::pow(prj_trace, -2);
    }

    Eigen::Vector3f p_nu = Eigen::Vector3f(sum_pointsj.getArray3fMap() * 1/sum_traces);
    double d_p = plane_normal.dot(p_nu);

    // Distance of every point to their projection on the best fitting plane
    Eigen::Vector3f point_prj, distance;
    double average_dist = 0;
    for(const PointT& pointj: submap_input.submap_pcl.points){
         point_prj =  Eigen::Vector3f(pointj.getArray3fMap()) - (Eigen::Vector3f(pointj.getArray3fMap()).dot(plane_normal) - d_p) * plane_normal;
         distance = Eigen::Vector3f(point_prj(0) - pointj.x,
                                    point_prj(1) - pointj.y,
                                    point_prj(2) - pointj.z);
         average_dist += distance.norm();
    }

    // Average of the distances
    average_dist = average_dist / submap_input.submap_pcl.size();

    // Filter out points closer than average
    int cnt = 0;
    submap submap_aux;
    for(const PointT& pointj: submap_input.submap_pcl.points){
        point_prj =  Eigen::Vector3f(pointj.getArray3fMap()) - (Eigen::Vector3f(pointj.getArray3fMap()).dot(plane_normal) - d_p) * plane_normal;
        distance = Eigen::Vector3f(point_prj(0) - pointj.x,
                                   point_prj(1) - pointj.y,
                                   point_prj(2) - pointj.z);
        if(distance.norm() >= average_dist*1.5){
            submap_aux.submap_pcl.push_back(pointj);
            submap_aux.pcl_covariances.push_back(submap_input.pcl_covariances.at(cnt));
        }
        ++cnt;
    }

    // Clear and store final output values
    submap_input.submap_pcl.clear();
    submap_input.pcl_covariances.clear();

    submap_input.submap_pcl = submap_aux.submap_pcl;
    submap_input.pcl_covariances = submap_aux.pcl_covariances;
}


void outlierFilter(submap &submap_input){

    // Build Kdtree
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    PointCloudT::Ptr pcl_ptr;
    pcl_ptr.reset(new PointCloudT(submap_input.submap_pcl));
    kdtree.setInputCloud(pcl_ptr);

    // Average distance between NN points in pcl
    int average_nn = 0;
    int K = 4;
    std::vector<int> pointIdxRadiusSearch(K);
    std::vector<float> pointRadiusSquaredDistance(K);

    double radius = 1;
    for(PointT point_i: submap_input.submap_pcl.points){
        average_nn += kdtree.radiusSearch (point_i, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance);
    }
    average_nn = average_nn / submap_input.submap_pcl.points.size();
    std::cout << "Average number of nn: " << average_nn << std::endl;

    // Filter out points farther from any other than average dist
    submap submap_aux;
    for(PointT point_i: submap_input.submap_pcl.points){
        if(kdtree.radiusSearch (point_i, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) >= average_nn*0.2){
            submap_aux.submap_pcl.push_back(point_i);
        }
    }

    submap_input.submap_pcl.clear();
    submap_input.submap_pcl = submap_aux.submap_pcl;
}

