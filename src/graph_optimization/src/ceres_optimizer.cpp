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

#include "graph_optimization/ceres_optimizer.hpp"

namespace ceres{
namespace optimizer {

// Constructs the nonlinear least squares optimization problem from the pose
// graph constraints.
void BuildOptimizationProblem(const VectorOfConstraints& constraints,
                              MapOfPoses* poses, ceres::Problem* problem,
                              int drConstraints) {
    CHECK(poses != NULL);
    CHECK(problem != NULL);
    if (constraints.empty()) {
        LOG(INFO) << "No constraints, no problem to optimize.";
        return;
    }

//    LossFunction* loss_function = new ceres::HuberLoss(1);
    LossFunction* loss_function = nullptr;
    SubsetParameterization* z_local_param = new SubsetParameterization(3, std::vector<int>{2});
    SubsetParameterization* roll_pitch_local_param = new SubsetParameterization(3, std::vector<int>{0,1});

    // DR constraints
    std::cout << "Adding DR " << drConstraints <<  " constraints" << std::endl;
    for (VectorOfConstraints::const_iterator constraints_iter = constraints.begin();
         constraints_iter != constraints.begin() + drConstraints; ++constraints_iter){

        const Constraint3d& constraint = *constraints_iter;

        MapOfPoses::iterator pose_begin = poses->find(constraint.id_begin);
        CHECK(pose_begin != poses->end())
            << "Pose with ID: " << constraint.id_begin << " not found.";
        MapOfPoses::iterator pose_end = poses->find(constraint.id_end);
        CHECK(pose_end != poses->end())
            << "Pose with ID: " << constraint.id_end << " not found.";

        // Cost function
        const Eigen::Matrix<double, 6, 6> sqrt_information =
            constraint.information.llt().matrixL();
        ceres::CostFunction* cost_function =
            PoseGraph3dErrorTerm::Create(constraint.t_be, sqrt_information);

        // Residual block
        problem->AddResidualBlock(cost_function, loss_function,
                                  pose_begin->second.p.data(),
                                  pose_begin->second.q.data(),
                                  pose_end->second.p.data(),
                                  pose_end->second.q.data());

        // Constraints in roll and pitch
        problem->SetParameterization(pose_begin->second.q.data(), roll_pitch_local_param);
        problem->SetParameterization(pose_end->second.q.data(), roll_pitch_local_param);

        // Constraints in z
        problem->SetParameterization(pose_begin->second.p.data(), z_local_param);
        problem->SetParameterization(pose_end->second.p.data(), z_local_param);

        // Set boundaries for x, y and yaw
//        double upp_constraint = 10;
//        double low_constraint = 10;
//        double yaw_upp_const = M_PI/10.0;
//        double yaw_low_const = M_PI/10.0;
//        problem->SetParameterLowerBound(pose_begin->second.p.data(), 0, pose_begin->second.p[0] - low_constraint);
//        problem->SetParameterLowerBound(pose_end->second.p.data(), 0, pose_end->second.p[0] - low_constraint );
//        problem->SetParameterUpperBound(pose_begin->second.p.data(), 0, pose_begin->second.p[0] + upp_constraint);
//        problem->SetParameterUpperBound(pose_end->second.p.data(), 0, pose_end->second.p[0] + upp_constraint);

//        problem->SetParameterLowerBound(pose_begin->second.p.data(), 1, pose_begin->second.p[1] - low_constraint);
//        problem->SetParameterLowerBound(pose_end->second.p.data(), 1, pose_end->second.p[1] - low_constraint);
//        problem->SetParameterUpperBound(pose_begin->second.p.data(), 1, pose_begin->second.p[1] + upp_constraint);
//        problem->SetParameterUpperBound(pose_end->second.p.data(), 1, pose_end->second.p[1] + upp_constraint);

//        problem->SetParameterLowerBound(pose_begin->second.q.data(), 2, pose_begin->second.q[2] - yaw_low_const);
//        problem->SetParameterLowerBound(pose_end->second.q.data(), 2, pose_end->second.q[2] - yaw_low_const);
//        problem->SetParameterUpperBound(pose_begin->second.q.data(), 2, pose_begin->second.q[2] + yaw_upp_const);
//        problem->SetParameterUpperBound(pose_end->second.q.data(), 2, pose_end->second.q[2] + yaw_upp_const);
    }

    // LC constraints
    if(constraints.size() - drConstraints != 0){
        std::cout << "Adding LC " << constraints.size() - drConstraints << " constraints" << std::endl;
        for (VectorOfConstraints::const_iterator constraints_iter = constraints.begin() + drConstraints+1;
             constraints_iter != constraints.end(); ++constraints_iter){

            const Constraint3d& constraint = *constraints_iter;

            MapOfPoses::iterator pose_begin = poses->find(constraint.id_begin);
            CHECK(pose_begin != poses->end())
                << "Pose with ID: " << constraint.id_begin << " not found.";
            MapOfPoses::iterator pose_end = poses->find(constraint.id_end);
            CHECK(pose_end != poses->end())
                << "Pose with ID: " << constraint.id_end << " not found.";

            const Eigen::Matrix<double, 6, 6> sqrt_information =
                constraint.information.llt().matrixL();

            // Cost function
            ceres::CostFunction* cost_function =
                PoseGraph3dErrorTerm::Create(constraint.t_be, sqrt_information);

            // Residual block
            problem->AddResidualBlock(cost_function, loss_function,
                                      pose_begin->second.p.data(),
                                      pose_begin->second.q.data(),
                                      pose_end->second.p.data(),
                                      pose_end->second.q.data());

            // Constraints in roll and pitch
            problem->SetParameterization(pose_begin->second.q.data(), roll_pitch_local_param);
            problem->SetParameterization(pose_end->second.q.data(), roll_pitch_local_param);

            // Constraints in z
            problem->SetParameterization(pose_begin->second.p.data(), z_local_param);
            problem->SetParameterization(pose_end->second.p.data(), z_local_param);

            // Set boundaries for x, y and yaw
//            double upp_constraint = 10;
//            double low_constraint = 10;
//            double yaw_upp_const = M_PI/100.0;
//            double yaw_low_const = M_PI/100.0;
    //        problem->SetParameterLowerBound(pose_begin->second.p.data(), 0, pose_begin->second.p[0] - low_constraint);
    //        problem->SetParameterLowerBound(pose_end->second.p.data(), 0, pose_end->second.p[0] - low_constraint );
    //        problem->SetParameterUpperBound(pose_begin->second.p.data(), 0, pose_begin->second.p[0] + upp_constraint);
    //        problem->SetParameterUpperBound(pose_end->second.p.data(), 0, pose_end->second.p[0] + upp_constraint);

    //        problem->SetParameterLowerBound(pose_begin->second.p.data(), 1, pose_begin->second.p[1] - low_constraint);
    //        problem->SetParameterLowerBound(pose_end->second.p.data(), 1, pose_end->second.p[1] - low_constraint);
    //        problem->SetParameterUpperBound(pose_begin->second.p.data(), 1, pose_begin->second.p[1] + upp_constraint);
    //        problem->SetParameterUpperBound(pose_end->second.p.data(), 1, pose_end->second.p[1] + upp_constraint);

//            problem->SetParameterLowerBound(pose_begin->second.q.data(), 2, pose_begin->second.q[2] - yaw_low_const);
//            problem->SetParameterLowerBound(pose_end->second.q.data(), 2, pose_end->second.q[2] - yaw_low_const);
//            problem->SetParameterUpperBound(pose_begin->second.q.data(), 2, pose_begin->second.q[2] + yaw_upp_const);
//            problem->SetParameterUpperBound(pose_end->second.q.data(), 2, pose_end->second.q[2] + yaw_upp_const);
        }
    }

    // Constrain the gauge freedom: set first AUV pose as anchor constant
    MapOfPoses::iterator pose_start_iter = poses->begin();
    CHECK(pose_start_iter != poses->end()) << "There are no poses.";
    problem->SetParameterBlockConstant(pose_start_iter->second.p.data());
    problem->SetParameterBlockConstant(pose_start_iter->second.q.data());
}

// Returns true if the solve was successful.
int SolveOptimizationProblem(ceres::Problem* problem) {
    CHECK(problem != NULL);

    ceres::Solver::Options options;
    options.max_num_iterations = 300;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
//    options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
//    options.linear_solver_type = ceres::SPARSE_SCHUR;
//    options.linear_solver_type = ceres::SPARSE_QR;
    //options.linear_solver_type = ceres::DENSE_QR;

    //  options.minimizer_type = ceres::LINE_SEARCH;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.use_inner_iterations = false;
    options.minimizer_progress_to_stdout = true;
    options.num_threads =4;
    options.eta = 1e-9;

    ceres::Solver::Summary summary;
    ceres::Solve(options, problem, &summary);
//    std::cout << summary.FullReport() << '\n';


    return summary.iterations.size();
}

// Output the poses to the file with format: id x y z q_x q_y q_z q_w.
bool OutputPoses(const std::string& filename, const MapOfPoses& poses) {
    std::fstream outfile;
    outfile.open(filename.c_str(), std::istream::out);
    if (!outfile) {
        LOG(ERROR) << "Error opening the file: " << filename;
        return false;
    }
    for (std::map<int, Pose3d, std::less<int>,
                Eigen::aligned_allocator<std::pair<const int, Pose3d> > >::
           const_iterator poses_iter = poses.begin();
       poses_iter != poses.end(); ++poses_iter) {
    const std::map<int, Pose3d, std::less<int>,
                   Eigen::aligned_allocator<std::pair<const int, Pose3d> > >::
        value_type& pair = *poses_iter;

    Eigen::Quaterniond q = eulerToQuat<double>(pair.second.q);
    outfile << pair.first << " " << pair.second.p.transpose() << " "
            << q.x() << " " << q.y() << " "
            << q.z() << " " << q.w() << '\n';
    }
    return true;
}

MapOfPoses ceresSolver(const std::string& outFilename, const int drConstraints){
    // Ceres solver
    ceres::optimizer::MapOfPoses poses;
    ceres::optimizer::VectorOfConstraints constraints;

    CHECK(ceres::optimizer::ReadG2oFile(outFilename, &poses, &constraints))
        << "Error reading the file: " << outFilename;

    CHECK(ceres::optimizer::OutputPoses("poses_corrupted.txt", poses))
        << "Error outputting to poses_corrupted.txt";

    std::cout << "Original poses output" << std::endl;

    ceres::Problem problem;
    ceres::optimizer::BuildOptimizationProblem(constraints, &poses, &problem, drConstraints);

    std::cout << "Ceres problem built" << std::endl;

    int iterations = ceres::optimizer::SolveOptimizationProblem(&problem);
//        << "The solve was not successful, exiting.";

    CHECK(ceres::optimizer::OutputPoses("poses_optimized.txt", poses))
        << "Error outputting to poses_optimized.txt";

    return poses;
}

void updateSubmapsCeres(const ceres::optimizer::MapOfPoses& poses, SubmapsVec& submaps_set){

    // Update pointclouds
    unsigned int i = 0;
    for(SubmapObj& submap: submaps_set){
        // Final pose of submap_i
        ceres::optimizer::Pose3d pose_i = poses.at(i);
        Eigen::Quaterniond q = Eigen::AngleAxisd(pose_i.q.x(), Eigen::Vector3d::UnitX())
                               * Eigen::AngleAxisd(pose_i.q.y(), Eigen::Vector3d::UnitY())
                               * Eigen::AngleAxisd(pose_i.q.z(), Eigen::Vector3d::UnitZ());

        Isometry3f final_tf = (Isometry3f) q.cast<float>();
        final_tf.translation() = pose_i.p.cast<float>();

        // Transform submap_i pcl and tf
        pcl::transformPointCloud(submap.submap_pcl_, submap.submap_pcl_,
                                 (final_tf * submap.submap_tf_.inverse()).matrix());
        submap.submap_tf_ = final_tf;
        i++;
    }
}


void saveOriginalTrajectory(SubmapsVec& submaps_set){

    covs covs_lc;
    GraphConstructor* graph = new GraphConstructor(covs_lc);

    for(SubmapObj& submap_i: submaps_set){
        graph->createNewVertex(submap_i);

        // Create DR edge i and store (skip submap 0)
        if(submap_i.submap_id_ != 0 ){
            graph->createDREdge(submap_i);
        }
    }

    string outFilename = "graph_original.g2o";
    graph->saveG2OFile(outFilename);
    ceres::optimizer::MapOfPoses poses;
    ceres::optimizer::VectorOfConstraints constraints;
    CHECK(ceres::optimizer::ReadG2oFile(outFilename, &poses, &constraints))
        << "Error reading the file: " << outFilename;

    CHECK(ceres::optimizer::OutputPoses("poses_original.txt", poses))
        << "Error outputting to poses_original.txt";
}


}  // namespace optimizer
}  // namespace ceres
