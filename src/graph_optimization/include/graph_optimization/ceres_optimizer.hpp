#ifndef CERES_OPTIMIZER_HPP
#define CERES_OPTIMIZER_HPP

#include <iostream>
#include <fstream>
#include <string>

#include "ceres/ceres.h"
#include "graph_optimization/read_g2o.h"
#include "graph_optimization/pose_graph_3d_error_term.h"
#include "graph_optimization/types.h"

#include "gflags/gflags.h"
#include "glog/logging.h"


namespace ceres{
namespace optimizer {

// Constructs the nonlinear least squares optimization problem from the pose
// graph constraints.
void BuildOptimizationProblem(const VectorOfConstraints& constraints,
                              MapOfPoses* poses, ceres::Problem* problem);

// Returns true if the solve was successful.
bool SolveOptimizationProblem(ceres::Problem* problem);

// Output the poses to the file with format: id x y z q_x q_y q_z q_w.
bool OutputPoses(const std::string& filename, const MapOfPoses& poses);

MapOfPoses ceresSolver(const std::string& outFilename);


}
}

#endif // CERES_OPTIMIZER_HPP
