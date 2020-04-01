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

#ifndef CERES_OPTIMIZER_HPP
#define CERES_OPTIMIZER_HPP

#include <iostream>
#include <fstream>
#include <string>

#include "ceres/ceres.h"
#include "graph_optimization/read_g2o.h"
#include "graph_optimization/pose_graph_3d_error_term.h"
#include "graph_optimization/types.h"
#include "graph_optimization/graph_construction.hpp"

#include "submaps_tools/submaps.hpp"

#include "gflags/gflags.h"
#include "glog/logging.h"


namespace ceres{
namespace optimizer {

// Constructs the nonlinear least squares optimization problem from the pose
// graph constraints.
void BuildOptimizationProblem(const VectorOfConstraints& constraints,
                              MapOfPoses* poses, ceres::Problem* problem);

// Returns true if the solve was successful.
int SolveOptimizationProblem(ceres::Problem* problem);

// Output the poses to the file with format: id x y z q_x q_y q_z q_w.
bool OutputPoses(const std::string& filename, const MapOfPoses& poses);

MapOfPoses ceresSolver(const std::string& outFilename, const int drConstraints);

void updateSubmapsCeres(const MapOfPoses &poses, SubmapsVec& submaps_set);

void saveOriginalTrajectory(SubmapsVec& submaps_set);

}
}

#endif // CERES_OPTIMIZER_HPP
