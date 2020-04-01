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

#ifndef GRAPH_CONSTRUCTION_HPP
#define GRAPH_CONSTRUCTION_HPP

#include "g2o/types/slam3d/vertex_se3.h"
#include "g2o/types/slam3d/edge_se3.h"
#include "g2o/core/factory.h"

#include <Eigen/Core>
#include <Eigen/Dense>

#include "submaps_tools/submaps.hpp"

#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>

#include "graph_optimization/utils_g2o.hpp"

using namespace std;
using namespace Eigen;
using namespace g2o;

typedef std::vector<Eigen::Isometry3d, Eigen::aligned_allocator<Eigen::Isometry3d> > tf_vec;

class GraphConstructor{

private:

public:

    vector<VertexSE3*> vertices_;
    vector<EdgeSE3*> drEdges_;
    vector<EdgeSE3*> lcEdges_;
    std::vector<Eigen::Matrix2d, Eigen::aligned_allocator<Eigen::Matrix2d> > covs_lc_;
    tf_vec drMeas_;
    tf_vec lcMeas_;
    tf_vec drChain_;
    int edge_covs_type_;

    GraphConstructor(std::vector<Eigen::Matrix2d, Eigen::aligned_allocator<Eigen::Matrix2d> > covs_lc);

    ~GraphConstructor();

    void createNewVertex(SubmapObj& submap);

    void saveG2OFile(std::string outFilename);

    void findLoopClosures(SubmapObj &submap_i,
                          const SubmapsVec& submaps_set, double info_thres);

    void createLCEdge(const SubmapObj& submap_from, const SubmapObj& submap_to);

    void createInitialEstimate(SubmapsVec &submaps_set);

    void createDREdge(const SubmapObj& submap);

    void addNoiseToGraph(GaussianGen& transSampler, GaussianGen& rotSampler);

};

#endif // GRAPH_CONSTRUCTION_HPP
