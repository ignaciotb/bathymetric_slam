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

#include "graph_optimization/graph_construction.hpp"

using namespace Eigen;
using namespace std;
using namespace g2o;

GraphConstructor::GraphConstructor(){

}

GraphConstructor::~GraphConstructor(){

}

void GraphConstructor::createNewVertex(SubmapObj& submap){

    VertexSE3* v = new VertexSE3;
    v->setId(submap.submap_id_);
    Eigen::Isometry3d t = submap.submap_tf_.cast<double>();
    v->setEstimate(t);
    vertices_.push_back(v);
}


void GraphConstructor::createDREdge(const SubmapObj& submap){

    // Generate dead reckoning edges
    VertexSE3* prev = vertices_[submap.submap_id_-1];
    VertexSE3* cur  = vertices_[submap.submap_id_];
    Eigen::Isometry3d t = prev->estimate().inverse() * cur->estimate();
    EdgeSE3* e = new EdgeSE3;
    e->setVertex(0, prev);
    e->setVertex(1, cur);
    e->setMeasurement(t);
    e->setInformation(submap.submap_info_);

    std::cout << "DR edge from " << submap.submap_id_ -1 << " to " << submap.submap_id_<< std::endl;
    drEdges_.push_back(e);
}


void GraphConstructor::createLCEdge(const SubmapObj& submap_from, const SubmapObj& submap_to){

    // Generate loop closure edges
    VertexSE3* from = vertices_[submap_from.submap_id_];
    VertexSE3* to   = vertices_[submap_to.submap_id_];
    Eigen::Isometry3d t = from->estimate().inverse() * submap_to.submap_tf_.cast<double>();
    EdgeSE3* e = new EdgeSE3;
    e->setVertex(0, from);
    e->setVertex(1, to);
    e->setMeasurement(t);

    // Information matrix for LC edges
    double info = computeInfoInSubmap(submap_from);

    // Submap covariance matrix
    Eigen::Matrix3f cov_matrix;
    Eigen::Vector4f xyz_centroid;
    pcl::compute3DCentroid(submap_from.submap_pcl_, xyz_centroid);
    pcl::computeCovarianceMatrix(submap_from.submap_pcl_, xyz_centroid, cov_matrix);

    // For testing
    Eigen::VectorXd info_diag(3);
    info_diag << 10000.0, 10000.0, 10000.0;
    Eigen::Matrix<double, 6, 6> information = Eigen::Matrix<double, 6, 6>::Zero();
    information.block<3,3>(0,0) = info_diag.asDiagonal()*0.1;
    information.block<3,3>(3,3) = info_diag.asDiagonal();
    e->setInformation(information);

    std::cout << "LC edge from " << submap_from.submap_id_ << " to " << submap_to.submap_id_ << std::endl;
    edges_.push_back(e);
}


void GraphConstructor::findLoopClosures(SubmapObj& submap_i, const SubmapsVec& submaps_set,
                                        double info_thres){

    // Check all submaps in overlaps vector
    for(SubmapObj submap_j: submaps_set){
        if(find(submap_i.overlaps_idx_.begin(), submap_i.overlaps_idx_.end(), submap_j.submap_id_)
                != submap_i.overlaps_idx_.end()){
            this->createLCEdge(submap_j, submap_i);
        }
    }
}


void GraphConstructor::createInitialEstimate(){

    // Concatenate all the odometry constraints to compute the initial state
    std::cout << "Number of DR edges " << drEdges_.size() << std::endl;
    for (size_t i =0; i < drEdges_.size(); ++i) {
        std::cout << "DR edge " << i << std::endl;
        EdgeSE3* e = drEdges_[i];
        VertexSE3* from = static_cast<VertexSE3*>(e->vertex(1));
        VertexSE3* to = static_cast<VertexSE3*>(e->vertex(0));
        HyperGraph::VertexSet aux;
        aux.insert(from);
        e->initialEstimate(aux, to);
    }
}


void GraphConstructor::saveG2OFile(std::string outFilename){

    // Save graph to output file
    ofstream fileOutputStream;
    if (outFilename == "-") {
        outFilename = "default_graph.g2o";
    }
    std::cout << "Saving graph to file " << outFilename << std::endl;

    fileOutputStream.open(outFilename.c_str());
    ostream& fout = outFilename != "-" ? fileOutputStream : std::cout;

    // Concatenate DR and LC edges (DR go first, according to g2o convention)
    edges_.insert(edges_.begin(), drEdges_.begin(), drEdges_.end());

    string vertexTag = Factory::instance()->tag(vertices_[0]);
    string edgeTag = Factory::instance()->tag(edges_[0]);

    for (size_t i = 0; i < vertices_.size(); ++i) {
      VertexSE3* v = vertices_[i];
      fout << vertexTag << " " << v->id() << " ";
      v->write(fout);
      fout << endl;
    }

    // Output
    for (size_t i = 0; i < edges_.size(); ++i) {
      EdgeSE3* e = edges_[i];
      VertexSE3* from = static_cast<VertexSE3*>(e->vertex(0));
      VertexSE3* to = static_cast<VertexSE3*>(e->vertex(1));
      fout << edgeTag << " " << from->id() << " " << to->id() << " ";
      e->write(fout);
      fout << endl;
    }
}

