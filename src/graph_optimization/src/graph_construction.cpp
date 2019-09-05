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

GraphConstructor::GraphConstructor(std::vector<Eigen::Matrix2d, Eigen::aligned_allocator<Eigen::Matrix2d> > covs_lc):
    covs_lc_(covs_lc){
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

    drEdges_.push_back(e);
    drMeas_.push_back(t);
}


void GraphConstructor::createLCEdge(const SubmapObj& submap_from, const SubmapObj& submap_to){

    std::cout << "LC edge from " << submap_from.submap_id_ << " to " << submap_to.submap_id_ << std::endl;

    // Generate loop closure edges
    VertexSE3* from = vertices_[submap_from.submap_id_];
    VertexSE3* to   = vertices_[submap_to.submap_id_];
    Eigen::Isometry3d t = from->estimate().inverse() * submap_to.submap_tf_.cast<double>();
    EdgeSE3* e = new EdgeSE3;
    e->setVertex(0, from);
    e->setVertex(1, to);
    e->setMeasurement(t);

    // Information matrix for LC edges
    Eigen::Array3f info = computeInfoInSubmap(submap_from);

    // Submap covariance matrix
    Eigen::Matrix3f cov_matrix;
    Eigen::Vector4f xyz_centroid;
    pcl::compute3DCentroid(submap_from.submap_pcl_, xyz_centroid);
    pcl::computeCovarianceMatrix(submap_from.submap_pcl_, xyz_centroid, cov_matrix);

    Eigen::Matrix<double, 6, 6> information = Eigen::Matrix<double, 6, 6>::Zero();

    // Info matrix proportional to variance in Z in the pointcloud
    if(covs_lc_.empty()){
        Eigen::VectorXd info_diag(3), info_diag_trans(3);
        double z = cov_matrix.normalized().inverse().cast<double>().row(2)(2);
        info_diag << 10000.0, 10000.0, 1000.0;
//        info_diag_trans << z, z, 10000.0;
        info_diag_trans << 1, 1, 10000.0;
        information.block<3,3>(0,0) = info_diag_trans.asDiagonal();
        information.block<3,3>(3,3) = info_diag.asDiagonal();
    }
    else{
        // Info matrix from NN training
        Eigen::Matrix2d cov_reg = this->covs_lc_.at(submap_from.submap_id_);
        std::cout << cov_reg << std::endl;
        Eigen::VectorXd info_diag(3), info_diag_trans(3);
        info_diag << 10000.0, 10000.0, 1000.0;
        information.topLeftCorner(2,2) = cov_reg.inverse();
        information(2,2) = 10000;
        information.block<3,3>(3,3) = info_diag.asDiagonal();
    }

    std::cout << information << std::endl;
    e->setInformation(information);

    // Check resulting COV is positive semi-definite
    Eigen::LLT<Eigen::MatrixXd> lltOfA(information);
    if(lltOfA.info() == Eigen::NumericalIssue){
        throw std::runtime_error("Non positive semi-definite CoV");
        std::exit(0);
    }

    lcEdges_.push_back(e);
    lcMeas_.push_back(t);
}


void GraphConstructor::findLoopClosures(SubmapObj& submap_i, const SubmapsVec& submaps_set,
                                        double info_thres){

    // Check all submaps in overlaps vector
    for(SubmapObj submap_j: submaps_set){
        if(find(submap_i.overlaps_idx_.begin(), submap_i.overlaps_idx_.end(), submap_j.submap_id_)
                != submap_i.overlaps_idx_.end()){
            this->createLCEdge(submap_i, submap_j);
        }
    }
}


void GraphConstructor::createInitialEstimate(SubmapsVec& submaps_set){

    // Concatenate all the odometry constraints to compute the initial kinematic chain
    for (size_t i =0; i < drEdges_.size(); i++) {
        Eigen::Isometry3d meas_i = drMeas_.at(i);
        EdgeSE3* e = drEdges_[i];
        VertexSE3* from = static_cast<VertexSE3*>(e->vertex(0));
        VertexSE3* to = static_cast<VertexSE3*>(e->vertex(1));

        Eigen::Isometry3d estimate_i = from->estimate() * meas_i;

//        std::cout << "From estimate " << from->estimate().translation().transpose() << std::endl;
//        std::cout << "Measurement " << meas_i.translation().transpose() << std::endl;
        Eigen::Matrix4f tf = (estimate_i.cast<float>() *
                              submaps_set.at(i+1).submap_tf_.cast<float>().inverse()).matrix();

//        std::cout << "Estimate " << i << " " << std::endl;
//        std::cout << tf << std::endl;
//        std::cout << "------" << std::endl;

        // Transform submap_i pcl and tf
        pcl::transformPointCloud(submaps_set.at(i+1).submap_pcl_, submaps_set.at(i+1).submap_pcl_,
                                 (estimate_i.cast<float>() * submaps_set.at(i+1).submap_tf_.cast<float>().inverse()).matrix());

        submaps_set.at(i+1).submap_tf_ = estimate_i.cast<float>();
        to->setEstimate(estimate_i);

        drChain_.push_back(estimate_i);
    }
}


void GraphConstructor::saveG2OFile(std::string outFilename){

    // Save graph to output file
    ofstream fileOutputStream;
    if (outFilename != "-") {
      cerr << "Writing into " << outFilename << endl;
      fileOutputStream.open(outFilename.c_str());
    } else {
      cerr << "writing to stdout" << endl;
    }
    ostream& fout = outFilename != "-" ? fileOutputStream : std::cout;

    // Concatenate DR and LC edges (DR go first, according to g2o convention)
    vector<EdgeSE3*> edges;
    edges.insert(edges.begin(), lcEdges_.begin(), lcEdges_.end());
    edges.insert(edges.begin(), drEdges_.begin(), drEdges_.end());
    tf_vec meas;
    meas.insert(meas.begin(), lcMeas_.begin(), lcMeas_.end());
    meas.insert(meas.begin(), drMeas_.begin(), drMeas_.end());

    string vertexTag = Factory::instance()->tag(vertices_[0]);
    string edgeTag = Factory::instance()->tag(edges[0]);

    for (size_t i = 0; i < vertices_.size(); ++i) {
      VertexSE3* v = vertices_[i];
      fout << vertexTag << " " << v->id() << " ";
      v->write(fout);
      fout << endl;
    }

    // Output
    for (size_t i = 0; i < edges.size(); ++i) {
      EdgeSE3* e = edges[i];
      VertexSE3* from = static_cast<VertexSE3*>(e->vertex(0));
      VertexSE3* to = static_cast<VertexSE3*>(e->vertex(1));
      fout << edgeTag << " " << from->id() << " " << to->id() << " ";
//      Vector7 meas=g2o::internal::toVectorQT(e->measurement());
      Vector7 meas_i=g2o::internal::toVectorQT(meas.at(i));
      for (int i=0; i<7; i++) fout  << meas_i[i] << " ";
      for (int i=0; i<e->information().rows(); i++){
        for (int j=i; j<e->information().cols(); j++) {
          fout <<  e->information()(i,j) << " ";
        }
      }
      fout << endl;
    }
}

