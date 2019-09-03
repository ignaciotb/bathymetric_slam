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

#include "graph_optimization/utils_g2o.hpp"

using namespace std;
using namespace g2o;
using namespace Eigen;

Matrix<double, 6,6> generateGaussianNoise(GaussianGen& transSampler,
                                          GaussianGen& rotSampler){

    bool randomSeed = false;
    std::vector<double> noiseTranslation;
    std::vector<double> noiseRotation;
    noiseTranslation.push_back(2);
    noiseTranslation.push_back(2);
    noiseTranslation.push_back(0.0001);
    noiseRotation.push_back(0.0001);
    noiseRotation.push_back(0.0001);
    noiseRotation.push_back(0.05);

    Eigen::Matrix3d transNoise = Eigen::Matrix3d::Zero();
    for (int i = 0; i < 3; ++i)
      transNoise(i, i) = std::pow(noiseTranslation[i], 2);

    Eigen::Matrix3d rotNoise = Eigen::Matrix3d::Zero();
    for (int i = 0; i < 3; ++i)
      rotNoise(i, i) = std::pow(noiseRotation[i], 2);

    // Information matrix of the distribution
    Eigen::Matrix<double, 6, 6> information = Eigen::Matrix<double, 6, 6>::Zero();
    information.block<3,3>(0,0) = transNoise.inverse();
    information.block<3,3>(3,3) = rotNoise.inverse();

    // Gaussian noise generators
    transSampler.setDistribution(transNoise);
    rotSampler.setDistribution(rotNoise);

    if (randomSeed) {
      std::random_device r;
      std::seed_seq seedSeq{r(), r(), r(), r(), r()};
      vector<int> seeds(2);
      seedSeq.generate(seeds.begin(), seeds.end());
      transSampler.seed(seeds[0]);
      rotSampler.seed(seeds[1]);
    }
    return information;
}

void addNoiseToSubmap(GaussianGen& transSampler,
                      GaussianGen& rotSampler,
                      SubmapObj& submap){
    // Noise to rotation
    Vector3d quatXYZ = rotSampler.generateSample();
    double qw = 1.0 - quatXYZ.norm();
    if (qw < 0) {
      qw = 0.;
      cerr << "x";
    }
    Quaterniond rot(qw, quatXYZ.x(), quatXYZ.y(), quatXYZ.z());
    rot.normalize();

    // Noise to translation
    Vector3d trans = transSampler.generateSample();

    // Transform point cloud and submap frame
    Isometry3f poseDR(Isometry3f(Translation3f(trans.cast<float>())) *
                      Isometry3f(rot.cast<float>()));

    pcl::transformPointCloud(submap.submap_pcl_, submap.submap_pcl_, poseDR.matrix());
    submap.submap_tf_ = poseDR * submap.submap_tf_;
}

void additiveNoiseToSubmap(GaussianGen& transSampler,
                           GaussianGen& rotSampler,
                           SubmapObj& submap_i,
                           SubmapObj& submap_i_1){

    Vector3f euler_i = submap_i.submap_tf_.linear().eulerAngles(0,1,2);
    float roll = 0, pitch = 0, yaw = /*euler_i(0) +*/ 0.5;
    Matrix3f m;
    m = AngleAxisf(roll, Vector3f::UnitX())
        * AngleAxisf(pitch, Vector3f::UnitY())
        * AngleAxisf(yaw, Vector3f::UnitZ());

    Eigen::Vector3f rel_trans = submap_i.submap_tf_.translation() - submap_i_1.submap_tf_.translation();
    rel_trans.z() = 0;  // Fix Z coordinate

    Eigen::Vector3f p_i = submap_i_1.submap_tf_.translation() + m * rel_trans;

    // Transform point cloud and submap frame
    Vector3f rel = m * rel_trans - (rel_trans);

//    Isometry3f drift(Isometry3f(Translation3f(rel))*
//                     Isometry3f((m).matrix()));

    pcl::transformPointCloud(submap_i.submap_pcl_, submap_i.submap_pcl_,
                             (Isometry3f(Translation3f(rel))).matrix());

    submap_i.submap_tf_.translation() = p_i;
    submap_i.submap_tf_.linear() = (submap_i.submap_tf_.linear().matrix() * m).matrix();
}


void addNoiseToGraph(GraphConstructor& graph_obj){

    GaussianGen transSampler, rotSampler;
    Matrix<double, 6,6> information = generateGaussianNoise(transSampler, rotSampler);

    // noise for all the edges
//    for (size_t i = 0; i < graph_obj.edges_.size(); ++i) {
////      EdgeSE3* e = graph_obj.edges_[i];
//      Eigen::Isometry3d meas_i = graph_obj.lcMeas_.at(i);
//      Eigen::Quaterniond gtQuat = (Eigen::Quaterniond)meas_i.linear();
//      Eigen::Vector3d gtTrans = meas_i.translation();

//      Eigen::Vector3d quatXYZ = rotSampler.generateSample();
//      double qw = 1.0 - quatXYZ.norm();
//      if (qw < 0) {
//        qw = 0.;
//        cerr << "x";
//      }
//      Eigen::Quaterniond rot(qw, quatXYZ.x(), quatXYZ.y(), quatXYZ.z());
//      rot.normalize();
//      Eigen::Vector3d trans = transSampler.generateSample();
//      rot = gtQuat * rot;
//      trans = gtTrans + trans;

//      Eigen::Isometry3d noisyMeasurement = (Eigen::Isometry3d) rot;
//      noisyMeasurement.translation() = trans;
////      e->setMeasurement(noisyMeasurement);
//      graph_obj.lcMeas_.at(i) = noisyMeasurement;
//    }

    // noise for all the edges
    for (size_t i = 0; i < graph_obj.drEdges_.size(); ++i) {
      Eigen::Isometry3d meas_i = graph_obj.drMeas_.at(i);
      Eigen::Quaterniond gtQuat = (Eigen::Quaterniond)meas_i.linear();
      Eigen::Vector3d gtTrans = meas_i.translation();

      Eigen::Vector3d quatXYZ = rotSampler.generateSample();
      double qw = 1.0 - quatXYZ.norm();
      if (qw < 0) {
        qw = 0.;
        cerr << "x";
      }
      Eigen::Quaterniond rot(qw, quatXYZ.x(), quatXYZ.y(), quatXYZ.z());
      rot.normalize();
      Eigen::Vector3d trans = transSampler.generateSample();
      rot = gtQuat * rot;
      trans = gtTrans + trans;

//      std::cout << "DR " << i << " " << gtTrans.transpose() << std::endl;
//      std::cout << gtQuat.coeffs() << std::endl;

      Eigen::Isometry3d noisyMeasurement = (Eigen::Isometry3d) rot;
      noisyMeasurement.translation() = trans;
      graph_obj.drMeas_.at(i) = noisyMeasurement;

    }
}

