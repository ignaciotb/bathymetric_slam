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

    bool randomSeed = true;
    std::vector<double> noiseTranslation;
    std::vector<double> noiseRotation;
    noiseTranslation.push_back(3);
    noiseTranslation.push_back(3);
    noiseTranslation.push_back(0.001);
    noiseRotation.push_back(0.0001);
    noiseRotation.push_back(0.0001);
    noiseRotation.push_back(0.001);

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

    Eigen::Quaterniond gtQuat = (Eigen::Quaterniond)submap.submap_tf_.linear().cast<double>();
    Eigen::Vector3d gtTrans = submap.submap_tf_.translation().cast<double>();

    Eigen::Vector3d quatXYZ = rotSampler.generateSample();
    double qw = 1.0 - quatXYZ.norm();
    if (qw < 0) {
    qw = 0.;
    cerr << "x";
    }
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0,0.1};

//    Eigen::Quaterniond rot(qw, quatXYZ.x(), quatXYZ.y(), quatXYZ.z());
    // Bias in yaw
    double roll = 0.0, pitch = 0.0, yaw = /*0.001*/ d(gen);
    Matrix3d m;
    m = AngleAxisd(roll, Vector3d::UnitX())
        * AngleAxisd(pitch, Vector3d::UnitY())
        * AngleAxisd(yaw, Vector3d::UnitZ());
    Eigen::Quaterniond rot(m);

    Eigen::Vector3d trans;
//    trans = transSampler.generateSample();
    trans.setZero();

    trans = gtTrans + trans;
    rot = gtQuat * rot;

    Eigen::Isometry3d noisyMeasurement = (Eigen::Isometry3d) rot;
    noisyMeasurement.translation() = trans;

    // Transform submap_i pcl and tf
    pcl::transformPointCloud(submap.submap_pcl_, submap.submap_pcl_,
                             (noisyMeasurement.cast<float>() * submap.submap_tf_.inverse()).matrix());

    submap.submap_tf_ = noisyMeasurement.cast<float>();
}

void addNoiseToMap(GaussianGen& transSampler,
                   GaussianGen& rotSampler,
                   SubmapsVec& submap_set){

    for (size_t i =1; i < submap_set.size(); i++){

        std::cout << i << " -------" << std::endl;
        std::cout << submap_set.at(i-1).submap_tf_.matrix() << std::endl;

    }
    // Noise for all the submaps
    for (size_t i =1; i < submap_set.size(); i++){
        Eigen::Isometry3f tf_prev = submap_set.at(i-1).submap_tf_;

        Eigen::Isometry3f meas_i = tf_prev.inverse() * submap_set.at(i).submap_tf_;
        Eigen::Quaternionf gtQuat = (Eigen::Quaternionf)meas_i.linear();
        Eigen::Vector3f gtTrans = meas_i.translation();

        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<> d{0,0.5};

        // Bias in yaw
        float roll = 0.0, pitch = 0.0, yaw = /*0.5*/ d(gen);
        Matrix3f m;
        m = AngleAxisf(roll, Vector3f::UnitX())
            * AngleAxisf(pitch, Vector3f::UnitY())
            * AngleAxisf(yaw, Vector3f::UnitZ());
        Eigen::Quaternionf rot(m);

        rot = gtQuat * rot;

        meas_i = (Eigen::Isometry3f) rot;
        meas_i.translation() = gtTrans;

        Eigen::Isometry3f estimate_i = tf_prev * meas_i;

        std::cout << i << " -------" << std::endl;
        std::cout << tf_prev.matrix() << std::endl;
        std::cout << meas_i.matrix() << std::endl;
        std::cout << estimate_i.matrix() << std::endl;

        // Transform submap_i pcl and tf
        pcl::transformPointCloud(submap_set.at(i).submap_pcl_, submap_set.at(i).submap_pcl_,
                                 (estimate_i * submap_set.at(i).submap_tf_.inverse()).matrix());

        submap_set.at(i).submap_tf_ = estimate_i.cast<float>();
    }
}
