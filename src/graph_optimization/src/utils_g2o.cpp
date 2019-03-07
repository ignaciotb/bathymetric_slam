#include "graph_optimization/utils_g2o.hpp"

using namespace std;
using namespace g2o;
using namespace Eigen;

Matrix<double, 6,6> generateGaussianNoise(GaussianGen& transSampler,
                                          GaussianGen& rotSampler){

    bool randomSeed = false;
    std::vector<double> noiseTranslation;
    std::vector<double> noiseRotation;
    noiseTranslation.push_back(0.1);
    noiseTranslation.push_back(0.1);
    noiseTranslation.push_back(0.0001);
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
                      const Matrix<double, 6,6>& information,
                      SubmapObj& submap,
                      Isometry3f& poseDRt){
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
    poseDRt = poseDRt * poseDR;

    pcl::transformPointCloud(submap.submap_pcl_, submap.submap_pcl_, poseDRt.matrix());
    submap.submap_tf_ = poseDRt * submap.submap_tf_;

    // Store information matrix
    submap.submap_info_ = information;


}
