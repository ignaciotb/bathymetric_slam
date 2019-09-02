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

#ifndef CERES_POSE_GRAPH_3D_ERROR_TERM_H_
#define CERES_POSE_GRAPH_3D_ERROR_TERM_H_

#include "Eigen/Core"
#include "ceres/autodiff_cost_function.h"

#include "graph_optimization/types.h"

namespace ceres {
namespace optimizer {

template <typename T> Eigen::Quaternion<T> eulerToQuat(const Eigen::Matrix<T,3,1>& vec){

    Eigen::Quaternion<T> q = Eigen::AngleAxis<T>(vec[0], Eigen::Matrix<T, 3, 1>::UnitX())
                            * Eigen::AngleAxis<T>(vec[1], Eigen::Matrix<T, 3, 1>::UnitY())
                            * Eigen::AngleAxis<T>(vec[2], Eigen::Matrix<T, 3, 1>::UnitZ());
    return q;
}


class PoseGraph3dErrorTerm {
 public:
  PoseGraph3dErrorTerm(const Pose3d& t_ab_measured, const Eigen::Matrix<double, 6, 6>& sqrt_information)
      : t_ab_measured_(t_ab_measured), sqrt_information_(sqrt_information) {}

  template <typename T>
  bool operator()(const T* const p_a_ptr, const T* const q_a_ptr,
                  const T* const p_b_ptr, const T* const q_b_ptr,
                  T* residuals_ptr) const
  {
    Eigen::Map<const Eigen::Matrix<T, 3, 1> > p_a(p_a_ptr);
    Eigen::Map<const Eigen::Matrix<T, 3, 1> > q_a(q_a_ptr);
    Eigen::Map<const Eigen::Matrix<T, 3, 1> > p_b(p_b_ptr);
    Eigen::Map<const Eigen::Matrix<T, 3, 1> > q_b(q_b_ptr);

    // Compute the relative transformation between the two frames.
    Eigen::Quaternion<T> q_a_quat  = eulerToQuat<T>(q_a);
    Eigen::Quaternion<T> q_b_quat  = eulerToQuat<T>(q_b);
    Eigen::Quaternion<T> q_a_inverse = q_a_quat.conjugate();
    Eigen::Quaternion<T> q_ab_estimated = q_a_inverse * q_b_quat;

    // Represent the displacement between the two frames in the A frame.
    Eigen::Matrix<T, 3, 1> p_ab_estimated = q_a_inverse * (p_b - p_a);

    // Compute the error between the two orientation estimates.
    Eigen::Matrix<T, 3, 1> t_ab_meas_q = t_ab_measured_.q.template cast<T>();
    Eigen::Quaternion<T> q_meas = eulerToQuat(t_ab_meas_q);
    Eigen::Quaternion<T> delta_q = q_meas * q_ab_estimated.conjugate();

    // Compute the residuals.
    Eigen::Map<Eigen::Matrix<T, 6, 1> > residuals(residuals_ptr);
    residuals.template block<3, 1>(0, 0) = p_ab_estimated - t_ab_measured_.p.template cast<T>();
    residuals.template block<3, 1>(3, 0) = T(2.0) * delta_q.vec();

    // Scale the residuals by the measurement uncertainty.
    residuals.applyOnTheLeft(sqrt_information_.template cast<T>());

    return true;
  }

  static ceres::CostFunction* Create(const Pose3d& t_ab_measured, const Eigen::Matrix<double, 6, 6>& sqrt_information)
  {
    return new ceres::AutoDiffCostFunction<PoseGraph3dErrorTerm, 6, 3, 3, 3, 3>(
        new PoseGraph3dErrorTerm(t_ab_measured, sqrt_information));
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

 private:
  // The measurement for the position of B relative to A in the A frame.
  const Pose3d t_ab_measured_;
  // The square root of the measurement information matrix.
  const Eigen::Matrix<double, 6, 6> sqrt_information_;
};

}  // namespace optimizer
}  // namespace ceres

#endif  // CERES_POSE_GRAPH_3D_ERROR_TERM_H_
