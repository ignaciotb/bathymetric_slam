#include "graph_optimization/types.h"

namespace ceres {
namespace optimizer {

std::istream& operator>>(std::istream& input, Pose3d& pose) {
    Eigen::Vector3d p;
    Eigen::Quaterniond q;

  input >> p.x() >> p.y() >> p.z() >> q.x() >>
      q.y() >> q.z() >> q.w();
  // Normalize the quaternion to account for precision loss due to
  // serialization.
  pose.p = p;
  pose.q = q.normalized().toRotationMatrix().eulerAngles(0,1,2);

  return input;
}

typedef std::map<int, Pose3d, std::less<int>,
                 Eigen::aligned_allocator<std::pair<const int, Pose3d> > >
    MapOfPoses;


std::istream& operator>>(std::istream& input, Constraint3d& constraint) {
  Pose3d& t_be = constraint.t_be;
  input >> constraint.id_begin >> constraint.id_end >> t_be;
//  t_be.q = constraint.t_be.q.toRotationMatrix().eulerAngles(0,1,2);

  for (int i = 0; i < 6 && input.good(); ++i) {
    for (int j = i; j < 6 && input.good(); ++j) {
      input >> constraint.information(i, j);
      if (i != j) {
        constraint.information(j, i) = constraint.information(i, j);
      }
    }
  }
  return input;
}

}  // namespace optimizer
}  // namespace ceres

