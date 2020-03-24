#ifndef MULTIBEAM_SIMPLE_HPP
#define MULTIBEAM_SIMPLE_HPP

#include <pcl/common/common.h>
#include <pcl/filters/voxel_grid.h>
#include "submaps_tools/submaps.hpp"

template <typename PointT>
class MultibeamSensor: public pcl::VoxelGrid<PointT>
{
protected:

    using pcl::VoxelGrid<PointT>::min_b_;
    using pcl::VoxelGrid<PointT>::max_b_;
    using pcl::VoxelGrid<PointT>::div_b_;
    using pcl::VoxelGrid<PointT>::leaf_size_;
    using pcl::VoxelGrid<PointT>::inverse_leaf_size_;

    using PointCloudFilt = typename pcl::Filter<PointT>::PointCloud;
    using PointCloudFiltPtr = typename PointCloudFilt::Ptr;
    using PointCloudFiltConstPtr = typename PointCloudFilt::ConstPtr;

public:
  /** \brief Empty constructor. */
  MultibeamSensor ()
  {
    initialized_ = false;
    this->setSaveLeafLayout (true);
  }

  /** \brief Destructor. */
  ~MultibeamSensor ()
  {
  }

  inline PointCloudFilt
  getFilteredPointCloud () { return filtered_cloud_; }


  /** \brief Returns the minimum bounding of coordinates of the voxel grid (x,y,z).
    * \return the minimum coordinates (x,y,z)
    */
  inline Eigen::Vector3f
  getMinBoundCoordinates () { return (b_min_.head<3> ()); }

  /** \brief Returns the maximum bounding of coordinates of the voxel grid (x,y,z).
    * \return the maximum coordinates (x,y,z)
    */
  inline Eigen::Vector3f
  getMaxBoundCoordinates () { return (b_max_.head<3> ()); }

  /** \brief Returns the corresponding centroid (x,y,z) coordinates
    * in the grid of voxel (i,j,k).
    * \param[in] ijk the coordinate (i, j, k) of the voxel
    * \return the (x,y,z) coordinate of the voxel centroid
    */
  inline Eigen::Vector4f
  getCentroidCoordinate (const Eigen::Vector3i& ijk)
  {
    int i,j,k;
    i = ((b_min_[0] < 0) ? (std::abs (min_b_[0]) + ijk[0]) : (ijk[0] - min_b_[0]));
    j = ((b_min_[1] < 0) ? (std::abs (min_b_[1]) + ijk[1]) : (ijk[1] - min_b_[1]));
    k = ((b_min_[2] < 0) ? (std::abs (min_b_[2]) + ijk[2]) : (ijk[2] - min_b_[2]));

    Eigen::Vector4f xyz;
    xyz[0] = b_min_[0] + (leaf_size_[0] * 0.5f) + (static_cast<float> (i) * leaf_size_[0]);
    xyz[1] = b_min_[1] + (leaf_size_[1] * 0.5f) + (static_cast<float> (j) * leaf_size_[1]);
    xyz[2] = b_min_[2] + (leaf_size_[2] * 0.5f) + (static_cast<float> (k) * leaf_size_[2]);
    xyz[3] = 0;
    return xyz;
  }

  inline float
  round (float d)
  {
    return static_cast<float> (std::floor (d + 0.5f));
  }

  // We use round here instead of std::floor due to some numerical issues.
  /** \brief Returns the corresponding (i,j,k) coordinates in the grid of point (x,y,z).
    * \param[in] x the X point coordinate to get the (i, j, k) index at
    * \param[in] y the Y point coordinate to get the (i, j, k) index at
    * \param[in] z the Z point coordinate to get the (i, j, k) index at
    */
  inline Eigen::Vector3i
  getGridCoordinatesRound (float x, float y, float z)
  {
    return Eigen::Vector3i (static_cast<int> (round (x * inverse_leaf_size_[0])),
                            static_cast<int> (round (y * inverse_leaf_size_[1])),
                            static_cast<int> (round (z * inverse_leaf_size_[2])));
  }

  // initialization flag
  bool initialized_;

  Eigen::Vector4f sensor_origin_;
  Eigen::Quaternionf sensor_orientation_;
  std::vector<Eigen::Quaternionf, Eigen::aligned_allocator<Eigen::Quaternionf> > beam_directions_;

  // minimum and maximum bounding box coordinates
  Eigen::Vector4f b_min_, b_max_;

  // voxel grid filtered cloud
  PointCloudFilt filtered_cloud_;

  void initializeVoxelGrid (SubmapObj submap_i){
      std::cout << "Number of points " << submap_i.submap_pcl_.points.size() << std::endl;
      // initialization set to true
      initialized_ = true;
      // create the voxel grid and store the output cloud
      PointCloudT::Ptr cloud_ptr(new PointCloudT);
      *cloud_ptr = submap_i.submap_pcl_;
      this->setInputCloud(cloud_ptr);
      this->filter (filtered_cloud_);
//      PointCloudT filtered = filtered_cloud_;

      // Get the minimum and maximum bounding box dimensions
      b_min_[0] = (static_cast<float> ( min_b_[0]) * leaf_size_[0]);
      b_min_[1] = (static_cast<float> ( min_b_[1]) * leaf_size_[1]);
      b_min_[2] = (static_cast<float> ( min_b_[2]) * leaf_size_[2]);
      b_max_[0] = (static_cast<float> ( (max_b_[0]) + 1) * leaf_size_[0]);
      b_max_[1] = (static_cast<float> ( (max_b_[1]) + 1) * leaf_size_[1]);
      b_max_[2] = (static_cast<float> ( (max_b_[2]) + 1) * leaf_size_[2]);
      std::cout << "Minb " << min_b_.transpose() << std::endl;
      std::cout << "Maxb " << max_b_.transpose() << std::endl;
  }


  void createMBES(float spam, float n_beams, Eigen::Isometry3f& sensor_tf){

      // set the sensor origin and sensor orientation
      sensor_origin_ << sensor_tf.translation(), 0.0;
      sensor_orientation_ = Eigen::Quaternionf(sensor_tf.linear());


      float roll_step = spam/n_beams;
      float pitch = 0.0, yaw = 0.0;
      Eigen::Quaternionf q_180;
      q_180 = Eigen::AngleAxisf(3.1415, Eigen::Vector3f::UnitX())
          * Eigen::AngleAxisf(pitch, Eigen::Vector3f::UnitY())
          * Eigen::AngleAxisf(yaw, Eigen::Vector3f::UnitZ());

      sensor_orientation_ *= q_180;
      Eigen::Vector3f z_or = sensor_orientation_.toRotationMatrix().col(2);
//      std::cout << "Frame origin " << sensor_origin_.transpose() << std::endl;
//      std::cout << "Frame direction " << z_or.transpose() << std::endl;

      for(int i = -n_beams/2; i<=n_beams/2; i++){
          Eigen::Quaternionf q = Eigen::AngleAxisf(roll_step*i, Eigen::Vector3f::UnitX())
                                  * Eigen::AngleAxisf(pitch, Eigen::Vector3f::UnitY())
                                  * Eigen::AngleAxisf(yaw, Eigen::Vector3f::UnitZ());
          beam_directions_.push_back(Eigen::Quaternionf(sensor_orientation_ * q));
      }

  }

  int pingComputation (PointCloudT& ping_pcl){
    if (!initialized_){
        PCL_ERROR ("Voxel grid not initialized; call initializeVoxelGrid () first! \n");
        return -1;
    }
    Eigen::Vector4f direction;

    // Check every beam
    for(Eigen::Quaternionf beam_n: beam_directions_){
        Eigen::Matrix3f rot_beam = beam_n.toRotationMatrix();
        direction << rot_beam.col(2), 0.0;
        direction.normalize ();
        // Estimate entry point into the voxel grid
        float tmin = rayBoxIntersection(sensor_origin_, direction);
        if (tmin != -1){
            float z_max = 300;
            for(float z=0; z<z_max; z+=1){
                // coordinate of the boundary of the voxel grid
                Eigen::Vector4f start = sensor_origin_ + (tmin+z) * direction;
                // i,j,k coordinate of the voxel were the ray enters the voxel grid
                Eigen::Vector3i ijk = getGridCoordinatesRound (start[0], start[1], start[2]);
                // centroid coordinate of the entry voxel
                Eigen::Vector4f voxel_max = getCentroidCoordinate (ijk);
                // if voxel is occluded
                int index = this->getCentroidIndexAt (ijk);
                if (index != -1){
                    ping_pcl.points.push_back(PointT(voxel_max[0], voxel_max[1], voxel_max[2]));
                    break;
                }
            }
        }
    }
    beam_directions_.clear();
    return 0;

  }

float rayBoxIntersection (const Eigen::Vector4f& origin, const Eigen::Vector4f& direction){

    float tmin, tmax, tymin, tymax, tzmin, tzmax;
    if (direction[0] >= 0){
        tmin = (b_min_[0] - origin[0]) / direction[0];
        tmax = (b_max_[0] - origin[0]) / direction[0];
    }
    else{
        tmin = (b_max_[0] - origin[0]) / direction[0];
        tmax = (b_min_[0] - origin[0]) / direction[0];
    }


    if (direction[1] >= 0){
        tymin = (b_min_[1] - origin[1]) / direction[1];
        tymax = (b_max_[1] - origin[1]) / direction[1];
    }
    else{
        tymin = (b_max_[1] - origin[1]) / direction[1];
        tymax = (b_min_[1] - origin[1]) / direction[1];
    }

    if ((tmin > tymax) || (tymin > tmax)){
        // PCL_ERROR ("no intersection with the bounding box \n");
        tmin = -1.0f;
        return tmin;
    }

    if (tymin > tmin)
        tmin = tymin;
    if (tymax < tmax)
        tmax = tymax;

    if (direction[2] >= 0){
        tzmin = (b_min_[2] - origin[2]) / direction[2];
        tzmax = (b_max_[2] - origin[2]) / direction[2];
    }
    else{
        tzmin = (b_max_[2] - origin[2]) / direction[2];
        tzmax = (b_min_[2] - origin[2]) / direction[2];
    }

    if ((tmin > tzmax) || (tzmin > tmax)){
        // PCL_ERROR ("no intersection with the bounding box \n");
        tmin = -1.0f;
        return tmin;
    }

    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;

    return tmin;
 }
};

#endif // MULTIBEAM_SIMPLE_HPP
