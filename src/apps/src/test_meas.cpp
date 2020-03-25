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

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <cereal/archives/binary.hpp>

#include "data_tools/std_data.h"
//#include "data_tools/benchmark.h"

#include "submaps_tools/cxxopts.hpp"
#include "submaps_tools/submaps.hpp"

#include "registration/utils_visualization.hpp"
#include "registration/gicp_reg.hpp"
#include "meas_models/multibeam_simple.hpp"

#define INTERACTIVE 0
#define VISUAL 1

using namespace Eigen;
using namespace std;
using namespace g2o;

int main(int argc, char** argv){

    // Inputs
    std::string track_str, map_str, output_str, original, simulation;
    cxxopts::Options options("MyProgram", "One line description of MyProgram");
    options.add_options()
        ("help", "Print help")
        ("output_cereal", "Output graph cereal", cxxopts::value(output_str))
        ("original", "Disturb original trajectory", cxxopts::value(original))
        ("simulation", "Simulation data from Gazebo", cxxopts::value(simulation))
        ("trajectory", "Input AUV GT data", cxxopts::value(track_str))
        ("map", "Localization map", cxxopts::value(map_str));

    auto result = options.parse(argc, argv);
    if (result.count("help")) {
        cout << options.help({ "", "Group" }) << endl;
        exit(0);
    }
    if(output_str.empty()){
        output_str = "output_cereal.cereal";
    }
    boost::filesystem::path output_path(output_str);
    string outFilename = "graph_corrupted.g2o";   // G2O output file

    // Parse submaps from cereal file
    SubmapsVec map_gt, traj_gt;
    boost::filesystem::path map_path(map_str);
    boost::filesystem::path auv_path(track_str);
    std::cout << "Map path " << boost::filesystem::basename(map_path) << std::endl;
    std::cout << "AUV path " << boost::filesystem::basename(auv_path) << std::endl;
    if(simulation == "yes"){
        map_gt = readSubmapsInDir(map_path.string());
    }
    else{
        if(original == "yes"){
            MapObj map_loc;
            Eigen::Isometry3d map_tf;
            std_data::pt_submaps ss = std_data::read_data<std_data::pt_submaps>(map_path);
            std::tie(map_loc, map_tf)= parseMapAUVlib(ss);
            map_gt.push_back(map_loc);

            std_data::mbes_ping::PingsT std_pings = std_data::read_data<std_data::mbes_ping::PingsT>(auv_path);
            std::cout << "Number of pings " << std_pings.size() << std::endl;
            traj_gt = parsePingsAUVlib(std_pings);
        }
        else{
            std::ifstream is(boost::filesystem::basename(map_path) + ".cereal", std::ifstream::binary);
            {
              cereal::BinaryInputArchive iarchive(is);
              iarchive(map_gt);
            }
        }
        // Filtering of submaps
        PointCloudT::Ptr cloud_ptr (new PointCloudT);
        pcl::UniformSampling<PointT> us_filter;
        us_filter.setInputCloud (cloud_ptr);
        us_filter.setRadiusSearch(1);   // 1 for Borno, 2 for Antarctica
        for(SubmapObj& submap_i: map_gt){
    //        std::cout << "before " << submap_i.submap_pcl_.size() << std::endl;
            *cloud_ptr = submap_i.submap_pcl_;
            us_filter.setInputCloud(cloud_ptr);
            us_filter.filter(*cloud_ptr);
            submap_i.submap_pcl_ = *cloud_ptr;
    //        std::cout << submap_i.submap_pcl_.size() << std::endl;
        }
    }

    // Visualization of initial map
    PCLVisualizer viewer ("Map viewer");
    SubmapsVisualizer* visualizer = new SubmapsVisualizer(viewer);
    visualizer->setVisualizer(map_gt, 2);
    while(!viewer.wasStopped ()){
        viewer.spinOnce ();
    }
    viewer.resetStoppedFlag();

    // Create voxel grid
    MultibeamSensor<PointT> vox_oc;
    vox_oc.setLeafSize(1,1,1);
    vox_oc.initializeVoxelGrid(map_gt.at(0));

    // Plot localization trajectory
    std::cout << "Number of pings " << traj_gt.size() << std::endl;
//    std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i> > occluded_voxels;
//    std::vector<int> idxs;
    PointCloudT ping_i;
    for(int j=0; j<traj_gt.size(); j=j+50){
//    for(SubmapObj& submap_i: traj_gt){
        // Add GT ping
        SubmapObj submap_i = traj_gt.at(j);
        submap_i.colors_ = Eigen::Vector3d(255,0, 0);
        visualizer->addSubmap(submap_i, 2);

//        visualizer->plotMBESPing(submap_i, 2.0944, 256, 2);

        // Compute simulated ping
        vox_oc.createMBES(1.74/2, 3, submap_i.submap_tf_);
        ping_i.clear();
        vox_oc.pingComputation(ping_i);
//        PointCloudT pcl_filtered = vox_oc.getFilteredPointCloud();

        SubmapObj simulated_ping;
        simulated_ping.submap_tf_ = submap_i.submap_tf_;
        simulated_ping.submap_pcl_ = ping_i;
        std::cout << "Voxels hit: " << simulated_ping.submap_pcl_.size() << std::endl;
        simulated_ping.colors_ = Eigen::Vector3d(0,0,255);
        visualizer->addSubmap(simulated_ping, 2);

        while(!viewer.wasStopped ()){
            viewer.spinOnce ();
        }
        viewer.resetStoppedFlag();
    }

    return 0;
}
