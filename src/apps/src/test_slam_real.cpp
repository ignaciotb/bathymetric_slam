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
#include "data_tools/benchmark.h"

#include "submaps_tools/cxxopts.hpp"
#include "submaps_tools/submaps.hpp"

#include "registration/utils_visualization.hpp"
#include "registration/gicp_reg.hpp"

#include "graph_optimization/utils_g2o.hpp"
#include "graph_optimization/graph_construction.hpp"
#include "graph_optimization/ceres_optimizer.hpp"
#include "graph_optimization/read_g2o.h"

#define INTERACTIVE 0
#define VISUAL 0

using namespace Eigen;
using namespace std;
using namespace g2o;

int main(int argc, char** argv){

    // Inputs
    std::string folder_str, path_str, output_str, original;
    cxxopts::Options options("MyProgram", "One line description of MyProgram");
    options.add_options()
        ("help", "Print help")
        ("covs_folder", "Input covs folder", cxxopts::value(folder_str))
        ("output_cereal", "Output graph cereal", cxxopts::value(output_str))
        ("original", "Disturb original trajectory", cxxopts::value(original))
        ("slam_cereal", "Input ceres file", cxxopts::value(path_str));

    auto result = options.parse(argc, argv);
    if (result.count("help")) {
        cout << options.help({ "", "Group" }) << endl;
        exit(0);
    }
    if(output_str.empty()){
        output_str = "output_cereal.cereal";
    }
    boost::filesystem::path output_path(output_str);
    string loopsFilename = "loop_closures.txt"; // LCs list output file
    string outFilename = "graph_corrupted.g2o";   // G2O output file

    // Parse submaps from cereal file
    SubmapsVec submaps_gt;
    boost::filesystem::path submaps_path(path_str);
    std::cout << "Input data " << boost::filesystem::basename(submaps_path) << std::endl;
    if(original == "yes"){
        std_data::pt_submaps ss = std_data::read_data<std_data::pt_submaps>(submaps_path);
        submaps_gt = parseSubmapsAUVlib(ss);
    }
    else{
        std::ifstream is(boost::filesystem::basename(submaps_path) + ".cereal", std::ifstream::binary);
        {
          cereal::BinaryInputArchive iarchive(is);
          iarchive(submaps_gt);
        }
    }

    // Read training covs from folder
    covs covs_lc;
    boost::filesystem::path folder(folder_str);
    if(boost::filesystem::is_directory(folder)) {
        covs_lc = readCovsFromFiles(folder);
    }

    // Filtering of submaps
    PointCloudT::Ptr cloud_ptr (new PointCloudT);
    pcl::UniformSampling<PointT> us_filter;
    us_filter.setInputCloud (cloud_ptr);
    us_filter.setRadiusSearch(1);   // 1 for Borno, 2 for Antarctica
    for(SubmapObj& submap_i: submaps_gt){
//        std::cout << "before " << submap_i.submap_pcl_.size() << std::endl;
        *cloud_ptr = submap_i.submap_pcl_;
        us_filter.setInputCloud(cloud_ptr);
        us_filter.filter(*cloud_ptr);
        submap_i.submap_pcl_ = *cloud_ptr;
//        std::cout << submap_i.submap_pcl_.size() << std::endl;
    }

    // Benchmark GT
    benchmark::track_error_benchmark benchmark("real_data");
    PointsT gt_map = pclToMatrixSubmap(submaps_gt);
    PointsT gt_track = trackToMatrixSubmap(submaps_gt);
    benchmark.add_ground_truth(gt_map, gt_track);
    ceres::optimizer::saveOriginalTrajectory(submaps_gt); // Save original trajectory to txt

    // Visualization
#if VISUAL == 1
    PCLVisualizer viewer ("Submaps viewer");
    viewer.loadCameraParameters("Antarctica7");
    SubmapsVisualizer* visualizer = new SubmapsVisualizer(viewer);
    visualizer->setVisualizer(submaps_gt, 1);
    while(!viewer.wasStopped ()){
        viewer.spinOnce ();
    }
    viewer.resetStoppedFlag();
//    viewer.saveCameraParameters("Antarctica7");
#endif

    // GICP reg for submaps
    SubmapRegistration* gicp_reg = new SubmapRegistration();

    // Graph constructor
    GraphConstructor* graph_obj = new GraphConstructor(covs_lc);
    ofstream fileOutputStream;
    fileOutputStream.open(loopsFilename, std::ofstream::out);

    // Noise generators
    GaussianGen transSampler, rotSampler;
    Matrix<double, 6,6> information = generateGaussianNoise(transSampler, rotSampler);

    SubmapObj submap_trg;
    SubmapsVec submaps_prev, submaps_reg;
    double info_thres = 0.1;
    for(SubmapObj& submap_i: submaps_gt){
        std::cout << " ----------- Submap " << submap_i.submap_id_ << ", swath "
                  << submap_i.swath_id_ << " ------------"<< std::endl;

        // Look for loop closures
        for(SubmapObj& submap_k: submaps_reg){
            // Don't look for overlaps between submaps of the same swath or the prev submap
            if(submap_k.submap_id_ != submap_i.submap_id_ - 1){
                submaps_prev.push_back(submap_k);
            }
        }
        submap_i.findOverlaps(submaps_prev);
        submaps_prev.clear();
        submaps_reg.push_back(submap_i); // Add submap_i to registered set (just for visualization here)

#if INTERACTIVE == 1
        // Update visualizer
        visualizer->updateVisualizer(submaps_reg);
        while(!viewer.wasStopped ()){
            viewer.spinOnce ();
        }
        viewer.resetStoppedFlag();
#endif
        // Create graph vertex i
        graph_obj->createNewVertex(submap_i);

        // Create DR edge i and store (skip submap 0)
        if(submap_i.submap_id_ != 0 ){
            std::cout << "DR edge from " << submap_i.submap_id_ -1 << " to " << submap_i.submap_id_<< std::endl;
            graph_obj->createDREdge(submap_i);
        }

        // If potential loop closure detected
        SubmapObj submap_final = submap_i;
        if(!submap_i.overlaps_idx_.empty()){
            // Save loop closure to txt
            if(fileOutputStream.is_open()){
                fileOutputStream << submap_i.submap_id_;
                for(unsigned int j=0; j<submap_i.overlaps_idx_.size(); j++){
                    fileOutputStream << " " << submap_i.overlaps_idx_.at(j);
                }
                fileOutputStream << "\n";
            }

            // Register overlapping submaps
            submap_trg = gicp_reg->constructTrgSubmap(submaps_reg, submap_i.overlaps_idx_);
            addNoiseToSubmap(transSampler, rotSampler, submap_i); // Add disturbance to source submap
            if(gicp_reg->gicpSubmapRegistration(submap_trg, submap_i)){
                submap_final = submap_i;
            }
            submap_trg.submap_pcl_.clear();

            // Create loop closures
            graph_obj->findLoopClosures(submap_final, submaps_reg, info_thres);
        }
        submaps_reg.pop_back(); // Remove unregistered submap_i
        submaps_reg.push_back(submap_final);    // Add registered submap_i

#if INTERACTIVE == 1
        // Update visualizer
        visualizer->updateVisualizer(submaps_reg);
        while(!viewer.wasStopped ()){
            viewer.spinOnce ();
        }
        viewer.resetStoppedFlag();
#endif
    }
    fileOutputStream.close();

#if VISUAL == 1
    // Update visualizer
    visualizer->updateVisualizer(submaps_reg);
    while(!viewer.wasStopped ()){
        viewer.spinOnce ();
    }
    viewer.resetStoppedFlag();
#endif

    // Add noise to edges on the graph
    addNoiseToGraph(transSampler, rotSampler, *graph_obj);

    // Create initial DR chain and visualize
    graph_obj->createInitialEstimate(submaps_reg);

#if VISUAL == 1
    visualizer->plotPoseGraphG2O(*graph_obj, submaps_reg);
    while(!viewer.wasStopped ()){
        viewer.spinOnce ();
    }
    viewer.resetStoppedFlag();
#endif

    // Save graph to output g2o file (optimization can be run with G2O)
    graph_obj->saveG2OFile(outFilename);

    // Benchmar corrupted
    PointsT reg_map = pclToMatrixSubmap(submaps_reg);
    PointsT reg_track = trackToMatrixSubmap(submaps_reg);
    benchmark.add_benchmark(reg_map, reg_track, "corrupted");

    // Optimize graph and save to cereal
    google::InitGoogleLogging(argv[0]);
    ceres::optimizer::MapOfPoses poses = ceres::optimizer::ceresSolver(outFilename, graph_obj->drEdges_.size());
    ceres::optimizer::updateSubmapsCeres(poses, submaps_reg);
    std::cout << "Output cereal: " << boost::filesystem::basename(output_path) << std::endl;
    std::ofstream os(boost::filesystem::basename(output_path) + ".cereal", std::ofstream::binary);
    {
        cereal::BinaryOutputArchive oarchive(os);
        oarchive(submaps_reg);
        os.close();
    }

#if VISUAL == 1
    // Visualize Ceres output
    visualizer->plotPoseGraphCeres(submaps_reg);
    while(!viewer.wasStopped ()){
        viewer.spinOnce ();
    }
    viewer.resetStoppedFlag();
    delete(visualizer);
#endif

    // Benchmark Optimized
    PointsT opt_map = pclToMatrixSubmap(submaps_reg);
    PointsT opt_track = trackToMatrixSubmap(submaps_reg);
    benchmark.add_benchmark(opt_map, opt_track, "optimized");
    benchmark.print_summary();

    std::string command_str = "./plot_results.py --initial_poses poses_original.txt --corrupted_poses poses_corrupted.txt --optimized_poses poses_optimized.txt";
    const char *command = command_str.c_str();
    system(command);

    delete(gicp_reg);
    delete(graph_obj);

    return 0;
}
