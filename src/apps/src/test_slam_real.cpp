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
#include <boost/filesystem.hpp>

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "submaps_tools/cxxopts.hpp"

#include "data_tools/std_data.h"
#include "data_tools/benchmark.h"

#include "submaps_tools/submaps.hpp"

#include "registration/utils_visualization.hpp"
#include "registration/gicp_reg.hpp"

#include "graph_optimization/utils_g2o.hpp"
#include "graph_optimization/graph_construction.hpp"
#include "graph_optimization/ceres_optimizer.hpp"
#include "graph_optimization/read_g2o.h"

#include <boost/algorithm/string.hpp>

#define INTERACTIVE 0
#define NOISE 0

using namespace Eigen;
using namespace std;
using namespace g2o;

/// Read loop closures from external file
std::vector<std::vector<int> > readGTLoopClosures(string& fileName, int submaps_nb){

    std::ifstream infile(fileName);
    std::vector<std::vector<int> > overlaps(static_cast<unsigned long>(submaps_nb));

    std::string line;
    while(std::getline(infile, line)){
        vector<string> result;
        boost::split(result, line, boost::is_any_of(" "));
        vector<int> result_d(result.size());
        std::transform(result.begin(), result.end(), result_d.begin(), [](const std::string& val){
            return std::stoi(val);
        });
        overlaps.at(result_d[0]).insert(overlaps.at(result_d[0]).end(), result_d.begin()+1, result_d.end());
    }

    return overlaps;
}

int main(int argc, char** argv){

    // Read submaps from folder
    string path_str;
    cxxopts::Options options("MyProgram", "One line description of MyProgram");
    options.add_options()
        ("help", "Print help")
        ("file", "Input ceres file", cxxopts::value(path_str));

    auto result = options.parse(argc, argv);
    if (result.count("help")) {
        cout << options.help({ "", "Group" }) << endl;
        exit(0);
    }
    boost::filesystem::path submaps_path(path_str);// / "submaps_full.cereal";
    std_data::pt_submaps ss = std_data::read_data<std_data::pt_submaps>(submaps_path);
    string loopsFilename = "loop_closures.txt";

    // Parse data structures
    SubmapsVec submaps_gt = parseSubmapsAUVlib(ss);

    // Benchmark GT
    benchmark::track_error_benchmark benchmark("real_data");
    PointsT gt_map = pclToMatrixSubmap(submaps_gt);
    PointsT gt_track = trackToMatrixSubmap(submaps_gt);
    benchmark.add_ground_truth(gt_map, gt_track);

    // Parameters for downsampling and filtering of submaps
    PointCloudT::Ptr cloud_ptr (new PointCloudT);
    pcl::UniformSampling<PointT> us_filter;
    us_filter.setInputCloud (cloud_ptr);
    us_filter.setRadiusSearch(0.2);
    pcl::StatisticalOutlierRemoval<pcl::PointXYZ> sor_filter;
    sor_filter.setMeanK (100);
    sor_filter.setStddevMulThresh(3);

    for(SubmapObj& submap_i: submaps_gt){
        std::cout << "Before preprocessing: " << submap_i.submap_pcl_.points.size() << std::endl;
        *cloud_ptr = submap_i.submap_pcl_;
        us_filter.setInputCloud(cloud_ptr);
        us_filter.filter(*cloud_ptr);
//        sor_filter.setInputCloud (cloud_ptr);
//        sor_filter.filter (*cloud_ptr);
        submap_i.submap_pcl_ = *cloud_ptr;
        std::cout << "After preprocessing: " << submap_i.submap_pcl_.points.size() << std::endl;
    }


    // Visualization
//    PCLVisualizer viewer ("Submaps viewer");
//    viewer.loadCameraParameters("Antarctica7");
//    SubmapsVisualizer* visualizer = new SubmapsVisualizer(viewer);
//    visualizer->setVisualizer(submaps_gt, 2);
//    while(!viewer.wasStopped ()){
//        viewer.spinOnce ();
//    }
//    viewer.resetStoppedFlag();
//    viewer.saveCameraParameters("Antarctica7");


#if NOISE == 1
    // Add Gaussian noise to vehicle's pose among submaps
    GaussianGen transSampler, rotSampler;
    Matrix<double, 6,6> information = generateGaussianNoise(transSampler, rotSampler);
    for(SubmapObj& submap_i: submaps_gt){
        addNoiseToSubmap(transSampler, rotSampler, submap_i);
//        additiveNoiseToSubmap(transSampler, rotSampler, submaps_gt.at(i), submaps_gt.at(i-1));
    }

    // Benchmark noisy
    PointsT corrupt_map = pclToMatrixSubmap(submaps_gt);
    PointsT corrupt_track = trackToMatrixSubmap(submaps_gt);
    benchmark.add_benchmark(corrupt_map, corrupt_track, "corrupted");

    // Get GT loop closures from external file
#endif
    std::vector<std::vector<int> > vecLCs = readGTLoopClosures(loopsFilename, submaps_gt.size());

    // GICP reg for submaps
    SubmapRegistration* gicp_reg = new SubmapRegistration();

    // Graph constructor
    GraphConstructor* graph_obj = new GraphConstructor();

    ofstream fileOutputStream;
    fileOutputStream.open(loopsFilename, std::ofstream::out);

    SubmapObj submap_trg;
    SubmapsVec submaps_prev;
    double info_thres = 0.1;
    SubmapsVec submaps_reg;
    for(SubmapObj& submap_i: submaps_gt){
        std::cout << " ----------- Submap " << submap_i.submap_id_ << ", swath "
                  << submap_i.swath_id_ << " ------------"<< std::endl;

        // Look for loop closures
#if NOISE == 0
//        for(SubmapObj& submap_k: submaps_reg){
//            // Don't look for overlaps between submaps of the same swath or the prev submap
//            if(submap_k.swath_id_ != submap_i.swath_id_ ||
//                    submap_k.submap_id_ != submap_i.submap_id_ - 1){
//                submaps_prev.push_back(submap_k);
//            }
//        }
//        submap_i.findOverlaps(submaps_prev);
//        submaps_prev.clear();
#else
#endif
        submap_i.overlaps_idx_ = vecLCs.at(submap_i.submap_id_);

        // Add submap_i to registered set (just for visualization here)
        submaps_reg.push_back(submap_i);

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
            graph_obj->createDREdge(submap_i);
        }

        // If potential loop closure detected
        SubmapObj submap_final = submap_i;
        if(!submap_i.overlaps_idx_.empty()){
            std::cout << "Overlaps of submap " << submap_i.submap_id_ << std::endl;
            if(fileOutputStream.is_open()){
                fileOutputStream << submap_i.submap_id_;
                for(unsigned int j=0; j<submap_i.overlaps_idx_.size(); j++){
                    fileOutputStream << " " << submap_i.overlaps_idx_.at(j);
                    std::cout << submap_i.overlaps_idx_.at(j) << std::endl;
                }
                fileOutputStream << "\n";
            }

            // Register overlapping submaps
            submap_trg = gicp_reg->constructTrgSubmap(submaps_reg, submap_i.overlaps_idx_);
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

    // Plot Pose graph
//    visualizer->updateVisualizer(submaps_reg);
//    visualizer->plotPoseGraphG2O(*graph_obj);
//    while(!viewer.wasStopped ()){
//        viewer.spinOnce ();
//    }
//    viewer.resetStoppedFlag();

    // Benchmark Registered
    PointsT reg_map = pclToMatrixSubmap(submaps_reg);
    PointsT reg_track = trackToMatrixSubmap(submaps_reg);
    benchmark.add_benchmark(reg_map, reg_track, "registered");

    // Create initial graph estimate
//    graph_obj->createInitialEstimate();

    // Save graph to output g2o file (optimization can be run with G2O)
    string outFilename = "graph.g2o";
    graph_obj->saveG2OFile(outFilename);

    // Optimize graph
    google::InitGoogleLogging(argv[0]);
    ceres::optimizer::MapOfPoses poses = ceres::optimizer::ceresSolver(outFilename);

    // Visualize and update submaps with Ceres output
    visualizer->plotPoseGraphCeres(poses, submaps_reg);
    while(!viewer.wasStopped ()){
        viewer.spinOnce ();
    }

    // Benchmark Optimized
    PointsT opt_map = pclToMatrixSubmap(submaps_reg);
    PointsT opt_track = trackToMatrixSubmap(submaps_reg);
    benchmark.add_benchmark(opt_map, opt_track, "optimized");
    benchmark.print_summary();

    delete(gicp_reg);
    delete(graph_obj);
//    delete(visualizer);

    return 0;
}
