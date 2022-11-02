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

#include "bathy_slam/bathy_slam.hpp"

#include <pcl/filters/voxel_grid.h>

#define INTERACTIVE 0
#define VISUAL 1

using namespace Eigen;
using namespace std;
using namespace g2o;

bool next_step = false;
int current_step = 0;


void add_benchmark(SubmapsVec& submaps, benchmark::track_error_benchmark& benchmark, const string& name, bool is_groundtruth=false) {
    PointsT map = pclToMatrixSubmap(submaps);
    PointsT track = trackToMatrixSubmap(submaps);
    if (is_groundtruth) {
        benchmark.add_ground_truth(map, track);
    } else {
        benchmark.add_benchmark(map, track, name);
    }
}

void benchmark_gt(SubmapsVec& submaps_gt, benchmark::track_error_benchmark& benchmark) {
    // Benchmark GT
    add_benchmark(submaps_gt, benchmark, "original", true);
    ceres::optimizer::saveOriginalTrajectory(submaps_gt); // Save original trajectory to txt
    std::cout << "Visualizing original survey, press space to continue" << std::endl;
}

SubmapsVec build_bathymetric_graph(GraphConstructor& graph_obj, SubmapsVec& submaps_gt,
GaussianGen& transSampler, GaussianGen& rotSampler, YAML::Node config) {

    // GICP reg for submaps
    SubmapRegistration gicp_reg(config);

    // Create SLAM solver and run offline
    std::cout << "Building bathymetric graph with GICP submap registration" << std::endl;
    BathySlam slam_solver(graph_obj, gicp_reg);
    SubmapsVec submaps_reg = slam_solver.runOffline(submaps_gt, transSampler, rotSampler, config);
    std::cout << "Done building graph, press space to continue" << std::endl;
    return submaps_reg;
}

// Create initial graph estimates, optionally add gaussian noise if add_gaussian_noise = true
void create_initial_graph_estimate(GraphConstructor& graph_obj, SubmapsVec& submaps_reg, GaussianGen& transSampler, GaussianGen& rotSampler, bool add_gaussian_noise) {
    std::cout << "Add gaussian noise = " << add_gaussian_noise << std::endl;
    if (add_gaussian_noise) {
        // Add noise to edges on the graph
        graph_obj.addNoiseToGraph(transSampler, rotSampler);
        std::cout << "Gaussian noise added to graph" << std::endl;
    }
    // Create initial DR chain and visualize
    graph_obj.createInitialEstimate(submaps_reg);
    std::cout << "Initial graph estimate constructed, press space to continue" << std::endl;

}

void optimize_graph(GraphConstructor& graph_obj, SubmapsVec& submaps_reg, std::string outFilename, char* argv0, boost::filesystem::path output_path) {
    // Save graph to output g2o file (optimization can be run with G2O)
    graph_obj.saveG2OFile(outFilename);

    // Optimize graph and save to cereal
    google::InitGoogleLogging(argv0);
    ceres::optimizer::MapOfPoses poses = ceres::optimizer::ceresSolver(outFilename, graph_obj.drEdges_.size());
    ceres::optimizer::updateSubmapsCeres(poses, submaps_reg);
    std::cout << "Output cereal: " << boost::filesystem::basename(output_path) << std::endl;
    std::ofstream os(boost::filesystem::basename(output_path) + ".cereal", std::ofstream::binary);
    {
        cereal::BinaryOutputArchive oarchive(os);
        oarchive(submaps_reg);
        os.close();
    }
    std::cout << "Graph optimized, press space to continue" << std::endl;
}

void print_benchmark_results(SubmapsVec& submaps_reg, benchmark::track_error_benchmark& benchmark) {
    benchmark.print_summary();

    std::string command_str = "python ../scripts/plot_results.py --initial_poses poses_original.txt --corrupted_poses poses_corrupted.txt --optimized_poses poses_optimized.txt";
    const char *command = command_str.c_str();
    system(command);
}

void keyboardEventOccurred(const pcl::visualization::KeyboardEvent& event, void* nothing) {
    if (event.getKeySym() == "space" && event.keyDown()) {
        next_step = true;
        current_step++;
    }
}

int main(int argc, char** argv){
    // Inputs
    std::string folder_str, path_str, output_str, simulation, config_path;
    cxxopts::Options options("MyProgram", "One line description of MyProgram");
    options.add_options()
        ("help", "Print help")
        ("simulation", "Simulation data from Gazebo", cxxopts::value(simulation))
        ("bathy_survey", "Input MBES pings in cereal file if simulation = no. If in simulation"
                          "input path to map_small folder", cxxopts::value(path_str))
        ("config", "YAML config file", cxxopts::value(config_path));

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

    YAML::Node config = YAML::LoadFile(config_path);
    std::cout << "Config file: " << config_path << std::endl;
    DRNoise dr_noise = loadDRNoiseFromFile(config);

    // Parse submaps from cereal file
    boost::filesystem::path submaps_path(path_str);
    std::cout << "Input data " << submaps_path << std::endl;

    SubmapsVec submaps_gt, submaps_reg;
    if(simulation == "yes"){
        submaps_gt = readSubmapsInDir(submaps_path.string(), dr_noise);
    }
    else{
        std_data::mbes_ping::PingsT std_pings = std_data::read_data<std_data::mbes_ping::PingsT>(submaps_path);
        std::cout << "Number of pings in survey " << std_pings.size() << std::endl;
        {
            SubmapsVec traj_pings = parsePingsAUVlib(std_pings, dr_noise);
            int submap_size = config["submap_size"].as<int>();
            submaps_gt = createSubmaps(traj_pings, submap_size, dr_noise);

            // Filtering of submaps
            PointCloudT::Ptr cloud_ptr (new PointCloudT);
            pcl::VoxelGrid<PointT> voxel_grid_filter;
            voxel_grid_filter.setInputCloud (cloud_ptr);
            voxel_grid_filter.setLeafSize(config["downsampling_leaf_x"].as<double>(),
                                          config["downsampling_leaf_y"].as<double>(),
                                          config["downsampling_leaf_z"].as<double>());
            for(SubmapObj& submap_i: submaps_gt){
                *cloud_ptr = submap_i.submap_pcl_;
                voxel_grid_filter.setInputCloud(cloud_ptr);
                voxel_grid_filter.filter(*cloud_ptr);
                submap_i.submap_pcl_ = *cloud_ptr;
            }
        }
    }
    std::cout << "Number of submaps " << submaps_gt.size() << std::endl;

    // Graph constructor
    // Read training covs from folder
    covs covs_lc;
    boost::filesystem::path folder(folder_str);
    if(boost::filesystem::is_directory(folder)) {
        covs_lc = readCovsFromFiles(folder);
    }
    GraphConstructor graph_obj(covs_lc);

    // Noise generators
    GaussianGen transSampler, rotSampler;
    Matrix<double, 6,6> information = generateGaussianNoise(transSampler, rotSampler);

    // flag for adding gaussian noise to submaps and graph
    bool add_gaussian_noise = config["add_gaussian_noise"].as<bool>();
    
    benchmark::track_error_benchmark benchmark("real_data", config["benchmark_nbr_rows"].as<int>(), config["benchmark_nbr_cols"].as<int>());
    std::cout << "Benchmark nbr rows and cols: " << benchmark.benchmark_nbr_rows << ", " << benchmark.benchmark_nbr_cols << std::endl;

#if VISUAL != 1
    submaps_reg = build_bathymetric_graph(graph_obj, submaps_gt, transSampler, rotSampler, add_gaussian_noise);
    create_initial_graph_estimate(graph_obj, submaps_reg, transSampler, rotSampler, add_gaussian_noise);
    optimize_graph(graph_obj, submaps_reg, outFilename, argv[0], output_path);
    benchmark_optimized(submaps_reg, benchmark);
#endif

    // Visualization
#if VISUAL == 1
    PCLVisualizer viewer ("Submaps viewer");
    viewer.registerKeyboardCallback(&keyboardEventOccurred, (void*) NULL);
    viewer.loadCameraParameters("Antarctica7");
    SubmapsVisualizer* visualizer = new SubmapsVisualizer(viewer);
    visualizer->setVisualizer(submaps_gt, 1);

    while (!viewer.wasStopped()) {
        viewer.spinOnce();
        if (next_step) {
            next_step = false;
            switch (current_step)
            {
            case 1:
                submaps_reg = build_bathymetric_graph(graph_obj, submaps_gt, transSampler, rotSampler, config);
                visualizer->updateVisualizer(submaps_reg);
                // Benchmark GT
                benchmark_gt(submaps_gt, benchmark);
                break;
            case 2:
                create_initial_graph_estimate(graph_obj, submaps_reg, transSampler, rotSampler, add_gaussian_noise);
                visualizer->plotPoseGraphG2O(graph_obj, submaps_reg);
                // Benchmark corrupted (or not corrupted if add_gaussian_noise = false)
                add_benchmark(submaps_reg, benchmark, "corrupted");
                break;
            case 3:
                optimize_graph(graph_obj, submaps_reg, outFilename, argv[0], output_path);
                // Visualize Ceres output
                visualizer->plotPoseGraphCeres(submaps_reg);
                // Benchmark Optimized
                add_benchmark(submaps_reg, benchmark, "optimized");
                break;
            default:
                break;
            }
        }
    }
    delete(visualizer);
    for (int i = 0; i < 10; i ++) {
        add_benchmark(submaps_gt, benchmark, "gt_" + std::to_string(i));
        add_benchmark(submaps_reg, benchmark, "reg_" + std::to_string(i));
    }
    print_benchmark_results(submaps_reg, benchmark);
#endif

    return 0;
}
