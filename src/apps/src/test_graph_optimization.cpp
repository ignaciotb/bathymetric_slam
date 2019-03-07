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

#include "registration/gicp_reg.hpp"
#include "registration/utils_visualization.hpp"

#include "submaps_tools/submaps.hpp"

#include "graph_optimization/utils_g2o.hpp"
#include "graph_optimization/graph_construction.hpp"
#include "graph_optimization/ceres_optimizer.hpp"
#include "graph_optimization/read_g2o.h"

#include "data_tools/benchmark.h"

using namespace Eigen;
using namespace std;
using namespace g2o;

int main(int, char** argv){
    string submaps_dir = argv[1];
    string outFilename = argv[2];

    // Read submaps pointclouds from folder
    SubmapsVec submaps_gt, submaps_init, submaps_reg;
    submaps_gt = readSubmapsInDir(submaps_dir);

    // Visualization
    PCLVisualizer viewer ("Submaps viewer");
    SubmapsVisualizer* visualizer = new SubmapsVisualizer(viewer);
    visualizer->setVisualizer(submaps_gt, 2);
    while(!viewer.wasStopped ()){
        viewer.spinOnce ();
    }
    viewer.resetStoppedFlag();

    // Benchmark GT
    benchmark::track_error_benchmark benchmark("real_data");
    PointsT gt_map = pclToMatrixSubmap(submaps_gt);
    PointsT gt_track = trackofSubmap(submaps_gt);
    benchmark.add_ground_truth(gt_map, gt_track);

    // Add additive Gaussian noise to vehicle's pose among submaps
    GaussianGen transSampler, rotSampler;
    Matrix<double, 6,6> information = generateGaussianNoise(transSampler, rotSampler);
    Isometry3f poseDR(Isometry3f(Translation3f(Vector3f(0,0,0))) *
                      Isometry3f(Quaternionf(1,0,0,0)));
    for(SubmapObj& submap_i: submaps_gt){
        addNoiseToSubmap(transSampler, rotSampler, information, submap_i, poseDR);
    }

    // Benchmar Initial
    PointsT init_map = pclToMatrixSubmap(submaps_gt);
    PointsT init_track = trackToMatrixSubmap(submaps_gt);
    benchmark.add_benchmark(init_map, init_track, "noisy");

    // GICP reg for submaps
    SubmapRegistration* gicp_reg = new SubmapRegistration();

    // Graph constructor
    GraphConstructor* graph_obj = new GraphConstructor();

    SubmapObj submap_trg;
    SubmapsVec submaps_prev;
    for(SubmapObj& submap_i: submaps_gt){
        // Skip first swath
        if(submap_i.swath_id_ > 0){
            // Find overlapping submaps only from previous swaths
            submaps_prev.push_back(submap_i);
            for(SubmapObj& submap_k: submaps_reg){
                if(submap_k.swath_id_ == submap_i.swath_id_-1){
                    submaps_prev.push_back(submap_k);
                }
            }
            submap_i.findOverlaps(submaps_prev);
        }
        // Add submap to registered set
        submaps_init.push_back(submap_i);
        submaps_reg.push_back(submap_i);
        // Update visualizer
        visualizer->updateVisualizer(submaps_reg);
        viewer.spinOnce ();

        // Create graph vertex i
        graph_obj->createNewVertex(submap_i);

        // Create DR edge i and store (skip submap 0)
        if(submap_i.submap_id_ != 0 ){
            graph_obj->createDREdge(submap_i);
        }

        // Potential loop closure detected
        if(!submap_i.overlaps_idx_.empty()){
            std::cout << "Overlaps of submap " << submap_i.submap_id_ << std::endl;
            for(unsigned int j=0; j<submap_i.overlaps_idx_.size(); j++){
                std::cout << submap_i.overlaps_idx_.at(j) << std::endl;
            }
            // Construct target submap
            submap_trg = gicp_reg->constructTrgSubmap(submaps_reg, submap_i.overlaps_idx_);

            // Registration
            if(gicp_reg->gicpSubmapRegistration(submap_trg, submap_i)){
                // If GICP converged, generate loop closure edges
                graph_obj->findLoopClosures(submap_i, submaps_reg, 0.0001);
            }
        }
        submaps_reg.pop_back();
        submaps_reg.push_back(submap_i);

        // Update visualizer
        visualizer->updateVisualizer(submaps_reg);
        viewer.spinOnce ();

        // Cleaning
        submap_trg.submap_pcl_.clear();
        submaps_prev.clear();
    }
    // Create initial graph estimate
    graph_obj->createInitialEstimate();

    // Plot Pose graph
    visualizer->plotPoseGraphG2O(*graph_obj);
    while(!viewer.wasStopped ()){
        viewer.spinOnce ();
    }
    viewer.resetStoppedFlag();

    // Save graph to output g2o file
    graph_obj->saveG2OFile(outFilename);

    // Optimize graph
    google::InitGoogleLogging(argv[0]);
    ceres::optimizer::MapOfPoses poses = ceres::optimizer::ceresSolver(outFilename);

    // Visualize
    visualizer->plotPoseGraphCeres(poses, submaps_reg);
    while(!viewer.wasStopped ()){
        viewer.spinOnce ();
    }

    PointsT opt_map = pclToMatrixSubmap(submaps_reg);
    benchmark.add_benchmark(opt_map, gt_track, "optimized");

    delete(gicp_reg);
    delete(graph_obj);
    delete(visualizer);

    return 0;
}
