#include "bathy_slam/bathy_slam.hpp"


BathySlam::BathySlam(GraphConstructor &graph_obj, SubmapRegistration &gicp_reg):
    graph_obj_(&graph_obj), gicp_reg_(&gicp_reg){

}

BathySlam::~BathySlam(){

}

SubmapsVec BathySlam::runOffline(SubmapsVec& submaps_gt, GaussianGen& transSampler, GaussianGen& rotSampler, YAML::Node config){
    DRNoise dr_noise = loadDRNoiseFromFile(config);
    SubmapObj submap_trg(dr_noise);
    SubmapsVec submaps_prev, submaps_reg;
    ofstream fileOutputStream;
    fileOutputStream.open("loop_closures.txt", std::ofstream::out);

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
	// Submaps in map_frame?
	bool submaps_in_map_tf = true;
        submap_i.findOverlaps(submaps_in_map_tf, submaps_prev, config["overlap_coverage"].as<double>());
        submaps_prev.clear();

    #if INTERACTIVE == 1
        // Update visualizer
        submaps_reg.push_back(submap_i); // Add submap_i to registered set (just for visualization here)
        visualizer->updateVisualizer(submaps_reg);
        while(!viewer.wasStopped ()){
            viewer.spinOnce ();
        }
        viewer.resetStoppedFlag();
        submaps_reg.pop_back();
    #endif
        // Create graph vertex i
        graph_obj_->createNewVertex(submap_i);

        // Create DR edge i and store (skip submap 0)
        if(submap_i.submap_id_ != 0 ){
            std::cout << "DR edge from " << submap_i.submap_id_ -1 << " to " << submap_i.submap_id_<< std::endl;
            graph_obj_->createDREdge(submap_i);
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
            submap_trg = gicp_reg_->constructTrgSubmap(submaps_reg, submap_i.overlaps_idx_, dr_noise);
            if (config["add_gaussian_noise"].as<bool>()) {
                addNoiseToSubmap(transSampler, rotSampler, submap_i); // Add disturbance to source submap
            }

            if(gicp_reg_->gicpSubmapRegistration(submap_trg, submap_i)){
                submap_final = submap_i;
            }
            submap_trg.submap_pcl_.clear();

            // Create loop closures
            graph_obj_->edge_covs_type_ = config["lc_edge_covs_type"].as<int>();
            graph_obj_->findLoopClosures(submap_final, submaps_reg, info_thres);
        }
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


    return submaps_reg;
}
