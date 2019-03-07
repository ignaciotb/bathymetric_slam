#ifndef GRAPH_CONSTRUCTION_HPP
#define GRAPH_CONSTRUCTION_HPP

#include "g2o/types/slam3d/vertex_se3.h"
#include "g2o/types/slam3d/edge_se3.h"
#include "g2o/core/factory.h"

#include <eigen3/Eigen/Core>

#include "submaps_tools/submaps.hpp"

#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>

using namespace std;
using namespace Eigen;
using namespace g2o;

class GraphConstructor{

private:

public:

    vector<VertexSE3*> vertices_;
    vector<EdgeSE3*> drEdges_;
    vector<EdgeSE3*> edges_;

    GraphConstructor();

    ~GraphConstructor();

    void createNewVertex(SubmapObj& submap);

    void saveG2OFile(std::string outFilename);

    void findLoopClosures(SubmapObj &submap_i,
                          const SubmapsVec& submaps_set, double info_thres);

    void createLCEdge(const SubmapObj& submap_from, const SubmapObj& submap_to);

    void createInitialEstimate();

    void createDREdge(const SubmapObj& submap);
};

#endif // GRAPH_CONSTRUCTION_HPP
