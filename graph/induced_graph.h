#ifndef GRAPH_INDUCEDGRAPH
#define GRAPH_INDUCEDGRAPH

#include <vector>

#include "utils/types.h"
#include "graph/graph.h"

class InducedGraph
{
public:
    const Graph* graph;
    std::vector<uint> v_map_;
    std::vector<std::pair<uint, uint>> e_lists_;

public:
    InducedGraph();
    InducedGraph(const Graph& g, uint v1, uint v2);
    InducedGraph(InducedGraph& g1, InducedGraph& g2, bool is_union);

    uint NumVertices() const { return v_map_.size(); }
    uint NumEdges() const { return e_lists_.size(); }

    uint GetDegree(uint v);
};

#endif //GRAPH_INDUCEDGRAPH
