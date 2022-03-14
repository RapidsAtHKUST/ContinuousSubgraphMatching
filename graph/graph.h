#ifndef GRAPH_GRAPH
#define GRAPH_GRAPH

#include <queue>
#include <tuple>
#include <vector>
#include "utils/types.h"
#include "utils/utils.h"

class Graph
{
protected:
    uint edge_count_;
    uint vlabel_count_;
    uint elabel_count_;
    std::vector<std::vector<uint>> neighbors_;
    std::vector<std::vector<uint>> elabels_;

public:
    std::queue<InsertUnit> updates_;
    std::vector<uint> vlabels_;

public:
    Graph();

    virtual uint NumVertices() const { return vlabels_.size(); }
    virtual uint NumEdges() const { return edge_count_; }
    uint NumVLabels() const { return vlabel_count_; }
    uint NumELabels() const { return elabel_count_; }
    uint GetDiameter() const;

    void AddVertex(uint id, uint label);
    void RemoveVertex(uint id);
    void AddEdge(uint v1, uint v2, uint label);
    void RemoveEdge(uint v1, uint v2);

    uint GetVertexLabel(uint u) const;
    const std::vector<uint>& GetNeighbors(uint v) const;
    const std::vector<uint>& GetNeighborLabels(uint v) const;
    uint GetDegree(uint v) const;
    std::tuple<uint, uint, uint> GetEdgeLabel(uint v1, uint v2) const;

    void LoadFromFile(const std::string &path);
    void LoadUpdateStream(const std::string &path);
    void PrintMetaData() const;
};

#endif //GRAPH_GRAPH
