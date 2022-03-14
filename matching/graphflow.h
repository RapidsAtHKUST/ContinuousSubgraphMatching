#ifndef MATCHING_GRAPHFLOW
#define MATCHING_GRAPHFLOW

#include <vector>

#include "utils/types.h"
#include "graph/graph.h"
#include "matching/matching.h"

class Graphflow : public matching
{
private:
    // a list of matching orders starting from each query edge
    // the first matching order also applies to the initial matching
    std::vector<std::vector<uint>> order_vs_;
    std::vector<std::vector<uint>> order_csrs_;
    std::vector<std::vector<uint>> order_offs_;

public:
    Graphflow(Graph& query_graph, Graph& data_graph, uint max_num_results,
            bool print_prep, bool print_enum, bool homo);
    ~Graphflow() override {};

    void Preprocessing() override;
    void InitialMatching() override;

    void AddEdge(uint v1, uint v2, uint label) override;
    void RemoveEdge(uint v1, uint v2) override;
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;
    
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;

private:
    void GenerateMatchingOrder();
    void FindMatches(uint order_index, uint depth,
            std::vector<uint> m, size_t &num_results);
};

#endif //MATCHING_GRAPHFLOW
