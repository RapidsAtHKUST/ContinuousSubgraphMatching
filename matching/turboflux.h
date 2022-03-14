#ifndef MATCHING_TURBOFLUX
#define MATCHING_TURBOFLUX

#include <queue>
#include <unordered_map>
#include <vector>

#include "graph/graph.h"
#include "matching/matching.h"

class TurboFlux : public matching
{
private:
    struct TreeNode {
        std::vector<uint> forwards_;
        std::vector<uint> forward_labels_;
        std::vector<uint> backwards_;
        std::vector<uint> backward_labels_;
        std::vector<uint> neighbors_;
    };

    std::vector<std::vector<uint>> eidx_;
    std::vector<TreeNode> treeNode_;
    uint q_root_;
    std::vector<uint> serialized_tree_;

    std::vector<std::vector<uint>> order_vs_;
    std::vector<std::vector<uint>> backward_vs_;

    std::vector<std::vector<std::vector<uint>>> join_check_vs_;
    std::vector<std::vector<std::vector<uint>>> join_check_labels_;

    std::vector<std::unordered_map<uint, std::vector<uint>>> DCS_;

    std::vector<std::unordered_map<uint, bool>> d1;
    std::vector<std::unordered_map<uint, bool>> d2;

    std::vector<std::unordered_map<uint, uint>> n1;
    std::vector<std::unordered_map<uint, uint>> np1;

    std::vector<std::unordered_map<uint, uint>> n2;
    std::vector<std::unordered_map<uint, uint>> nc2;
    
    std::queue<std::pair<uint, uint>> Q1;
    std::queue<std::pair<uint, uint>> Q2;

public:

    TurboFlux(Graph& query_graph, Graph& data_graph, uint max_num_results,
            bool print_prep, bool print_enum, bool homo);
    ~TurboFlux() override {};

    void Preprocessing() override;
    void InitialMatching() override;
    
    void AddEdge(uint v1, uint v2, uint label) override;
    void RemoveEdge(uint v1, uint v2) override;
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;
    
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;

private:
    void BuildDAG();
    void BuildDCS();
    void GenerateMatchingOrder();
    void CountDownwards(uint u, uint v, std::vector<uint>& num_explicit_pathes,
        std::vector<std::unordered_map<uint, uint>>& num_dp);
    
    void InsertionTopDown(uint u, uint u_c, uint v, uint v_c);
    void InsertionBottomUp(uint u, uint u_p, uint v, uint v_p);
    void DeletionTopDown(uint u, uint u_c, uint v, uint v_c);
    void DeletionBottomUp(uint u, uint u_p, uint v, uint v_p);

    void FindMatches(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results);
};
#endif //MATCHING_TURBOFLUX
