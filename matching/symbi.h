#ifndef MATCHING_SYMBI
#define MATCHING_SYMBI

#include <queue>
#include <unordered_map>
#include <vector>

#include "graph/graph.h"
#include "matching/matching.h"

class SymBi : public matching
{
private:
    struct TreeNode {
        std::vector<uint> forwards_;
        std::vector<uint> forward_labels_;
        std::vector<uint> backwards_;
        std::vector<uint> backward_labels_;
        std::vector<uint> neighbors_;
    };
    struct ExtendableVertex {
        uint E;
        uint matched_nbrs;
        uint u_min;

        ExtendableVertex()
        : E(NOT_EXIST), matched_nbrs(0u), u_min(NOT_EXIST) {}
        ExtendableVertex(uint E_arg, uint matched_nbrs_arg, uint u_min_arg)
        : E(E_arg), matched_nbrs(matched_nbrs_arg), u_min(u_min_arg) {}
    };

    std::vector<std::vector<uint>> eidx_;
    std::vector<TreeNode> treeNode_;
    uint q_root_;
    std::vector<uint> serialized_tree_;
    std::vector<std::vector<uint>> pre_defined_order_;
    std::vector<std::vector<uint>> pre_defined_backward_nbr_;

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

    SymBi(Graph& query_graph, Graph& data_graph, uint max_num_results,
            bool print_prep, bool print_enum, bool homo, 
            std::vector<std::vector<uint>> orders);
    ~SymBi() override {};

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
    
    void InsertionTopDown(uint u, uint u_c, uint v, uint v_c);
    void InsertionBottomUp(uint u, uint u_p, uint v, uint v_p);
    void DeletionTopDown(uint u, uint u_c, uint v, uint v_c);
    void DeletionBottomUp(uint u, uint u_p, uint v, uint v_p);

    void FindMatches(uint depth, std::vector<uint>& m, 
            std::vector<ExtendableVertex>& extendable, size_t &num_results);
    void FindMatches(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results);
};
#endif //MATCHING_SYMBI
