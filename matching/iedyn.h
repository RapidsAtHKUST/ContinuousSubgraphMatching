#ifndef MATCHING_IEDYN
#define MATCHING_IEDYN

#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "graph/graph.h"
#include "matching/matching.h"

class IEDyn : public matching
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

    std::vector<uint> order_vs_;
    std::vector<uint> backward_vs_;

    std::vector<std::unordered_map<uint, std::vector<uint>>> DCS_;

    std::vector<std::unordered_map<uint, std::vector<uint>>> DCS_CD_;
    std::vector<std::unordered_map<uint, std::vector<uint>>> DCS_CD_delta_;
    std::unordered_set<uint> root_d2_indices_;
    std::unordered_set<uint> root_d2_indices_delta_;

    std::vector<bool> is_delta_node_;
    uint match_node_;

    std::vector<std::unordered_map<uint, bool>> d2;

    std::vector<std::unordered_map<uint, uint>> n2;
    std::vector<std::unordered_map<uint, uint>> nc2;
    
    std::queue<std::pair<uint, uint>> Q2;

public:

    IEDyn(Graph& query_graph, Graph& data_graph, uint max_num_results,
            bool print_prep, bool print_enum, bool homo);
    ~IEDyn() override {};

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
    
    void InsertionTopDown(uint u, uint u_c, uint v, uint v_c);
    void InsertionBottomUp(uint u, uint u_p, uint v, uint v_p);
    void DeletionTopDown(uint u, uint u_c, uint v, uint v_c);
    void DeletionBottomUp(uint u, uint u_p, uint v, uint v_p);

    void CombineDeltaDCS();

    void FindMatchesInitial(uint depth, std::vector<uint>& m, size_t &num_results);
    void FindMatchesIncremental(uint depth, std::vector<uint>& m, size_t &num_results);
};
#endif //MATCHING_IEDYN
