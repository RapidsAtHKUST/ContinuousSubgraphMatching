#ifndef MATCHING_SJ_TREE
#define MATCHING_SJ_TREE

#include <unordered_map>
#include <vector>

#include "utils/types.h"
#include "graph/graph.h"
#include "graph/induced_graph.h"
#include "matching/matching.h"

class SJTree : public matching
{
private:
    struct TreeNode {
        uint parent_ = UINT_MAX;
        // do not need to store children, because only
        // bottom-up traverse is performed
        uint sibling_ = UINT_MAX;

        InducedGraph intersection_; // intersection of child nodes
        InducedGraph graph_;

        std::unordered_map<
            std::pair<uint, uint>, 
            std::vector<uint>,
            pair_hash
        > matches_;
    };
    struct EdgeWithLabel {
        uint v1; uint v2; uint e_label;

        EdgeWithLabel() : v1(), v2(), e_label() {}

        EdgeWithLabel(uint v1_arg, uint v2_arg, uint e_label_arg)
        : v1(v1_arg), v2(v2_arg), e_label(e_label_arg)
        {}
    };
    
    std::vector<EdgeWithLabel> order_es_;
    std::vector<TreeNode> treeNodes_;

public:
    SJTree(Graph& query_graph, Graph& data_graph, uint max_num_results,
            bool print_prep, bool print_enum, bool homo);
    ~SJTree() override {};

    void Preprocessing() override;
    void InitialMatching() override;

    void AddEdge(uint v1, uint v2, uint label) override;
    void RemoveEdge(uint v1, uint v2) override;
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;

    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;

private:
    void GenerateMatchingOrder();
    void BuildSJTree();
    void AddSingleMatch(const std::vector<uint>& match, 
            TreeNode &node, size_t& num_results);
};

#endif //MATCHING_SJ_TREE
