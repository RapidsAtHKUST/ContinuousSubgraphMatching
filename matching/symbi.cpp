#include <algorithm>
#include <iostream>
#include <queue>
#include <unordered_set>
#include <vector>

#include "utils/types.h"
#include "utils/globals.h"
#include "graph/graph.h"
#include "matching/symbi.h"

SymBi::SymBi(Graph& query_graph, Graph& data_graph, 
        uint max_num_results,
        bool print_prep, 
        bool print_enum, 
        bool homo,
        std::vector<std::vector<uint>> orders)
: matching(
    query_graph, data_graph, max_num_results, 
    print_prep, print_enum, homo)
, eidx_(query_.NumVertices())
, treeNode_(query_.NumVertices())
, q_root_(0u)
, serialized_tree_(query_.NumVertices())
, pre_defined_order_()
, pre_defined_backward_nbr_()

, DCS_(query_.NumEdges() * 2)
, d1(query_.NumVertices())
, d2(query_.NumVertices())
, n1(query_.NumEdges() * 2)
, np1(query_.NumVertices())
, n2(query_.NumEdges() * 2)
, nc2(query_.NumVertices())
, Q1{}
, Q2{}
{
    uint edge_pos = 0;
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        eidx_[i].resize(query_.NumVertices());
        auto& q_nbrs = query_.GetNeighbors(i);

        for (uint j = 0; j < q_nbrs.size(); j++)
        {
            const uint u_other = q_nbrs[j];
            eidx_[i][u_other] = edge_pos++;
        }
    }
    if (orders.size() == query_.NumEdges())
    {
        pre_defined_order_.resize(query_.NumEdges() * 2);

        for (uint i = 0; i < orders.size(); i++)
        {
            pre_defined_order_[eidx_
                [std::min(orders[i][0], orders[i][1])]
                [std::max(orders[i][0], orders[i][1])]
            ] = orders[i];
        }
    }
    else if (orders.size() == query_.NumEdges() + 1)
    {
        std::vector<std::vector<bool>> is_tree_edge (query_.NumVertices());
        for (uint i = 0; i < query_.NumVertices(); i++)
        {
            is_tree_edge[i].resize(query_.NumVertices(), false);
        }
        for (uint i = 0; i < orders[0].size(); i += 2)
        {
            is_tree_edge[orders[0][i]][orders[0][i + 1]] = true;
            is_tree_edge[orders[0][i + 1]][orders[0][i]] = true;
        }

        pre_defined_order_.resize(query_.NumEdges() * 2);
        pre_defined_backward_nbr_.resize(query_.NumEdges() * 2);

        for (uint i = 1; i < orders.size(); i++)
        {
            pre_defined_order_[eidx_
                [std::min(orders[i][0], orders[i][1])]
                [std::max(orders[i][0], orders[i][1])]
            ] = orders[i];

            auto& backwards = pre_defined_backward_nbr_[eidx_
                [std::min(orders[i][0], orders[i][1])]
                [std::max(orders[i][0], orders[i][1])]
            ];
            backwards.resize(query_.NumVertices());

            for (uint j = 1; j < orders[i].size(); j++)
            {
                for (uint k = 0; k < j; k++)
                {
                    if (
                        std::get<2>(query_.GetEdgeLabel(orders[i][j], orders[i][k])) != NOT_EXIST &&
                        is_tree_edge[orders[i][j]][orders[i][k]]
                    ) {
                        backwards[orders[i][j]] = orders[i][k];
                    }
                }
            }
        }
    }
}

void SymBi::Preprocessing()
{
    BuildDAG();
    BuildDCS();
}

void SymBi::BuildDAG()
{
    struct EdgeFreq{
        uint v1; uint v2; uint freq;
        EdgeFreq()
        : v1(NOT_EXIST), v2(NOT_EXIST), freq(NOT_EXIST)
        {}
        EdgeFreq(uint v1_arg, uint v2_arg, uint freq_arg)
        : v1(v1_arg), v2(v2_arg), freq(freq_arg)
        {}
    };

    // get the frequency of each query edge
    std::vector<std::vector<size_t>> edge_freqs(query_.NumVertices());
    for (auto& vec: edge_freqs) vec.resize(query_.NumVertices(), 0ul);
    size_t min_freq = ULONG_MAX;
    std::array<uint, 2> min_vs {};

    for (uint u = 0; u < query_.NumVertices(); u++)
    {
        const auto& q_nbrs = query_.GetNeighbors(u);
        const auto& q_nbr_labels = query_.GetNeighborLabels(u);

        for (uint i = 0; i < q_nbrs.size(); i++)
        {
            const auto& u_other = q_nbrs[i];
            const auto& u_other_label = q_nbr_labels[i];
            
            if (u > u_other)
            {
                edge_freqs[u][u_other] = edge_freqs[u_other][u];
                continue;
            }

            for (size_t v = 0; v < data_.NumVertices(); v++)
            if (data_.GetVertexLabel(v) != NOT_EXIST)
            {
                const auto& d_nbrs = data_.GetNeighbors(v);
                const auto& d_nbr_labels = data_.GetNeighborLabels(v);

                for (uint j = 0; j < d_nbrs.size(); j++)
                {
                    const auto& v_other = d_nbrs[j];
                    const auto& v_other_label = d_nbr_labels[j];
                    
                    if (
                        data_.GetVertexLabel(v) == query_.GetVertexLabel(u) &&
                        data_.GetVertexLabel(v_other) == query_.GetVertexLabel(u_other) &&
                        v_other_label == u_other_label
                    ) {
                        edge_freqs[u][u_other]++;
                    }
                }
            }
            if (edge_freqs[u][u_other] < min_freq)
            {
                min_freq = edge_freqs[u][u_other];
                min_vs[0] = u;
                min_vs[1] = u_other;
            }
        }
    }

    // find vertex with low frequency
    std::array<uint, 2> u_freq {};
    for (size_t v = 0; v < data_.NumVertices(); v++)
    if (data_.GetVertexLabel(v) != NOT_EXIST)
    {
        if (data_.GetVertexLabel(v) == query_.GetVertexLabel(min_vs[0]))
            u_freq[0] ++;
        if (data_.GetVertexLabel(v) == query_.GetVertexLabel(min_vs[1]))
            u_freq[1] ++;
    }
    if (u_freq[0] > u_freq[1])
        q_root_ = min_vs[1];
    else if (u_freq[0] < u_freq[1])
        q_root_ = min_vs[0];
    else if (query_.GetDegree(min_vs[0]) > query_.GetDegree(min_vs[1]))
        q_root_ = min_vs[0];
    else
        q_root_ = min_vs[1];

    // build a spaning tree with the BFS order
    std::vector<bool> visited(query_.NumVertices(), false);
    std::vector<bool> exist_on_tree(query_.NumVertices(), false);
    std::queue<uint> bfs_queue;

    bfs_queue.push(q_root_);
    exist_on_tree[q_root_] = true;
    uint order_pos = 0u;
    serialized_tree_[order_pos++] = q_root_;

    while (!bfs_queue.empty())
    {
        const uint u = bfs_queue.front();
        bfs_queue.pop();
        visited[u] = true;

        const auto& q_nbrs = query_.GetNeighbors(u);
        const auto& q_nbr_labels = query_.GetNeighborLabels(u);

        for (uint j = 0; j < q_nbrs.size(); j++)
        {
            const uint u_other = q_nbrs[j];
            const uint u_other_label = q_nbr_labels[j];
            
            if (!exist_on_tree[u_other])
            {
                bfs_queue.push(u_other);
                exist_on_tree[u_other] = true;
                serialized_tree_[order_pos++] = u_other;
            }
            if (!visited[u_other])
            {
                treeNode_[u].forwards_.push_back(u_other);
                treeNode_[u].forward_labels_.push_back(u_other_label);
            }
            else
            {
                treeNode_[u].backwards_.push_back(u_other);
                treeNode_[u].backward_labels_.push_back(u_other_label);
            }
            treeNode_[u].neighbors_.push_back(u_other);
        }
    }  
    if (print_preprocessing_results_)
    {
        std::cout << "DAG: " << std::endl;
        for (uint i = 0; i < query_.NumVertices(); ++i)
        {
            std::cout << serialized_tree_[i] << ": (backwards: ";
            for (auto j: treeNode_[serialized_tree_[i]].backwards_)
                std::cout << j << " ";

            std::cout << ") (forwards: ";
            for (auto j: treeNode_[serialized_tree_[i]].forwards_)
                std::cout << j << " ";

            std::cout << ")" << std::endl;
        }
        std::cout << std::endl;
    }
}

void SymBi::BuildDCS()
{
    // top-down
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const uint u = serialized_tree_[i];
        const uint query_label = query_.GetVertexLabel(u);
        
        for (size_t j = 0; j < data_.NumVertices(); j++)
        if (data_.GetVertexLabel(j) == query_label)
        {
            // add j to hash maps
            const auto& q_nbrs = query_.GetNeighbors(u);
            for (uint k = 0; k < q_nbrs.size(); k++)
            {
                DCS_[eidx_[u][q_nbrs[k]]][j] = {};
            }
            for (uint k = 0; k < treeNode_[u].backwards_.size(); k++)
            {
                const uint u_other = treeNode_[u].backwards_[k];
                const uint u_other_elabel = treeNode_[u].backward_labels_[k];

                const auto& d_nbrs = data_.GetNeighbors(j);
                const auto& d_elabels = data_.GetNeighborLabels(j);

                for (uint m = 0; m < d_nbrs.size(); m++)
                {
                    const uint v_other = d_nbrs[m];
                    const uint v_other_elabel = d_elabels[m];
                    if (
                        query_.GetVertexLabel(u_other) == data_.GetVertexLabel(v_other) &&
                        u_other_elabel == v_other_elabel
                    ) {
                        DCS_[eidx_[u][u_other]][j].push_back(v_other);
                        if (d1[u_other][v_other] > 0)
                            n1[eidx_[u][u_other]][j] += 1;
                    }
                }
                if (n1[eidx_[u][u_other]][j])
                {
                    np1[u][j] += 1;
                }
            }
            if (np1[u][j] == treeNode_[u].backwards_.size())
            {
                d1[u][j] = 1;
            }
        }
    }

    // bottom-up
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const uint u = serialized_tree_[query_.NumVertices() - 1 - i];
        const uint query_label = query_.GetVertexLabel(u);
        
        for (size_t j = 0; j < data_.NumVertices(); j++)
        if (data_.GetVertexLabel(j) == query_label)
        {
            for (uint k = 0; k < treeNode_[u].forwards_.size(); k++)
            {
                const uint u_other = treeNode_[u].forwards_[k];
                const uint u_other_elabel = treeNode_[u].forward_labels_[k];

                const auto& d_nbrs = data_.GetNeighbors(j);
                const auto& d_elabels = data_.GetNeighborLabels(j);

                for (uint m = 0; m < d_nbrs.size(); m++)
                {
                    const uint v_other = d_nbrs[m];
                    const uint v_other_elabel = d_elabels[m];
                    if (
                        query_.GetVertexLabel(u_other) == data_.GetVertexLabel(v_other) &&
                        u_other_elabel == v_other_elabel
                    ) {
                        DCS_[eidx_[u][u_other]][j].push_back(v_other);
                        if (d2[u_other][v_other] > 0)
                            n2[eidx_[u][u_other]][j] += 1;
                    }
                }
                if (n2[eidx_[u][u_other]][j])
                {
                    nc2[u][j] += 1;
                }
            }
            if (nc2[u][j] == treeNode_[u].forwards_.size() && d1[u][j] == 1)
            {
                d2[u][j] = 1;
            }
        }
    }

    // top-down modify n2
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const uint u = serialized_tree_[i];
        for (uint k = 0; k < treeNode_[u].backwards_.size(); k++)
        {
            const uint u_other = treeNode_[u].backwards_[k];
            
            for (auto & [v, d_nbrs]: DCS_[eidx_[u][u_other]])
            {
                for (uint m = 0; m < d_nbrs.size(); m++)
                {
                    const uint v_other = d_nbrs[m];

                    if (d2[u_other][v_other] > 0) {
                        n2[eidx_[u][u_other]][v] += 1;
                    }
                }
            }
        }
    }
}

void SymBi::InitialMatching()
{
    std::vector<uint> m (query_.NumVertices(), UNMATCHED);
    
    for (auto & [v, d_nbrs]: DCS_[eidx_[q_root_][query_.GetNeighbors(q_root_).front()]])
    {
        if (d2[q_root_][v] == false) continue;

        m[q_root_] = v;
        visited_[v] = true;

        std::vector<ExtendableVertex> extendable(query_.NumVertices());
        for (auto& u_other: treeNode_[q_root_].neighbors_)
        {
            if (m[u_other] != UNMATCHED) continue;

            if (n2[eidx_[q_root_][u_other]][v] < extendable[u_other].E)
            {
                extendable[u_other].E = n2[eidx_[q_root_][u_other]][v];
                extendable[u_other].u_min = q_root_;
            }
            extendable[u_other].matched_nbrs ++;
        }
        FindMatches(1, m, extendable, num_initial_results_);

        visited_[v] = false;
        m[q_root_] = UNMATCHED;
    }
}

void SymBi::InsertionTopDown(uint u, uint u_c, uint v, uint v_c)
{
    if (n1[eidx_[u_c][u]][v_c] == 0)
    {
        np1[u_c][v_c] += 1;
        if (np1[u_c][v_c] == treeNode_[u_c].backwards_.size())
        {
            d1[u_c][v_c] = 1;
            Q1.emplace(u_c, v_c);
            if (nc2[u_c][v_c] == treeNode_[u_c].forwards_.size())
            {
                d2[u_c][v_c] = 1;
                Q2.emplace(u_c, v_c);
            }
        }
    }
    n1[eidx_[u_c][u]][v_c] += 1;
}

void SymBi::InsertionBottomUp(uint u, uint u_p, uint v, uint v_p)
{
    if (n2[eidx_[u_p][u]][v_p] == 0)
    {
        nc2[u_p][v_p] += 1;
        if (d1[u_p][v_p] && nc2[u_p][v_p] == treeNode_[u_p].forwards_.size())
        {
            d2[u_p][v_p] = 1;
            Q2.emplace(u_p, v_p);
        }
    }
    n2[eidx_[u_p][u]][v_p] += 1;
}

void SymBi::DeletionTopDown(uint u, uint u_c, uint v, uint v_c)
{
    n1[eidx_[u_c][u]][v_c] -= 1;
    if (n1[eidx_[u_c][u]][v_c] == 0)
    {
        if (d1[u_c][v_c] == 1)
        {
            if (d2[u_c][v_c] == 1)
            {
                Q2.emplace(u_c, v_c);
                d2[u_c][v_c] = 0;
            }
            Q1.emplace(u_c, v_c);
            d1[u_c][v_c] = 0;
        }
        np1[u_c][v_c] -= 1;
    }
}

void SymBi::DeletionBottomUp(uint u, uint u_p, uint v, uint v_p)
{
    n2[eidx_[u_p][u]][v_p] -= 1;
    if (n2[eidx_[u_p][u]][v_p] == 0)
    {
        if (d2[u_p][v_p] == 1)
        {
            Q2.emplace(u_p, v_p);
            d2[u_p][v_p] = 0;
        }
        nc2[u_p][v_p] -= 1;
    }
}

void SymBi::FindMatches(uint depth, std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, size_t &num_results)
{
    if (reach_time_limit) return;

    uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
    uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;

    uint u;
    for (uint i = 0; i < extendable.size(); i++)
    {
        if (m[i] != UNMATCHED) continue;

        // check if a extendable query vertex is isolated or not
        if (extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
        {
            if (extendable[i].E < isolate_minE)
            {
                isolate_minE = extendable[i].E;
                isolate_u = i;
            }
        }
        else
        {
            if (extendable[i].E < non_isolate_minE)
            {
                non_isolate_minE = extendable[i].E;
                non_isolate_u = i;
            }
        }
    }

    // if no non-isolated vertex exists, then choose an isolated one
    if (non_isolate_minE == NOT_EXIST)
        u = isolate_u;
    else
        u = non_isolate_u;
    uint u_min = extendable[u].u_min;
    extendable[u] = {};

    // enumerate each neighbor of m[u_min]
    bool candidate_empty = true;
    for (auto& v: DCS_[eidx_[u_min][u]][m[u_min]])
    {
        // 1. check index
        num_intermediate_results_before_index_check_++;
        if (d2[u][v] == 0) continue;
        num_intermediate_results_after_index_check_++;

        // 2. check if joinable
        bool joinable = true;
        for (auto& u_other: treeNode_[u].neighbors_)
        {
            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
            if (
                it == DCS_[eidx_[u_other][u]][m[u_other]].end() ||
                *it != v
            ) {
                joinable = false;
                break;
            }
        }
        if (!joinable) continue;
        num_intermediate_results_after_joinability_check_++;

        candidate_empty = false;

        // 3. check if visited
        if (!homomorphism_ && visited_[v]) continue;
        num_intermediate_results_after_visit_check_++;

        // 4. add a vertex mapping
        m[u] = v;
        visited_[v] = true;
        std::vector<ExtendableVertex> temp_extendable(extendable);
        for (auto& u_other: treeNode_[u].neighbors_)
        {
            if (m[u_other] != UNMATCHED) continue;

            if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E)
            {
                temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
                temp_extendable[u_other].u_min = u;
            }
            temp_extendable[u_other].matched_nbrs ++;
        }

        if (depth == query_.NumVertices() - 1)
        {
            num_results++;
            if (print_enumeration_results_)
            {
                for (auto j: m)
                    std::cout << j << " ";
                std::cout << std::endl;
            }
        }
        else
        {
            size_t num_results_before_recursion = num_results;
            FindMatches(depth + 1, m, temp_extendable, num_results);
            if (num_results == num_results_before_recursion)
            {
                num_intermediate_results_without_results_++;
            }
        }

        visited_[v] = false;
        m[u] = UNMATCHED;
        if (num_results >= max_num_results_) return;
        if (reach_time_limit) return;
    }
    if (candidate_empty) num_intermediate_results_with_empty_candidate_set_++;
}

void SymBi::FindMatches(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results)
{
    if (reach_time_limit) return;

    uint u = pre_defined_order_[order_index][depth];

    uint u_min = NOT_EXIST;
    if (pre_defined_backward_nbr_.empty())
    {
        uint u_min_size = UINT_MAX;
        for (auto& u_other: treeNode_[u].neighbors_)
        {
            if (m[u_other] == UNMATCHED) continue;
            if (DCS_[eidx_[u_other][u]][m[u_other]].size() < u_min_size)
            {
                u_min_size = DCS_[eidx_[u_other][u]][m[u_other]].size();
                u_min = u_other;
            }
        }
    }
    else
    {
        u_min = pre_defined_backward_nbr_[order_index][u];
    }

    // enumerate each neighbor of m[u_min]
    bool candidate_empty = true;
    for (auto& v: DCS_[eidx_[u_min][u]][m[u_min]])
    {
        // 1. check index
        num_intermediate_results_before_index_check_++;
        if (d2[u][v] == 0) continue;
        num_intermediate_results_after_index_check_++;

        // 2. check if joinable
        bool joinable = true;
        for (auto& u_other: treeNode_[u].neighbors_)
        {
            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
            if (
                it == DCS_[eidx_[u_other][u]][m[u_other]].end() ||
                *it != v
            ) {
                joinable = false;
                break;
            }
        }
        if (!joinable) continue;
        num_intermediate_results_after_joinability_check_++;

        candidate_empty = false;

        // 3. check if visited
        if (!homomorphism_ && visited_[v]) continue;
        num_intermediate_results_after_visit_check_++;
        
        // 4. add a vertex mapping
        m[u] = v;
        visited_[v] = true;

        if (depth == query_.NumVertices() - 1)
        {
            num_results++;
            if (print_enumeration_results_)
            {
                for (auto j: m)
                    std::cout << j << " ";
                std::cout << std::endl;
            }
        }
        else
        {
            size_t num_results_before_recursion = num_results;
            FindMatches(order_index, depth + 1, m, num_results);
            if (num_results == num_results_before_recursion)
            {
                num_intermediate_results_without_results_++;
            }
        }
        
        visited_[v] = false;
        m[u] = UNMATCHED;
        if (num_results >= max_num_results_) return;
        if (reach_time_limit) return;
    }
    if (candidate_empty) num_intermediate_results_with_empty_candidate_set_++;
}

void SymBi::AddEdge(uint v1, uint v2, uint label)
{
    data_.AddEdge(v1, v2, label);

    // enumerate all query edges that matches v1 --> v2
    for (uint u1 = 0; u1 < query_.NumVertices(); u1++)
    if (data_.GetVertexLabel(v1) == query_.GetVertexLabel(u1))
    {
    for (uint u2 = 0; u2 < query_.NumVertices(); u2++)
    if (data_.GetVertexLabel(v2) == query_.GetVertexLabel(u2))
    {
        if (std::get<2>(query_.GetEdgeLabel(u1, u2)) != label) continue;
        
        bool reversed = false;
        if (std::find(treeNode_[u1].backwards_.begin(), treeNode_[u1].backwards_.end(), u2) != treeNode_[u1].backwards_.end())
        {
            std::swap(u1, u2);
            std::swap(v1, v2);
            reversed = true;
        }
        if (std::find(treeNode_[u2].backwards_.begin(), treeNode_[u2].backwards_.end(), u1) != treeNode_[u2].backwards_.end())
        {
            auto it = std::lower_bound(DCS_[eidx_[u1][u2]][v1].begin(), DCS_[eidx_[u1][u2]][v1].end(), v2);
            DCS_[eidx_[u1][u2]][v1].insert(it, v2);
            it = std::lower_bound(DCS_[eidx_[u2][u1]][v2].begin(), DCS_[eidx_[u2][u1]][v2].end(), v1);
            DCS_[eidx_[u2][u1]][v2].insert(it, v1);

            bool old_p_d1 = d1[u1][v1], old_p_d2 = d2[u1][v1], old_c_d2 = d2[u2][v2];

            if (old_p_d1)
                InsertionTopDown(u1, u2, v1, v2);
            if (old_c_d2)
                InsertionBottomUp(u2, u1, v2, v1);
            if (old_p_d2)
                n2[eidx_[u2][u1]][v2] += 1;

            while (!Q1.empty())
            {
                auto [u_queue, v_queue] = Q1.front();
                Q1.pop();
                for (auto& u_c_queue : treeNode_[u_queue].forwards_)
                for (auto& v_c_queue : DCS_[eidx_[u_queue][u_c_queue]][v_queue])
                {
                    InsertionTopDown(u_queue, u_c_queue, v_queue, v_c_queue);
                    if (reach_time_limit) return;
                }
            }
            while (!Q2.empty())
            {
                auto [u_queue, v_queue] = Q2.front();
                Q2.pop();
                for (auto& u_p_queue : treeNode_[u_queue].backwards_)
                for (auto& v_p_queue : DCS_[eidx_[u_queue][u_p_queue]][v_queue])
                {
                    InsertionBottomUp(u_queue, u_p_queue, v_queue, v_p_queue);
                    if (reach_time_limit) return;
                }
                for (auto& u_c_queue : treeNode_[u_queue].forwards_)
                for (auto& v_c_queue : DCS_[eidx_[u_queue][u_c_queue]][v_queue])
                {
                    n2[eidx_[u_c_queue][u_queue]][v_c_queue] += 1;
                    if (reach_time_limit) return;
                }
            }
        }
        if (reversed)
        {
            std::swap(u1, u2);
            std::swap(v1, v2);
        }
    }
    }
    if (max_num_results_ == 0) return;
    
    std::vector<uint> m(query_.NumVertices(), UNMATCHED);

    // enumerate all query edges that matches v1 --> v2
    size_t num_results = 0ul;
    for (uint u1 = 0; u1 < query_.NumVertices(); u1++)
    if (data_.GetVertexLabel(v1) == query_.GetVertexLabel(u1))
    {
    for (uint u2 = 0; u2 < query_.NumVertices(); u2++)
    if (data_.GetVertexLabel(v2) == query_.GetVertexLabel(u2))
    {
        if (std::get<2>(query_.GetEdgeLabel(u1, u2)) != label) continue;
        
        bool reversed = false;
        if (std::find(treeNode_[u1].backwards_.begin(), treeNode_[u1].backwards_.end(), u2) != treeNode_[u1].backwards_.end())
        {
            std::swap(u1, u2);
            std::swap(v1, v2);
            reversed = true;
        }
        if (std::find(treeNode_[u2].backwards_.begin(), treeNode_[u2].backwards_.end(), u1) != treeNode_[u2].backwards_.end()
            && d2[u1][v1] == 1 && d2[u2][v2] == 1)
        {
            m[u1] = v1;
            m[u2] = v2;
            visited_[v1] = true;
            visited_[v2] = true;

            std::vector<ExtendableVertex> extendable(query_.NumVertices());
            for (auto u: {u1, u2})
            {
                for (auto& u_other: treeNode_[u].neighbors_)
                {
                    if (m[u_other] != UNMATCHED) continue;

                    if (n2[eidx_[u][u_other]][m[u]] < extendable[u_other].E)
                    {
                        extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
                        extendable[u_other].u_min = u;
                    }
                    extendable[u_other].matched_nbrs ++;
                }
            }
            if (pre_defined_order_.empty())
                FindMatches(2, m, extendable, num_results);
            else
                FindMatches(eidx_[std::min(u1, u2)][std::max(u1, u2)], 2, m, 
                        num_results);

            visited_[v1] = false;
            visited_[v2] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;
            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) goto END_ENUMERATION;
        }
        if (reversed)
        {
            std::swap(u1, u2);
            std::swap(v1, v2);
        }
    }
    }
    END_ENUMERATION:
    num_positive_results_ += num_results;
}

void SymBi::RemoveEdge(uint v1, uint v2)
{
    std::vector<uint> m(query_.NumVertices(), UNMATCHED);

    size_t num_results = 0ul;
    if (max_num_results_ == 0) goto END_ENUMERATION;

    // enumerate all query edges that matches v1 --> v2
    for (uint u1 = 0; u1 < query_.NumVertices(); u1++)
    if (data_.GetVertexLabel(v1) == query_.GetVertexLabel(u1))
    {
    for (uint u2 = 0; u2 < query_.NumVertices(); u2++)
    if (data_.GetVertexLabel(v2) == query_.GetVertexLabel(u2))
    {
        auto it = std::lower_bound(DCS_[eidx_[u1][u2]][v1].begin(), DCS_[eidx_[u1][u2]][v1].end(), v2);
        if (
            it == DCS_[eidx_[u1][u2]][v1].end() ||
            *it != v2
        ) continue;
        
        bool reversed = false;
        if (std::find(treeNode_[u1].backwards_.begin(), treeNode_[u1].backwards_.end(), u2) != treeNode_[u1].backwards_.end())
        {
            std::swap(u1, u2);
            std::swap(v1, v2);
            reversed = true;
        }
        if (std::find(treeNode_[u2].backwards_.begin(), treeNode_[u2].backwards_.end(), u1) != treeNode_[u2].backwards_.end()
            && d2[u1][v1] == 1 && d2[u2][v2] == 1)
        {
            m[u1] = v1;
            m[u2] = v2;
            visited_[v1] = true;
            visited_[v2] = true;

            std::vector<ExtendableVertex> extendable(query_.NumVertices());
            for (auto u: {u1, u2})
            {
                for (auto& u_other: treeNode_[u].neighbors_)
                {
                    if (m[u_other] != UNMATCHED) continue;

                    if (n2[eidx_[u][u_other]][m[u]] < extendable[u_other].E)
                    {
                        extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
                        extendable[u_other].u_min = u;
                    }
                    extendable[u_other].matched_nbrs ++;
                }
            }
            if (pre_defined_order_.empty())
                FindMatches(2, m, extendable, num_results);
            else
                FindMatches(eidx_[std::min(u1, u2)][std::max(u1, u2)], 2, m, 
                        num_results);

            visited_[v1] = false;
            visited_[v2] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;
            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) return;
        }
        if (reversed)
        {
            std::swap(u1, u2);
            std::swap(v1, v2);
        }
    }
    }

    END_ENUMERATION:
    num_negative_results_ += num_results;

    // enumerate all query edges that matches v1 --> v2
    for (uint u1 = 0; u1 < query_.NumVertices(); u1++)
    if (data_.GetVertexLabel(v1) == query_.GetVertexLabel(u1))
    {
    for (uint u2 = 0; u2 < query_.NumVertices(); u2++)
    if (data_.GetVertexLabel(v2) == query_.GetVertexLabel(u2))
    {
        auto it = std::lower_bound(DCS_[eidx_[u1][u2]][v1].begin(), DCS_[eidx_[u1][u2]][v1].end(), v2);
        if (
            it == DCS_[eidx_[u1][u2]][v1].end() ||
            *it != v2
        ) continue;
        
        bool reversed = false;
        if (std::find(treeNode_[u1].backwards_.begin(), treeNode_[u1].backwards_.end(), u2) != treeNode_[u1].backwards_.end())
        {
            std::swap(u1, u2);
            std::swap(v1, v2);
            reversed = true;
        }
        if (std::find(treeNode_[u2].backwards_.begin(), treeNode_[u2].backwards_.end(), u1) != treeNode_[u2].backwards_.end())
        {
            auto it = std::lower_bound(DCS_[eidx_[u1][u2]][v1].begin(), DCS_[eidx_[u1][u2]][v1].end(), v2);
            DCS_[eidx_[u1][u2]][v1].erase(it);
            it = std::lower_bound(DCS_[eidx_[u2][u1]][v2].begin(), DCS_[eidx_[u2][u1]][v2].end(), v1);
            DCS_[eidx_[u2][u1]][v2].erase(it);

            bool old_p_d1 = d1[u1][v1], old_p_d2 = d2[u1][v1], old_c_d2 = d2[u2][v2];

            if (old_c_d2)
                DeletionBottomUp(u2, u1, v2, v1);
            if (old_p_d2)
                n2[eidx_[u2][u1]][v2] -= 1;
            if (old_p_d1)
                DeletionTopDown(u1, u2, v1, v2);

            while (!Q1.empty())
            {
                auto [u_queue, v_queue] = Q1.front();
                Q1.pop();
                for (auto& u_c_queue : treeNode_[u_queue].forwards_)
                for (auto& v_c_queue : DCS_[eidx_[u_queue][u_c_queue]][v_queue])
                {
                    DeletionTopDown(u_queue, u_c_queue, v_queue, v_c_queue);
                    if (reach_time_limit) return;
                }
            }
            while (!Q2.empty())
            {
                auto [u_queue, v_queue] = Q2.front();
                Q2.pop();
                for (auto& u_p_queue : treeNode_[u_queue].backwards_)
                for (auto& v_p_queue : DCS_[eidx_[u_queue][u_p_queue]][v_queue])
                {
                    DeletionBottomUp(u_queue, u_p_queue, v_queue, v_p_queue);
                    if (reach_time_limit) return;
                }
                for (auto& u_c_queue : treeNode_[u_queue].forwards_)
                for (auto& v_c_queue : DCS_[eidx_[u_queue][u_c_queue]][v_queue])
                {
                    n2[eidx_[u_c_queue][u_queue]][v_c_queue] -= 1;
                    if (reach_time_limit) return;
                }
            }
        }
        if (reversed)
        {
            std::swap(u1, u2);
            std::swap(v1, v2);
        }
    }
    }

    data_.RemoveEdge(v1, v2);
}

void SymBi::AddVertex(uint id, uint label)
{
    data_.AddVertex(id, label);

    visited_.resize(id + 1, false);

    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const uint u = serialized_tree_[i];
        if (data_.GetVertexLabel(id) == query_.GetVertexLabel(u))
        {
            // add j to hash maps
            const auto& q_nbrs = query_.GetNeighbors(u);
            for (uint k = 0; k < q_nbrs.size(); k++)
            {
                DCS_[eidx_[u][q_nbrs[k]]][id] = {};
            }
            if (treeNode_[u].backwards_.empty())
                d1[u][id] = 1;
            if (d1[u][id] == 1 && treeNode_[u].forwards_.empty())
                d2[u][id] = 1;
        }
    }
}

void SymBi::RemoveVertex(uint id)
{
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const uint u = serialized_tree_[i];
        if (data_.GetVertexLabel(id) == query_.GetVertexLabel(u))
        {
            // add j to hash maps
            const auto& q_nbrs = query_.GetNeighbors(u);
            for (uint k = 0; k < q_nbrs.size(); k++)
            {
                DCS_[eidx_[u][q_nbrs[k]]].erase(id);
                if (n1[eidx_[u][q_nbrs[k]]].find(id) != n1[eidx_[u][q_nbrs[k]]].end())
                    n1[eidx_[u][q_nbrs[k]]].erase(id);
                if (n2[eidx_[u][q_nbrs[k]]].find(id) != n2[eidx_[u][q_nbrs[k]]].end())
                    n2[eidx_[u][q_nbrs[k]]].erase(id);
            }
            if (np1[u].find(id) != np1[u].end())
                np1[u].erase(id);
            if (nc2[u].find(id) != nc2[u].end())
                nc2[u].erase(id);
            if (d1[u].find(id) != d1[u].end())
                d1[u].erase(id);
            if (d2[u].find(id) != d2[u].end())
                d2[u].erase(id);
        }
    }

    data_.RemoveVertex(id);
}

void SymBi::GetMemoryCost(size_t &num_edges, size_t &num_vertices)
{
    num_edges = 0ul;
    num_vertices = 0ul;

    std::vector<std::unordered_set<uint>> all_vertices_in_index(query_.NumVertices());
    std::vector<std::unordered_set<uint>> d1_vertices_in_index(query_.NumVertices());
    std::vector<std::unordered_set<uint>> d2_vertices_in_index(query_.NumVertices());

    std::vector<size_t> num_all_edges(query_.NumEdges() * 2, 0);
    std::vector<size_t> num_d1_edges(query_.NumEdges() * 2, 0);
    std::vector<size_t> num_d2_edges(query_.NumEdges() * 2, 0);

    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const auto& children = treeNode_[i].forwards_;

        if (!children.empty())
        {
            for (auto& u_other: children)
            {
                const uint eidx = eidx_[i][u_other];
                for (const auto& [v, d_nbrs]: DCS_[eidx])
                {
                    if (d_nbrs.empty()) continue;

                    all_vertices_in_index[i].insert(v);
                    if (d1[i][v])
                    {
                        d1_vertices_in_index[i].insert(v);
                    }
                    if (d2[i][v])
                    {
                        d2_vertices_in_index[i].insert(v);
                    }

                    for (uint j = 0; j < d_nbrs.size(); j++)
                    {
                        num_all_edges[eidx] += 1;
                        if (d1[u_other][d_nbrs[j]])
                        {
                            num_d1_edges[eidx]++;
                        }
                        if (d2[i][v])
                        {
                            num_d2_edges[eidx]++;
                        }
                    }
                }
            }
        }
        else
        {
            const auto& parents = treeNode_[i].backwards_;
            
            for (auto& u_other: parents)
            {
                const uint eidx = eidx_[i][u_other];
                for (const auto& [v, d_nbrs]: DCS_[eidx])
                {
                    if (d_nbrs.empty()) continue;

                    all_vertices_in_index[i].insert(v);
                    if (d1[i][v])
                        d1_vertices_in_index[i].insert(v);
                    if (d2[i][v])
                        d2_vertices_in_index[i].insert(v);
                }
            }
        }
    }

    size_t total_num_candidate_vs = 0ul;
    size_t total_num_valid_candidate_vs = 0ul;
    size_t total_num_candidate_es = 0ul;
    size_t total_num_valid_candidate_es = 0ul;

    std::cout << "# vertices in index: "; 
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        std::cout << i << ": " << all_vertices_in_index[i].size() << " ";
        num_vertices += all_vertices_in_index[i].size();
    }
    std::cout << "\n# d1 vertex in index: "; 
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        std::cout << i << ": " << d1_vertices_in_index[i].size() << " ";
        total_num_candidate_vs += d1_vertices_in_index[i].size();
    }
    std::cout << "\n# d2 vertex in index: "; 
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        std::cout << i << ": " << d2_vertices_in_index[i].size() << " ";
        total_num_valid_candidate_vs += d2_vertices_in_index[i].size();
    }

    std::cout << "\n# edges in index: ";
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const auto& children = treeNode_[i].forwards_;

        for (auto& u_other: children)
        {
            const uint eidx = eidx_[i][u_other];
            std::cout << i << "-" << u_other << ": " <<
                num_all_edges[eidx] << " ";
            
            num_edges += num_all_edges[eidx];
        }
    }
    std::cout << "\n# d1 edges in index: ";
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const auto& children = treeNode_[i].forwards_;

        for (auto& u_other: children)
        {
            const uint eidx = eidx_[i][u_other];
            std::cout << i << "-" << u_other << ": " <<
                num_d1_edges[eidx] << " ";

            total_num_candidate_es += num_d1_edges[eidx];
        }
    }
    std::cout << "\n# d2 edges in index: ";
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const auto& children = treeNode_[i].forwards_;

        for (auto& u_other: children)
        {
            const uint eidx = eidx_[i][u_other];
            std::cout << i << "-" << u_other << ": " <<
                num_d2_edges[eidx] << " ";

            total_num_valid_candidate_es += num_d2_edges[eidx];
        }
    }
    std::cout << "\n# candidates vertices: " << total_num_candidate_vs;
    std::cout << "\n# valid candidates vertices: " << total_num_valid_candidate_vs;
    std::cout << "\n# candidates edges: " << total_num_candidate_es;
    std::cout << "\n# valid candidates edges: " << total_num_valid_candidate_es << std::endl;
}
