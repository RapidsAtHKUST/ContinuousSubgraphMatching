#include <algorithm>
#include <iostream>
#include <queue>
#include <stack>
#include <unordered_set>
#include <vector>

#include "utils/types.h"
#include "utils/globals.h"
#include "graph/graph.h"
#include "matching/iedyn.h"

IEDyn::IEDyn(Graph& query_graph, Graph& data_graph, 
        uint max_num_results,
        bool print_prep, 
        bool print_enum, 
        bool homo)
: matching(
    query_graph, data_graph, max_num_results, 
    print_prep, print_enum, homo)
, eidx_(query_.NumVertices())
, treeNode_(query_.NumVertices())
, q_root_(0u)
, serialized_tree_(query_.NumVertices())

, DCS_(query_.NumEdges() * 2)
, DCS_CD_(query_.NumEdges() * 2)
, DCS_CD_delta_(query_.NumEdges() * 2)
, root_d2_indices_{}
, root_d2_indices_delta_{}

, is_delta_node_(query_.NumEdges() * 2)
, match_node_(UINT_MAX)

, d2(query_.NumVertices())
, n2(query_.NumEdges() * 2)
, nc2(query_.NumVertices())
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
}

void IEDyn::Preprocessing()
{
    BuildDAG();
    BuildDCS();
    GenerateMatchingOrder();
}

void IEDyn::BuildDAG()
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

void IEDyn::BuildDCS()
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
                    }
                }
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
                        {
                            n2[eidx_[u][u_other]][j] += 1;
                            DCS_CD_[eidx_[u][u_other]][j].push_back(v_other);
                        }
                    }
                }
                if (n2[eidx_[u][u_other]][j])
                {
                    nc2[u][j] += 1;
                }
            }
            if (nc2[u][j] == treeNode_[u].forwards_.size())
            {
                d2[u][j] = 1;
                if (u == q_root_)
                    root_d2_indices_.insert(j);
            }
        }
    }

}

void IEDyn::GenerateMatchingOrder()
{
    std::stack<uint> stack;
    std::vector<uint> progress (query_.NumVertices(), 0);
    uint cur = q_root_;
    while (cur != UINT_MAX || !stack.empty())
    {
        while (cur != UINT_MAX)
        {
            // process the current node
            order_vs_.push_back(cur);
            if (!treeNode_[cur].backwards_.empty())
                backward_vs_.push_back(treeNode_[cur].backwards_[0]);
            else
                backward_vs_.push_back(UINT_MAX);
            stack.push(cur);
            //progress[cur] += 1;
            if (treeNode_[cur].forwards_.size() > 0)
                cur = treeNode_[cur].forwards_[0];
            else
                cur = UINT_MAX;
        }
        cur = stack.top();
        progress[cur] += 1;
        if (progress[cur] < treeNode_[cur].forwards_.size())
        {
            cur = treeNode_[cur].forwards_[progress[cur]];
        }
        else
        {
            stack.pop();
            cur = UINT_MAX;
        }
    }
    std::cout << "matching order:" << std::endl;
    std::cout << "-vertex(backward neighbors)-" << std::endl;
    for (size_t i = 0ul; i < query_.NumVertices(); i++)
    {
        if (i == 0ul)
            std::cout << order_vs_[i];
        else
            std::cout << order_vs_[i] << "(" << backward_vs_[i] << ")";
        if (i < query_.NumVertices() - 1)
            std::cout << "-";
    }
    std::cout << std::endl;
}

void IEDyn::InitialMatching()
{
    std::vector<uint> m (query_.NumVertices(), UNMATCHED);
    
    for (auto & v: root_d2_indices_)
    {
        m[q_root_] = v;
        visited_[v] = true;

        FindMatchesInitial(1, m, num_initial_results_);

        visited_[v] = false;
        m[q_root_] = UNMATCHED;
    }
}

void IEDyn::InsertionBottomUp(uint u, uint u_p, uint v, uint v_p)
{
    if (n2[eidx_[u_p][u]][v_p] == 0)
    {
        nc2[u_p][v_p] += 1;
    }
    n2[eidx_[u_p][u]][v_p] += 1;

    if (nc2[u_p][v_p] == treeNode_[u_p].forwards_.size())
    {
        is_delta_node_[eidx_[u_p][u]] = true;
        d2[u_p][v_p] = 1;
        if (u_p != q_root_)
        {
            DCS_CD_delta_[eidx_[u_p][treeNode_[u_p].backwards_[0]]][v_p] = 
                DCS_[eidx_[u_p][treeNode_[u_p].backwards_[0]]][v_p];
            for (auto& nbr: DCS_[eidx_[u_p][treeNode_[u_p].backwards_[0]]][v_p])
            {
                auto& temp_nbrs = DCS_CD_delta_[eidx_[treeNode_[u_p].backwards_[0]][u_p]][nbr];
                auto it = std::lower_bound(temp_nbrs.begin(), temp_nbrs.end(), v_p);
                temp_nbrs.insert(it, v_p);
            }
        }
        else
        {
            root_d2_indices_delta_.insert(v_p);
        }
        
        Q2.emplace(u_p, v_p);
    }
}

void IEDyn::FindMatchesInitial(uint depth, std::vector<uint>& m, size_t &num_results)
{
    if (reach_time_limit) return;

    uint u = order_vs_[depth];
    uint u_min = backward_vs_[depth];

    // enumerate each neighbor of m[u_min]
    bool candidate_empty = true;

    for (auto& v: DCS_CD_[eidx_[u_min][u]].at(m[u_min]))
    {
        // 1. check index
        num_intermediate_results_before_index_check_++;
        if (d2[u][v] == 0) continue;
        num_intermediate_results_after_index_check_++;

        // 2. check if joinable, skipped
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
            FindMatchesInitial(depth + 1, m, num_results);
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

void IEDyn::FindMatchesIncremental(uint depth, std::vector<uint>& m, size_t &num_results)
{
    if (reach_time_limit) return;

    uint u = order_vs_[depth];
    uint u_min = backward_vs_[depth];

    // enumerate each neighbor of m[u_min]
    bool candidate_empty = true;

    auto& DCS_CD = is_delta_node_[eidx_[u_min][u]] ? 
        DCS_CD_delta_ :
        DCS_CD_;
    for (auto& v: DCS_CD[eidx_[u_min][u]].at(m[u_min]))
    {
        // 1. check index
        if (u != match_node_) num_intermediate_results_before_index_check_++;
        if (d2[u][v] == 0) continue;
        if (u != match_node_) num_intermediate_results_after_index_check_++;

        // 2. check if joinable, skipped
        if (u != match_node_) num_intermediate_results_after_joinability_check_++;

        candidate_empty = false;

        // 3. check if visited
        if (!homomorphism_ && visited_[v]) continue;
        if (u != match_node_) num_intermediate_results_after_visit_check_++;
        
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
            if (u != match_node_) 
            {
                size_t num_results_before_recursion = num_results;
                FindMatchesIncremental(depth + 1, m, num_results);
                if (num_results == num_results_before_recursion)
                {
                    num_intermediate_results_without_results_++;
                }
            }
            else
            {
                FindMatchesIncremental(depth + 1, m, num_results);
            }
        }
        
        visited_[v] = false;
        m[u] = UNMATCHED;
        if (num_results >= max_num_results_) return;
        if (reach_time_limit) return;
    }
    if (candidate_empty) num_intermediate_results_with_empty_candidate_set_++;
}

void IEDyn::CombineDeltaDCS()
{
    for (size_t i = 0ul; i < DCS_CD_.size(); i++)
    {
        for (auto& [v, nbrs]: DCS_CD_delta_[i])
        {
            if (DCS_CD_[i].find(v) == DCS_CD_[i].end())
            {
                DCS_CD_[i][v] = nbrs;
            }
            else
            {
                for (const uint& nbr: nbrs)
                {
                    auto it = std::lower_bound(
                        DCS_CD_[i][v].begin(),
                        DCS_CD_[i][v].end(),
                        nbr
                    );
                    if (it == DCS_CD_[i][v].end() || *it != nbr)
                    {
                        DCS_CD_[i][v].insert(it, nbr);
                    }
                }
            }
        }
        DCS_CD_delta_[i].clear();
        is_delta_node_[i] = false;
    }
    for (size_t i: root_d2_indices_delta_)
    {
        root_d2_indices_.insert(i);
    }
    root_d2_indices_delta_.clear();
}

void IEDyn::AddEdge(uint v1, uint v2, uint label)
{
    data_.AddEdge(v1, v2, label);

    size_t num_results = 0ul;
    std::vector<uint> m(query_.NumVertices(), UNMATCHED);

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

            is_delta_node_[eidx_[u1][u2]] = true;
            match_node_ = u2;

            bool old_c_d2 = d2[u2][v2];

            if (old_c_d2)
            {
                it = std::lower_bound(DCS_CD_delta_[eidx_[u1][u2]][v1].begin(), DCS_CD_delta_[eidx_[u1][u2]][v1].end(), v2);
                DCS_CD_delta_[eidx_[u1][u2]][v1].insert(it, v2);
                it = std::lower_bound(DCS_CD_delta_[eidx_[u2][u1]][v2].begin(), DCS_CD_delta_[eidx_[u2][u1]][v2].end(), v1);
                DCS_CD_delta_[eidx_[u2][u1]][v2].insert(it, v1);

                InsertionBottomUp(u2, u1, v2, v1);
            }

            while (!Q2.empty())
            {
                auto [u_queue, v_queue] = Q2.front();
                Q2.pop();
                for (auto& u_p_queue : treeNode_[u_queue].backwards_)
                for (auto& v_p_queue : DCS_CD_delta_[eidx_[u_queue][u_p_queue]][v_queue])
                {
                    InsertionBottomUp(u_queue, u_p_queue, v_queue, v_p_queue);
                    if (reach_time_limit) return;
                }
            }

            if (num_results >= max_num_results_) goto END_CURRENT_ENUMERATION1;

            for (auto & v: root_d2_indices_delta_)
            {
                m[q_root_] = v;
                visited_[v] = true;

                FindMatchesIncremental(1, m, num_results);

                visited_[v] = false;
                m[q_root_] = UNMATCHED;
            }
            
            END_CURRENT_ENUMERATION1:
            if (reach_time_limit) goto END_ENUMERATION;
            CombineDeltaDCS();
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

void IEDyn::RemoveEdge(uint v1, uint v2)
{
    data_.RemoveEdge(v1, v2);
}

void IEDyn::AddVertex(uint id, uint label)
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
            if (treeNode_[u].forwards_.empty())
                d2[u][id] = 1;
        }
    }
}

void IEDyn::RemoveVertex(uint id)
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
                if (n2[eidx_[u][q_nbrs[k]]].find(id) != n2[eidx_[u][q_nbrs[k]]].end())
                    n2[eidx_[u][q_nbrs[k]]].erase(id);
            }
            if (nc2[u].find(id) != nc2[u].end())
                nc2[u].erase(id);
            if (d2[u].find(id) != d2[u].end())
                d2[u].erase(id);
        }
    }

    data_.RemoveVertex(id);
}

void IEDyn::GetMemoryCost(size_t &num_edges, size_t &num_vertices)
{
    num_edges = 0ul;
    num_vertices = 0ul;

    std::vector<std::unordered_set<uint>> all_vertices_in_index(query_.NumVertices());
    std::vector<std::unordered_set<uint>> d2_vertices_in_index(query_.NumVertices());

    std::vector<size_t> num_all_edges(query_.NumEdges() * 2, 0);
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
                    if (d2[i][v])
                    {
                        d2_vertices_in_index[i].insert(v);
                    }

                    for (uint j = 0; j < d_nbrs.size(); j++)
                    {
                        num_all_edges[eidx] += 1;
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
                    if (d2[i][v])
                        d2_vertices_in_index[i].insert(v);
                }
            }
        }
    }

    size_t total_num_valid_candidate_vs = 0ul;
    size_t total_num_valid_candidate_es = 0ul;

    std::cout << "# vertices in index: "; 
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        std::cout << i << ": " << all_vertices_in_index[i].size() << " ";
        num_vertices += all_vertices_in_index[i].size();
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
    std::cout << "\n# valid candidates vertices: " << total_num_valid_candidate_vs;
    std::cout << "\n# valid candidates edges: " << total_num_valid_candidate_es << std::endl;
}
