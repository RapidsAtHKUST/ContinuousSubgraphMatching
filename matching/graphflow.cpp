#include <algorithm>
#include <iostream>
#include <vector>

#include "utils/types.h"
#include "utils/globals.h"
#include "utils/utils.h"
#include "graph/graph.h"
#include "matching/graphflow.h"

Graphflow::Graphflow(Graph& query_graph, Graph& data_graph, 
        uint max_num_results,
        bool print_prep, 
        bool print_enum, 
        bool homo)
: matching(
    query_graph, data_graph, max_num_results, 
    print_prep, print_enum, homo)
, order_vs_(query_.NumEdges())
, order_csrs_(query_.NumEdges())
, order_offs_(query_.NumEdges())
{
    for (uint i = 0; i < query_.NumEdges(); ++i)
    {
        order_vs_[i].resize(query_.NumVertices());
        order_csrs_[i].resize(query_.NumEdges() + 1);
        order_offs_[i].resize(query_.NumVertices(), 0);
    }
}

void Graphflow::Preprocessing()
{
    GenerateMatchingOrder();
}

void Graphflow::GenerateMatchingOrder()
{
    // generate the initial matching order, order_*s_[0]
    std::vector<bool> visited(query_.NumVertices(), false);
    uint max_degree = 0u;
    for (size_t i = 0; i < query_.NumVertices(); i++)
    {
        if (query_.GetDegree(i) > max_degree)
        {
            max_degree = query_.GetDegree(i);
            order_vs_[0][0] = i;
        }
    }
    visited[order_vs_[0][0]] = true;

    // loop over all remaining positions of the order
    for (uint i = 1; i < query_.NumVertices(); ++i)
    {
        uint max_adjacent = 0;
        uint max_adjacent_u = NOT_EXIST;
        for (size_t j = 0; j < query_.NumVertices(); j++)
        {
            uint cur_adjacent = 0u;
            if (visited[j]) continue;

            auto& q_nbrs = query_.GetNeighbors(j);
            for (auto& other : q_nbrs)
                if (visited[other])
                    cur_adjacent++;

            if (!cur_adjacent) continue;
            if (
                max_adjacent_u == NOT_EXIST ||
                (cur_adjacent == max_adjacent &&
                    query_.GetDegree(j) > query_.GetDegree(max_adjacent_u)) ||
                cur_adjacent > max_adjacent
            ) {
                max_adjacent = cur_adjacent;
                max_adjacent_u = j;
            }
        }
        order_vs_[0][i] = max_adjacent_u;
        visited[max_adjacent_u] = true;

        order_offs_[0][i] = order_offs_[0][i - 1];
        auto& q_nbrs = query_.GetNeighbors(max_adjacent_u);
        for (auto &other: q_nbrs)
            if (visited[other])
                order_csrs_[0][order_offs_[0][i]++] = other;
    }

    // generate other incremental matching orders
    for (uint i = 1; i < query_.NumEdges(); ++i)
    {
        std::vector<bool> visited(query_.NumVertices(), false);
        
        // get the first edge
        std::vector<uint>::iterator it = std::lower_bound(
            order_offs_[0].begin(), order_offs_[0].end(), i + 1
        );
        order_vs_[i][0] = *(order_vs_[0].begin() + std::distance(order_offs_[0].begin(), it));
        order_vs_[i][1] = order_csrs_[0][i];
        order_csrs_[i][0] = order_vs_[i][0];

        visited[order_vs_[i][0]] = true;
        visited[order_vs_[i][1]] = true;

        order_offs_[i][2] = order_offs_[i][1] = 1;
        for (uint j = 2; j < query_.NumVertices(); ++j)
        {
            uint max_adjacent = 0;
            uint max_adjacent_u = NOT_EXIST;
            for (size_t k = 0; k < query_.NumVertices(); k++)
            {
                uint cur_adjacent = 0u;
                if (visited[k]) continue;
                
                auto& q_nbrs = query_.GetNeighbors(k);
                for (auto& other : q_nbrs)
                    if (visited[other])
                        cur_adjacent++;

                if (!cur_adjacent) continue;
                if (
                    max_adjacent_u == NOT_EXIST ||
                    (cur_adjacent == max_adjacent &&
                        query_.GetDegree(k) > query_.GetDegree(max_adjacent_u)) ||
                    cur_adjacent > max_adjacent
                ) {
                    max_adjacent = cur_adjacent;
                    max_adjacent_u = k;
                }
            }
            order_vs_[i][j] = max_adjacent_u;
            visited[max_adjacent_u] = true;

            order_offs_[i][j] = order_offs_[i][j - 1];
            auto& q_nbrs = query_.GetNeighbors(max_adjacent_u);
            for (auto &other: q_nbrs)
                if (visited[other])
                    order_csrs_[i][order_offs_[i][j]++] = other;
        }
    }
    if (print_preprocessing_results_)
    {
        std::cout << "matching order: " << std::endl;
        std::cout << "-vertex(backward neighbors)-\n";
        for (uint i = 0; i < query_.NumEdges(); ++i)
        {
            std::cout << "#" << i << ": ";
            for (uint j = 0; j < query_.NumVertices(); ++j)
            {
                std::cout << order_vs_[i][j];
                if (j == 0)
                {
                    std::cout << "-";
                    continue;
                }

                for (uint k = order_offs_[i][j - 1]; k < order_offs_[i][j]; k++)
                {
                    if (k == order_offs_[i][j - 1]) std::cout << "(";
                    std::cout << order_csrs_[i][k];
                    if (k != order_offs_[i][j] - 1) std::cout << ",";
                    else std::cout << ")";
                }
                if (j != query_.NumVertices() - 1)
                    std::cout << "-";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void Graphflow::FindMatches(uint order_index, uint depth, std::vector<uint> m, size_t &num_results)
{
    if (reach_time_limit) return;

    uint u = order_vs_[order_index][depth];
    uint u_min = NOT_EXIST;
    uint u_min_label = NOT_EXIST;
    uint u_min_size = UINT_MAX;

    // find u_min
    const auto& q_nbrs = query_.GetNeighbors(u);
    const auto& q_nbr_labels = query_.GetNeighborLabels(u);

    for (uint i = 0u; i < q_nbrs.size(); i++)
    {
        const uint u_other = q_nbrs[i];
        const uint u_other_label = q_nbr_labels[i];

        if (m[u_other] == UNMATCHED) continue;

        const uint cur_can_size = data_.GetNeighbors(m[u_other]).size();
        if (cur_can_size < u_min_size)
        {
            u_min_size = cur_can_size;
            u_min = u_other;
            u_min_label = u_other_label;
        }
    }

    const auto& u_min_nbrs = data_.GetNeighbors(m[u_min]);
    const auto& u_min_nbr_labels = data_.GetNeighborLabels(m[u_min]);

    bool candidate_empty = true;
    for (uint i = 0u; i < u_min_nbrs.size(); i++)
    {
        const uint v = u_min_nbrs[i];

        // 1. check labels
        num_intermediate_results_before_index_check_++;
        if (
            data_.GetVertexLabel(v) != query_.GetVertexLabel(u) ||
            u_min_nbr_labels[i] != u_min_label
        ) continue;
        num_intermediate_results_after_index_check_++;

        // 2. check if joinable
        bool joinable = true;
        for (uint j = 0u; j < q_nbrs.size(); j++)
        {
            const uint u_other = q_nbrs[j];
            const uint u_other_labels = q_nbr_labels[j];

            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(data_.GetNeighbors(m[u_other]).begin(), data_.GetNeighbors(m[u_other]).end(), v);
            uint dis = std::distance(data_.GetNeighbors(m[u_other]).begin(), it);
            if (
                it == data_.GetNeighbors(m[u_other]).end() ||
                *it != v ||
                data_.GetNeighborLabels(m[u_other])[dis] != u_other_labels
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

void Graphflow::InitialMatching()
{
    std::vector<uint> m(query_.NumVertices(), UNMATCHED);
    
    for (size_t i = 0; i < data_.NumVertices(); i++)
    if (data_.GetVertexLabel(i) != NOT_EXIST)
    if (query_.GetVertexLabel(order_vs_[0][0]) == data_.GetVertexLabel(i))
    {
        m[order_vs_[0][0]] = i;
        visited_[i] = true;

        FindMatches(0, 1, m, num_initial_results_);

        visited_[i] = false;
        m[order_vs_[0][0]] = UNMATCHED;
    }
}

void Graphflow::AddEdge(uint v1, uint v2, uint label)
{
    data_.AddEdge(v1, v2, label);

    std::vector<uint> m(query_.NumVertices(), UNMATCHED);

    if (max_num_results_ == 0) return;

    size_t num_results = 0ul;
    for (uint i = 0; i < query_.NumEdges(); i++)
    {
        uint u1 = order_vs_[i][0], u2 = order_vs_[i][1];
        auto temp_q_labels = query_.GetEdgeLabel(u1, u2);
        
        // check if any query edge match (v1 --> v2)
        if (
            std::get<0>(temp_q_labels) == data_.GetVertexLabel(v1) &&
            std::get<1>(temp_q_labels) == data_.GetVertexLabel(v2) &&
            std::get<2>(temp_q_labels) == label
        ) {
            m[u1] = v1;
            m[u2] = v2;
            visited_[v1] = true;
            visited_[v2] = true;

            FindMatches(i, 2, m, num_results);

            visited_[v1] = false;
            visited_[v2] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;

            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) goto END_ENUMERATION;
        }
        // check if any query edge match (v2 --> v1)
        if (
            std::get<0>(temp_q_labels) == data_.GetVertexLabel(v2) &&
            std::get<1>(temp_q_labels) == data_.GetVertexLabel(v1) &&
            std::get<2>(temp_q_labels) == label
        ) {
            m[u1] = v2;
            m[u2] = v1;
            visited_[v2] = true;
            visited_[v1] = true;

            FindMatches(i, 2, m, num_results);

            visited_[v2] = false;
            visited_[v1] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;

            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) goto END_ENUMERATION;
        }
    }
    END_ENUMERATION:
    num_positive_results_ += num_results;
}

void Graphflow::RemoveEdge(uint v1, uint v2)
{
    std::vector<uint> m(query_.NumVertices(), UNMATCHED);

    std::tuple labels = data_.GetEdgeLabel(v1, v2);
    
    size_t num_results = 0ul;
    if (max_num_results_ == 0) goto END_ENUMERATION;

    for (uint i = 0; i < query_.NumEdges(); i++)
    {
        uint u1 = order_vs_[i][1], u2 = order_csrs_[i][0];
        auto temp_q_labels = query_.GetEdgeLabel(u1, u2);

        if (
            std::get<0>(temp_q_labels) == std::get<0>(labels) &&
            std::get<1>(temp_q_labels) == std::get<1>(labels) &&
            std::get<2>(temp_q_labels) == std::get<2>(labels)
        ) {
            m[u1] = v1;
            m[u2] = v2;
            visited_[v1] = true;
            visited_[v2] = true;

            FindMatches(i, 2, m, num_results);

            visited_[v1] = false;
            visited_[v2] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;

            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) goto END_ENUMERATION;
        }
        if (
            std::get<1>(temp_q_labels) == std::get<0>(labels) &&
            std::get<0>(temp_q_labels) == std::get<1>(labels) &&
            std::get<2>(temp_q_labels) == std::get<2>(labels)
        ) {
            m[u1] = v2;
            m[u2] = v1;
            visited_[v2] = true;
            visited_[v1] = true;

            FindMatches(i, 2, m, num_results);

            visited_[v2] = false;
            visited_[v1] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;

            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) goto END_ENUMERATION;
        }
    }
    END_ENUMERATION:

    num_negative_results_ += num_results;
    data_.RemoveEdge(v1, v2);
}

void Graphflow::AddVertex(uint id, uint label)
{
    data_.AddVertex(id, label);
    
    visited_.resize(id + 1, false);
}

void Graphflow::RemoveVertex(uint id)
{
    data_.RemoveVertex(id);
}

void Graphflow::GetMemoryCost(size_t &num_edges, size_t &num_vertices)
{
    num_edges = 0ul;
    num_vertices = 0ul;
}
