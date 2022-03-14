#include <algorithm>
#include <array>
#include <iostream>
#include <unordered_set>
#include <vector>

#include "utils/types.h"
#include "utils/globals.h"
#include "graph/induced_graph.h"
#include "matching/sj_tree.h"

SJTree::SJTree(Graph& query_graph, Graph& data_graph, 
        uint max_num_results,
        bool print_prep,
        bool print_enum,
        bool homo)
: matching(
    query_graph, data_graph, max_num_results, 
    print_prep, print_enum, homo)
, order_es_(query_.NumEdges())
, treeNodes_(query_.NumEdges() * 2 - 1)
{}

void SJTree::Preprocessing()
{
    GenerateMatchingOrder();
    BuildSJTree();
}

void SJTree::GenerateMatchingOrder()
{
    struct EdgeFreq{
        uint v1; uint v2; uint v1_label;
        uint v2_label; uint e_label; uint freq;
        EdgeFreq(
            uint v1_arg, uint v2_arg, uint v1_label_arg, 
            uint v2_label_arg, uint e_label_arg, uint freq_arg)
        : v1(v1_arg), v2(v2_arg), v1_label(v1_label_arg)
        , v2_label(v2_label_arg), e_label(e_label_arg), freq(freq_arg)
        {}
    };

    // get the frequency of each query edge
    std::vector<EdgeFreq> edge_freqs;

    for (size_t u = 0; u < query_.NumVertices(); u++)
    {
        const auto& q_nbrs = query_.GetNeighbors(u);
        const auto& q_nbr_labels = query_.GetNeighborLabels(u);

        for (uint i = 0; i < q_nbrs.size(); i++)
        {
            const uint u_other = q_nbrs[i];
            const uint u_other_label = q_nbr_labels[i];
            
            if (u > u_other) continue;

            uint freq = 0u;

            for (size_t v = 0; v < data_.NumVertices(); v++)
            if (data_.GetVertexLabel(v) != NOT_EXIST)
            {
                const auto& d_nbrs = data_.GetNeighbors(v);
                const auto& d_nbr_labels = data_.GetNeighborLabels(v);

                for (uint j = 0; j < d_nbrs.size(); j++)
                {
                    const uint v_other = d_nbrs[j];
                    const uint v_other_label = d_nbr_labels[j];
                    
                    if (
                        data_.GetVertexLabel(v) == query_.GetVertexLabel(u) &&
                        data_.GetVertexLabel(v_other) == query_.GetVertexLabel(u_other) &&
                        v_other_label == u_other_label
                    ) {
                        freq ++;
                    }
                }
            }

            edge_freqs.emplace_back(
                u, u_other, query_.GetVertexLabel(u),
                query_.GetVertexLabel(u_other), u_other_label, freq
            );
        }
    }
    
    // find edge with min frequency
    auto it = std::min_element(edge_freqs.begin(), edge_freqs.end(),
    [](EdgeFreq& t1, EdgeFreq& t2) {
        return t1.freq < t2.freq;
    });

    order_es_[0] = {it->v1, it->v2, it->e_label};
    std::unordered_set<uint> visited;
    visited.insert(it->v1);
    visited.insert(it->v2);

    edge_freqs.erase(it);

    for (uint i = 1; i < query_.NumEdges(); i++)
    {
        uint min_freq = UINT_MAX;
        uint min_index = 0;
        for (uint j = 0; j < edge_freqs.size(); j++)
        {
            auto& edge_freq = edge_freqs[j];
            if (
                edge_freq.freq < min_freq &&
                (visited.find(edge_freq.v1) != visited.end() || 
                visited.find(edge_freq.v2) != visited.end())
            ) {
                min_freq = edge_freq.freq;
                min_index = j;
            }
        }
        order_es_[i] = {edge_freqs[min_index].v1, edge_freqs[min_index].v2, edge_freqs[min_index].e_label};
        visited.insert(edge_freqs[min_index].v1);
        visited.insert(edge_freqs[min_index].v2);
        
        edge_freqs.erase(edge_freqs.begin() + min_index);
    }
}

void SJTree::BuildSJTree()
{
    for (uint i = 0; i < query_.NumEdges(); i++)
    {
        if (i == 0)
        {
            treeNodes_[0].graph_ = InducedGraph(query_, order_es_[i].v1, order_es_[i].v2);
        }
        else
        {
            uint cur_leaf_id = 2 * i;
            uint parent_id = cur_leaf_id - 1;
            uint sibling_id = cur_leaf_id == 2 ? 0 : cur_leaf_id - 3;

            treeNodes_[cur_leaf_id].parent_ = parent_id;
            treeNodes_[cur_leaf_id].sibling_ = sibling_id;
            treeNodes_[cur_leaf_id].graph_ = InducedGraph(query_, order_es_[i].v1, order_es_[i].v2);

            treeNodes_[sibling_id].parent_ = parent_id;
            treeNodes_[sibling_id].sibling_ = cur_leaf_id;

            treeNodes_[parent_id].intersection_ = InducedGraph(treeNodes_[sibling_id].graph_, treeNodes_[cur_leaf_id].graph_, false);
            treeNodes_[parent_id].graph_ = InducedGraph(treeNodes_[sibling_id].graph_, treeNodes_[cur_leaf_id].graph_, true);
        }
    }

    if (print_preprocessing_results_)
    {
        std::cout << "SJ-Tree Nodes: " << std::endl;
        for (uint i = 0; i < query_.NumEdges() * 2 - 1; ++i)
        {
            std::cout << "#" << i << ": v:";
            for (uint j = 0; j < treeNodes_[i].graph_.v_map_.size(); j++)
            {
                std::cout << treeNodes_[i].graph_.v_map_[j];
                if (j != treeNodes_[i].graph_.v_map_.size() - 1)
                    std::cout << ",";
            }
            std::cout << " e:";
            for (uint j = 0; j < treeNodes_[i].graph_.e_lists_.size(); j++)
            {
                std::cout << treeNodes_[i].graph_.e_lists_[j].first << 
                    "-" << treeNodes_[i].graph_.e_lists_[j].second;
                if (j != treeNodes_[i].graph_.e_lists_.size() - 1)
                    std::cout << ",";
            }
            for (uint j = 0; j < treeNodes_[i].intersection_.v_map_.size(); j++)
            {
                if (j == 0) std::cout << " keys: ";
                std::cout << treeNodes_[i].intersection_.v_map_[j];
                if (j != treeNodes_[i].intersection_.v_map_.size() - 1)
                    std::cout << ",";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void SJTree::AddSingleMatch(const std::vector<uint>& match, 
        TreeNode &node, size_t& num_results)
{
    if (reach_time_limit) return;
    
    uint parent_id = node.parent_;
    if (parent_id == UINT_MAX)
    {
        // node.matches[0][0].emplace_back(match);
        num_results++;
        if (print_enumeration_results_)
        {
            for (auto j: match)
                std::cout << j << " ";
            std::cout << std::endl;
        }
        return;
    }
    InducedGraph& intersect = treeNodes_[parent_id].intersection_;
    std::array<uint, 2> key {};
    
    // recursive join
    if (intersect.v_map_.size() > 1)
    {
        // The intersection graph contains two vertex.
        // This means one child node represent a query edge,
        // and the other child node contains all query vertices
        // in the parent node.

        // get the two join keys
        for (uint i = 0; i < node.graph_.v_map_.size(); i++)
        {
            if (node.graph_.v_map_[i] == intersect.v_map_[0])
                key[0] = match[i];

            if (node.graph_.v_map_[i] == intersect.v_map_[1])
                key[1] = match[i];
        }
        for (auto& m: match)
        {
            node.matches_[{key[0], key[1]}].emplace_back(m);
        }

        auto& s_matches = treeNodes_[node.sibling_].matches_;
        if (s_matches.find({key[0], key[1]}) != s_matches.end())
        {
            const uint sibling_size = treeNodes_[node.sibling_].graph_.v_map_.size();
            // if the sibling node represent an edge, then the parent match
            // is generated from current match
            if (sibling_size == 2)
            {
                AddSingleMatch(match, treeNodes_[parent_id], num_results);
                if (reach_time_limit) return;
            }
            // if the current node represent an edge, then each parent match
            // is generated from each match on the sibling node
            else
            {
                uint num_matches = s_matches.at({key[0], key[1]}).size() / sibling_size;
                for (uint i = 0; i < num_matches; i++)
                {
                    std::vector<uint> s_match(
                        s_matches[{key[0], key[1]}].begin() + i * sibling_size,
                        s_matches[{key[0], key[1]}].begin() + (i + 1) * sibling_size);
                    AddSingleMatch(s_match, treeNodes_[parent_id], num_results);
                    if (reach_time_limit) return;
                }
            }
        }
    }
    else
    {
        // The intersection graph contains one vertex.

        // get the join key
        for (uint i = 0; i < node.graph_.v_map_.size(); i++)
        {
            if (node.graph_.v_map_[i] == intersect.v_map_.front())
                key[0] = match[i];
        }
        for (auto& m: match)
        {
            node.matches_[{key[0], 0}].emplace_back(m);
        }

        // find the matches on the sibling node with the same join key
        auto& s_matches = treeNodes_[node.sibling_].matches_;
        if (s_matches.find({key[0], 0}) != s_matches.end())
        {
            const uint sibling_size = treeNodes_[node.sibling_].graph_.v_map_.size();

            uint num_matches = s_matches.at({key[0], key[1]}).size() / sibling_size;
            for (uint i = 0; i < num_matches; i++)
            {
                std::vector<uint> s_match(
                    s_matches[{key[0], key[1]}].begin() + i * sibling_size,
                    s_matches[{key[0], key[1]}].begin() + (i + 1) * sibling_size);

                // build a match after join
                uint p_size = sibling_size + match.size() - 1;
                std::vector<uint> p_match(p_size);
                uint cur_pos = 0, s_pos = 0, write_pos = 0;
                auto& s_vs = treeNodes_[node.sibling_].graph_.v_map_;
                auto& cur_vs = node.graph_.v_map_;
                
                while (write_pos < p_size)
                {
                    if (s_pos >= s_match.size())
                    {
                        p_match[write_pos] = match[cur_pos++];
                    }
                    else if (cur_pos >= match.size())
                    {
                        p_match[write_pos] = s_match[s_pos++];
                    }
                    else if (cur_vs[cur_pos] < s_vs[s_pos])
                    {
                        p_match[write_pos] = match[cur_pos++];
                    }
                    else if (cur_vs[cur_pos] > s_vs[s_pos])
                    {
                        p_match[write_pos] = s_match[s_pos++];
                    }
                    else
                    {
                        p_match[write_pos] = match[cur_pos++];
                        s_pos ++;
                    }
                    write_pos ++;
                }
                if (homomorphism_)
                {
                    AddSingleMatch(p_match, treeNodes_[parent_id], num_results);
                }
                else
                {
                    // ensure no two query vertices map the same data vertex
                    std::vector<uint> sorted_p_match(p_match);
                    std::sort(sorted_p_match.begin(), sorted_p_match.end());
                    if (std::unique(sorted_p_match.begin(), sorted_p_match.end()) == sorted_p_match.end())
                        AddSingleMatch(p_match, treeNodes_[parent_id], num_results);
                }
                if (reach_time_limit) return;
            }
        }
    }
}

void SJTree::InitialMatching()
{
    for (uint i = 0; i < data_.NumVertices(); i++)
    if (data_.GetVertexLabel(i) != NOT_EXIST)
    {
        auto& d_nbrs = data_.GetNeighbors(i);
        auto& d_nbr_labels = data_.GetNeighborLabels(i);

        for (uint j = 0; j < d_nbrs.size(); j++)
        {
            const uint v_other = d_nbrs[j];
            const uint v_other_label = d_nbr_labels[j];

            for (uint k = 0; k < order_es_.size(); k++)
            {
                if (
                    data_.GetVertexLabel(v_other) == query_.GetVertexLabel(order_es_[k].v1) &&
                    data_.GetVertexLabel(i) == query_.GetVertexLabel(order_es_[k].v2) &&
                    v_other_label == order_es_[k].e_label
                ) {
                    std::vector<uint> temp_match {v_other, i};
                    AddSingleMatch(temp_match, treeNodes_[k * 2], num_initial_results_);
                }
            }
        }
    }
}

void SJTree::AddEdge(uint v1, uint v2, uint label)
{
    data_.AddEdge(v1, v2, label);

    for (uint k = 0; k < order_es_.size(); k++)
    {
        // check if any query edge match (v1 --> v2)
        if (
            data_.GetVertexLabel(v1) == query_.GetVertexLabel(order_es_[k].v1) &&
            data_.GetVertexLabel(v2) == query_.GetVertexLabel(order_es_[k].v2) &&
            label == order_es_[k].e_label
        ) {
            std::vector<uint> temp_match {v1, v2};
            AddSingleMatch(temp_match, treeNodes_[k * 2], num_positive_results_);
            if (reach_time_limit) return;
        }
        // check if any query edge match (v2 --> v1)
        if (
            data_.GetVertexLabel(v2) == query_.GetVertexLabel(order_es_[k].v1) &&
            data_.GetVertexLabel(v1) == query_.GetVertexLabel(order_es_[k].v2) &&
            label == order_es_[k].e_label
        ) {
            std::vector<uint> temp_match {v2, v1};
            AddSingleMatch(temp_match, treeNodes_[k * 2], num_positive_results_);
            if (reach_time_limit) return;
        }
    }
}

void SJTree::RemoveEdge(uint v1, uint v2)
{
    data_.RemoveEdge(v1, v2);
}

void SJTree::AddVertex(uint id, uint label)
{
    data_.AddVertex(id, label);

    visited_.resize(id + 1, false);
}

void SJTree::RemoveVertex(uint id)
{
    data_.RemoveVertex(id);
}

void SJTree::GetMemoryCost(size_t &num_edges, size_t &num_vertices)
{
    num_edges = 0ul;
    num_vertices = 0ul;
    for (uint i = 0; i < 2 * query_.NumEdges() - 1; i++)
    {
        for (const auto& [k1, v]: treeNodes_[i].matches_)
        {
            num_vertices += v.size();
            num_edges += v.size() / 
                treeNodes_[i].graph_.NumVertices() * 
                treeNodes_[i].graph_.NumEdges();
        }
    }
}
