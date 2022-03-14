#include <algorithm>
#include <cmath>
#include <vector>

#include "utils/types.h"
#include "graph/graph.h"
#include "graph/induced_graph.h"

InducedGraph::InducedGraph()
: graph()
, v_map_{}
, e_lists_{}
{}

InducedGraph::InducedGraph(const Graph& g, uint v1, uint v2)
: graph(&g)
, v_map_{std::min(v1, v2), std::max(v1, v2)}
, e_lists_{{std::min(v1, v2), std::max(v1, v2)}}
{}

InducedGraph::InducedGraph(InducedGraph& g1, InducedGraph& g2, bool is_union)
: graph(g1.graph)
, v_map_(g1.v_map_.size() + g2.v_map_.size())
, e_lists_(g1.e_lists_.size() + g2.e_lists_.size())
{
    if (is_union)
    {
        std::vector<uint>::iterator it1 = std::set_union(
            g1.v_map_.begin(), g1.v_map_.end(),
            g2.v_map_.begin(), g2.v_map_.end(),
            v_map_.begin()
        );
        v_map_.resize(it1 - v_map_.begin());

        std::vector<std::pair<uint, uint>>::iterator it2 = std::set_union(
            g1.e_lists_.begin(), g1.e_lists_.end(),
            g2.e_lists_.begin(), g2.e_lists_.end(),
            e_lists_.begin()
        );
        e_lists_.resize(it2 - e_lists_.begin());
    }
    else
    {
        std::vector<uint>::iterator it1 = std::set_intersection(
            g1.v_map_.begin(), g1.v_map_.end(),
            g2.v_map_.begin(), g2.v_map_.end(),
            v_map_.begin()
        );
        v_map_.resize(it1 - v_map_.begin());

        std::vector<std::pair<uint, uint>>::iterator it2 = std::set_intersection(
            g1.e_lists_.begin(), g1.e_lists_.end(),
            g2.e_lists_.begin(), g2.e_lists_.end(),
            e_lists_.begin()
        );
        e_lists_.resize(it2 - e_lists_.begin());
    }
}

uint InducedGraph::GetDegree(uint v)
{
    uint degree = 0u;

    for (auto& e: e_lists_)
        if (e.first == v || e.second == v)
            degree++;

    return degree;
}