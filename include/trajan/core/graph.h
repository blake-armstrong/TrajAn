#pragma once
#include <ankerl/unordered_dense.h>
#include <fmt/core.h>
#include <queue>
#include <trajan/core/log.h>
#include <vector>

namespace trajan::core::graph {

template <typename NodeType, typename EdgeType> class ConnectedComponent;

template <typename NodeType, typename EdgeType> class Graph {
public:
  using NodeID = size_t;
  using ConnectedComponent = ConnectedComponent<NodeType, EdgeType>;
  using AdjacencyList = ankerl::unordered_dense::map<
      NodeID, ankerl::unordered_dense::map<NodeID, EdgeType>>;

  Graph() = default;
  Graph(const std::vector<NodeType> &nodes) : m_nodes(nodes) {}

  inline const std::vector<NodeType> &nodes() const { return m_nodes; }

  inline size_t num_nodes() const { return m_nodes.size(); }

  inline std::vector<ConnectedComponent> find_connected_components() const {
    std::vector<ConnectedComponent> components;
    components.clear();
    components.reserve(m_nodes.size());
    if (m_adjacency_list.empty()) {
      trajan::log::critical("Graph structure not built.");
      return components;
    }
    ankerl::unordered_dense::set<NodeID> visited;

    for (size_t i = 0; i < m_nodes.size(); i++) {
      NodeID node_id = get_node_id(i);

      if (visited.contains(node_id)) {
        continue;
      }

      ConnectedComponent component;

      std::queue<NodeID> queue;
      queue.push(node_id);
      visited.insert(node_id);

      while (!queue.empty()) {
        NodeID current_id = queue.front();
        queue.pop();
        component.add_node(current_id, get_node_by_id(current_id));

        if (!m_adjacency_list.contains(current_id)) {
          continue;
        }
        for (const auto &[neighbour_id, edge_data] :
             m_adjacency_list.at(current_id)) {
          component.add_edge(current_id, neighbour_id, edge_data);

          if (visited.contains(neighbour_id)) {
            continue;
          }
          queue.push(neighbour_id);
          visited.insert(neighbour_id);
        }
      }

      components.push_back(component);
    }

    return components;
  }

  inline const NodeType &get_node_by_id(NodeID id) const {
    for (const auto &node : m_nodes) {
      if (get_node_id_from_node(node) == id) {
        return node;
      }
    }
    throw std::runtime_error("Node ID not found");
  }

  inline const AdjacencyList &get_adjacency_list() const {
    return m_adjacency_list;
  }

  virtual NodeID get_node_id_from_node(const NodeType &node) const = 0;

  inline void add_edge(NodeID idx1, NodeID idx2, EdgeType &edge) {
    m_adjacency_list[idx1][idx2] = edge;
    m_adjacency_list[idx2][idx1] = edge;
  }
  inline void remove_edge(NodeID idx1, NodeID idx2) {
    auto source_it1 = m_adjacency_list.find(idx1);
    if (source_it1 != m_adjacency_list.end()) {
      source_it1->second.erase(idx2);
      if (source_it1->second.empty()) {
        m_adjacency_list.erase(source_it1);
      }
    }
    auto source_it2 = m_adjacency_list.find(idx2);
    if (source_it2 != m_adjacency_list.end()) {
      source_it2->second.erase(idx1);
      if (source_it2->second.empty()) {
        m_adjacency_list.erase(source_it2);
      }
    }
  }
  inline void clear_edges() { m_adjacency_list.clear(); }

protected:
  std::vector<NodeType> m_nodes;
  AdjacencyList m_adjacency_list;

  NodeID get_node_id(size_t index) const {
    return get_node_id_from_node(m_nodes[index]);
  }
};

struct PairHash {
  template <class T1, class T2>
  size_t operator()(const std::pair<T1, T2> &p) const {
    auto h1 = std::hash<T1>{}(p.first);
    auto h2 = std::hash<T2>{}(p.second);
    return h1 ^ (h2 << 1);
  }
};

template <typename NodeType, typename EdgeType> class ConnectedComponent {
public:
  using NodeID = typename Graph<NodeType, EdgeType>::NodeID;
  using Nodes = ankerl::unordered_dense::map<NodeID, NodeType>;
  using Edges = ankerl::unordered_dense::map<std::pair<NodeID, NodeID>,
                                             EdgeType, PairHash>;

  inline void add_node(NodeID id, const NodeType &node) { m_nodes[id] = node; }
  inline void add_edge(NodeID from, NodeID to, const EdgeType &edge_data) {
    m_edges[{from, to}] = edge_data;
  }
  inline const Nodes &get_nodes() const { return m_nodes; }
  inline const Edges &get_edges() const { return m_edges; }

  // inline std::vector<NodeID> get_node_ids() const {
  //   std::vector<NodeID> ids;
  //   for (const auto &[id, _] : m_nodes) {
  //     ids.push_back(id);
  //   }
  //   return ids;
  // }
  //
  // inline std::vector<NodeType> get_node_types() const {
  //   std::vector<NodeType> types;
  //   for (const auto &[_, type] : m_nodes) {
  //     types.push_back(type);
  //   }
  //   return types;
  // }

private:
  Nodes m_nodes;
  Edges m_edges;
};

}; // namespace trajan::core::graph
