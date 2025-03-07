#pragma once
#include <functional>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace trajan::core {

template <typename NodeType, typename EdgeType> class ConnectedComponent;

template <typename NodeType, typename EdgeType> class Graph {
public:
  using NodeId = int;
  using EdgePredicate = std::function<std::optional<EdgeType>(
      const NodeType &, const NodeType &)>;

  Graph(const std::vector<NodeType> &nodes, EdgePredicate edge_predicate)
      : m_nodes(nodes), m_edge_predicate(edge_predicate) {}

  void build_graph() {
    for (size_t i = 0; i < m_nodes.size(); i++) {
      NodeId id1 = get_node_id(i);

      for (size_t j = i + 1; j < m_nodes.size(); j++) {
        NodeId id2 = get_node_id(j);
        std::optional<EdgeType> edge_data =
            m_edge_predicate(m_nodes[i], m_nodes[j]);
        if (!edge_data) {
          continue;
        }
        m_adjacency_list[id1][id2] = *edge_data;
        m_adjacency_list[id2][id1] = *edge_data;
      }
    }
    m_init = true;
  }

  std::vector<ConnectedComponent<NodeType, EdgeType>>
  find_connected_components() const {
    if (!m_init) {
      throw std::runtime_error("Graph not built. Call build_graph().");
    }

    std::vector<ConnectedComponent<NodeType, EdgeType>> components;
    std::unordered_set<NodeId> visited;

    for (size_t i = 0; i < m_nodes.size(); i++) {
      NodeId node_id = get_node_id(i);

      if (visited.find(node_id) != visited.end()) {
        continue;
      }

      ConnectedComponent<NodeType, EdgeType> component;

      std::queue<NodeId> queue;
      queue.push(node_id);
      visited.insert(node_id);

      while (!queue.empty()) {
        NodeId current_id = queue.front();
        queue.pop();
        component.add_node(current_id, get_node_by_id(current_id));

        if (m_adjacency_list.find(current_id) != m_adjacency_list.end()) {
          for (const auto &[neighbor_id, edge_data] :
               m_adjacency_list.at(current_id)) {
            component.add_edge(current_id, neighbor_id, edge_data);

            if (visited.find(neighbor_id) == visited.end()) {
              queue.push(neighbor_id);
              visited.insert(neighbor_id);
            }
          }
        }
      }

      components.push_back(component);
    }

    return components;
  }

  const NodeType &get_node_by_id(NodeId id) const {
    for (const auto &node : m_nodes) {
      if (get_node_id_from_node(node) == id) {
        return node;
      }
    }
    throw std::runtime_error("Node ID not found");
  }

  const std::unordered_map<NodeId, std::unordered_map<NodeId, EdgeType>> &
  get_adjacency_list() const {
    return m_adjacency_list;
  }

  virtual NodeId get_node_id_from_node(const NodeType &node) const = 0;

protected:
  std::vector<NodeType> m_nodes;
  EdgePredicate m_edge_predicate;
  std::unordered_map<NodeId, std::unordered_map<NodeId, EdgeType>>
      m_adjacency_list;
  bool m_init = false;

  NodeId get_node_id(size_t index) const {
    return get_node_id_from_node(m_nodes[index]);
  }
};

template <typename NodeType, typename EdgeType> class ConnectedComponent {
public:
  using NodeId = typename Graph<NodeType, EdgeType>::NodeId;

  void add_node(NodeId id, const NodeType &node) { m_nodes[id] = node; }

  void add_edge(NodeId from, NodeId to, const EdgeType &edge_data) {
    m_edges[{from, to}] = edge_data;
  }

  const std::unordered_map<NodeId, NodeType> &get_nodes() const {
    return m_nodes;
  }

  struct PairHash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const {
      auto h1 = std::hash<T1>{}(p.first);
      auto h2 = std::hash<T2>{}(p.second);
      return h1 ^ (h2 << 1);
    }
  };

  const std::unordered_map<std::pair<NodeId, NodeId>, EdgeType, PairHash> &
  get_edges() const {
    return m_edges;
  }

  std::vector<NodeId> get_node_ids() const {
    std::vector<NodeId> ids;
    for (const auto &[id, _] : m_nodes) {
      ids.push_back(id);
    }
    return ids;
  }

  std::vector<NodeType> get_node_types() const {
    std::vector<NodeType> types;
    for (const auto &[_, type] : m_nodes) {
      types.push_back(type);
    }
    return types;
  }

private:
  std::unordered_map<NodeId, NodeType> m_nodes;
  std::unordered_map<std::pair<NodeId, NodeId>, EdgeType, PairHash> m_edges;
};

}; // namespace trajan::core
