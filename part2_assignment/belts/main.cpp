#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>
#include <iomanip>
#include "nlohmann/json.hpp"
#include "ortools/graph/max_flow.h"

// Use nlohmann::json for convenience
using json = nlohmann::json;
// Use Google OR-Tools
namespace operations_research {

// Helper structure for original edges
struct Edge {
    std::string from;
    std::string to;
    double lower_bound;
    double upper_bound;
    int arc_index = -1; // Index in the OR-Tools graph
};

// Main solver function for the belts problem
void solve_belts() {
    json input;
    try {
        // Read all stdin into the json object
        std::cin >> input;
    } catch (json::parse_error& e) {
        std::cerr << "JSON parse error: " << e.what() << std::endl;
        return;
    }

    // --- 1. Data Parsing and Graph Transformation Setup ---

    SimpleMaxFlow max_flow_solver;
    
    std::map<std::string, double> sources = input["sources"];
    std::string sink_node = input["sink"];
    std::map<std::string, double> node_caps;
    if (input.contains("node_caps")) {
        node_caps = input["node_caps"].get<std::map<std::string, double>>();
    }
    std::vector<Edge> original_edges;
    for (const auto& edge_data : input["edges"]) {
        original_edges.push_back({
            edge_data["from"],
            edge_data["to"],
            edge_data["lower_bound"],
            edge_data["upper_bound"]
        });
    }

    std::set<std::string> all_node_names;
    for (const auto& [name, supply] : sources) {
        all_node_names.insert(name);
    }
    all_node_names.insert(sink_node);
    for (const auto& edge : original_edges) {
        all_node_names.insert(edge.from);
        all_node_names.insert(edge.to);
    }
    for (const auto& [name, cap] : node_caps) {
        all_node_names.insert(name);
    }

    // Node mapping:
    // We map string names to integer indices.
    // For capped nodes, we split: v -> v_in, v_out
    int node_index_counter = 0;
    std::map<std::string, int> node_to_in_index;
    std::map<std::string, int> node_to_out_index;
    std::vector<std::string> index_to_name_map;

    auto get_name_from_index = [&](int index) {
        if (index >= 0 && index < index_to_name_map.size()) {
            return index_to_name_map[index];
        }
        return std::string("INVALID_INDEX");
    };

    // We need a super-source (s_star) and super-sink (t_star)
    // for the feasibility circulation problem.
    const int s_star = node_index_counter++;
    index_to_name_map.push_back("s_star");
    const int t_star = node_index_counter++;
    index_to_name_map.push_back("t_star");

    for (const auto& name : all_node_names) {
        int in_idx = node_index_counter++;
        index_to_name_map.push_back(name); // Map both indices back to the same name
        node_to_in_index[name] = in_idx;

        if (node_caps.count(name) && !sources.count(name) && name != sink_node) {
            // This is a capped intermediate node. Split it.
            int out_idx = node_index_counter++;
            index_to_name_map.push_back(name);
            node_to_out_index[name] = out_idx;
            // Add the capacity edge from v_in -> v_out
            max_flow_solver.AddArcWithCapacity(in_idx, out_idx, node_caps[name]);
        } else {
            // Not capped, or is a source/sink. In and Out are the same.
            node_to_out_index[name] = in_idx;
        }
    }

    // --- 2. Build the Feasibility Graph (Circulation with Demands) ---

    // B[i] = net demand at node i.
    // B[i] > 0 -> demand, needs flow from s_star
    // B[i] < 0 -> supply, sends flow to t_star
    std::map<int, double> node_balance;
    double total_supply = 0.0;
    double total_demand_at_sink = 0.0;
    const double tolerance = 1e-9;

    // a) Add imbalances from supplies and sink demand
    for (const auto& [name, supply] : sources) {
        node_balance[node_to_out_index[name]] -= supply;
        total_supply += supply;
    }
    // Sink demand is total supply
    node_balance[node_to_in_index[sink_node]] += total_supply;
    total_demand_at_sink = total_supply;

    // b) Add imbalances from edge lower bounds
    for (auto& edge : original_edges) {
        int u = node_to_out_index[edge.from];
        int v = node_to_in_index[edge.to];
        double lo = edge.lower_bound;
        double hi = edge.upper_bound;

        if (hi - lo < -tolerance) {
            // Infeasible: upper < lower
            json result = {
                {"status", "infeasible"},
                {"cut_reachable", {}},
                {"deficit", {
                    {"demand_balance", total_demand_at_sink},
                    {"tight_nodes", {}},
                    {"tight_edges", {
                        {"from", edge.from},
                        {"to", edge.to},
                        {"flow_needed", hi - lo}
                    }}
                }}
            };
            std::cout << result.dump(2) << std::endl;
            return;
        }
        
        // Add residual capacity arc
        edge.arc_index = max_flow_solver.AddArcWithCapacity(u, v, hi - lo);

        // Adjust node balances
        node_balance[u] -= lo;
        node_balance[v] += lo;
    }

    // c) Add s_star and t_star edges to balance the demands
    double total_demand_from_s_star = 0.0;
    for (const auto& [node_idx, balance] : node_balance) {
        if (balance > tolerance) {
            // Net demand at node_idx
            max_flow_solver.AddArcWithCapacity(s_star, node_idx, balance);
            total_demand_from_s_star += balance;
        } else if (balance < -tolerance) {
            // Net supply at node_idx
            max_flow_solver.AddArcWithCapacity(node_idx, t_star, -balance);
        }
    }

    // --- 3. Solve Max-Flow and Check Feasibility ---

    if (max_flow_solver.Solve(s_star, t_star) != SimpleMaxFlow::OPTIMAL) {
        json result = {
            {"status", "infeasible"},
            {"cut_reachable", {}},
            {"deficit", {
                {"demand_balance", total_demand_from_s_star},
                {"tight_nodes", {}},
                {"tight_edges", {}}
            }}
        };
        std::cout << result.dump(2) << std::endl;
        return;
    }

    double max_flow = max_flow_solver.OptimalFlow();

    if (max_flow < total_demand_from_s_star - tolerance) {
        // Infeasible: Could not satisfy all demands
        json result;
        result["status"] = "infeasible";
        result["deficit"] = {
            {"demand_balance", total_demand_from_s_star - max_flow}
        };

        std::vector<int> cut_indices;
        max_flow_solver.GetSourceSideMinCut(&cut_indices);

        std::set<std::string> reachable_nodes;
        for (int idx : cut_indices) {
            if (idx != s_star) {
                reachable_nodes.insert(get_name_from_index(idx));
            }
        }
        result["cut_reachable"] = reachable_nodes;
        
        // Find tight constraints (this is illustrative)
        std::set<std::string> tight_nodes_set;
        json tight_edges_json = json::array();

        for (int u = 0; u < max_flow_solver.NumNodes(); ++u) {
            if (reachable_nodes.count(get_name_from_index(u))) {
                for (int arc_idx = 0; arc_idx < max_flow_solver.NumArcs(); ++arc_idx) {
                    if (max_flow_solver.Tail(arc_idx) == u) {
                        int v = max_flow_solver.Head(arc_idx);
                        if (!reachable_nodes.count(get_name_from_index(v))) {
                            // This arc crosses the cut
                            // Check if it's an original edge or a node cap
                            bool is_original = false;
                            for (const auto& edge : original_edges) {
                                if (edge.arc_index == arc_idx) {
                                    tight_edges_json.push_back({
                                        {"from", edge.from},
                                        {"to", edge.to},
                                        {"flow_needed", "at capacity"}
                                    });
                                    is_original = true;
                                    break;
                                }
                            }
                            if (!is_original) {
                                // Must be a node cap
                                std::string node_name = get_name_from_index(u);
                                if (node_caps.count(node_name)) {
                                    tight_nodes_set.insert(node_name);
                                }
                            }
                        }
                    }
                }
            }
        }
        result["deficit"]["tight_nodes"] = tight_nodes_set;
        result["deficit"]["tight_edges"] = tight_edges_json;

        std::cout << result.dump(2) << std::endl;
        return;
    }

    // --- 4. Format Success Output ---
    json result;
    result["status"] = "ok";
    result["max_flow_per_min"] = total_supply; // Total flow *into* the sink
    
    json flows = json::array();
    for (const auto& edge : original_edges) {
        // f_original = f_residual + lower_bound
        double residual_flow = max_flow_solver.Flow(edge.arc_index);
        double final_flow = residual_flow + edge.lower_bound;
        
        if (final_flow > tolerance) {
            flows.push_back({
                {"from", edge.from},
                {"to", edge.to},
                {"flow", final_flow}
            });
        }
    }
    result["flows"] = flows;

    std::cout << std::setprecision(10) << result.dump(2) << std::endl;
}

} // namespace operations_research

int main() {
    try {
        operations_research::solve_belts();
    } catch (const std::exception& e) {
        std::cerr << "Unhandled exception: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}