#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>
#include <iomanip>
#include "nlohmann/json.hpp"
#include "ortools/linear_solver/linear_solver.h"

// Use nlohmann::json for convenience
using json = nlohmann::json;
// Use Google OR-Tools
namespace operations_research {

// Helper structure to hold recipe details
struct Recipe {
    std::string name;
    std::string machine;
    double time_s;
    std::map<std::string, double> in;
    std::map<std::string, double> out;
    
    // Calculated effective values
    double eff_crafts_per_min = 0.0;
    double machine_cost_per_craft = 0.0; // 1.0 / eff_crafts_per_min
    double prod_mult = 1.0;
};

// Main solver function for the factory problem
void solve_factory() {
    json input;
    try {
        // Read all stdin into the json object
        std::cin >> input;
    } catch (json::parse_error& e) {
        std::cerr << "JSON parse error: " << e.what() << std::endl;
        return;
    }

    // --- 1. Data Parsing and Pre-computation ---

    std::map<std::string, Recipe> recipes;
    std::set<std::string> all_items;
    std::set<std::string> raw_items;
    std::set<std::string> intermediate_items;
    std::map<std::string, std::string> recipe_machine_map;
    std::map<std::string, double> machine_caps;
    std::map<std::string, double> raw_caps;
    std::string target_item = input["target"]["item"];
    double requested_target_rate = input["target"]["rate_per_min"];

    // Initialize raw items from supply limits
    for (auto const& [item, cap] : input["limits"]["raw_supply_per_min"].items()) {
        raw_items.insert(item);
        all_items.insert(item);
        raw_caps[item] = cap;
    }

    // Initialize machine caps
    for (auto const& [machine, cap] : input["limits"]["max_machines"].items()) {
        machine_caps[machine] = cap;
    }

    // Process recipes
    for (auto const& [name, data] : input["recipes"].items()) {
        Recipe r;
        r.name = name;
        r.machine = data["machine"];
        r.time_s = data["time_s"];
        recipe_machine_map[name] = r.machine;

        if (data.contains("in")) {
            for (auto const& [item, qty] : data["in"].items()) {
                r.in[item] = qty;
                all_items.insert(item);
            }
        }
        if (data.contains("out")) {
            for (auto const& [item, qty] : data["out"].items()) {
                r.out[item] = qty;
                all_items.insert(item);
            }
        }

        // Get machine and module data
        double base_speed = input["machines"][r.machine]["crafts_per_min"];
        double mod_speed = 0.0;
        double mod_prod = 0.0;
        if (input["modules"].contains(r.machine)) {
            if (input["modules"][r.machine].contains("speed")) {
                mod_speed = input["modules"][r.machine]["speed"];
            }
            if (input["modules"][r.machine].contains("prod")) {
                mod_prod = input["modules"][r.machine]["prod"];
            }
        }

        r.prod_mult = 1.0 + mod_prod;
        r.eff_crafts_per_min = base_speed * (1.0 + mod_speed) * 60.0 / r.time_s;
        
        // Avoid division by zero if time_s or eff_crafts is 0
        if (r.eff_crafts_per_min > 1e-9) {
            r.machine_cost_per_craft = 1.0 / r.eff_crafts_per_min;
        } else {
            // This recipe is impossible to craft, give it an infinite cost
            r.machine_cost_per_craft = 1e30; 
        }

        recipes[name] = r;
    }

    // Classify items
    all_items.erase(target_item);
    for (const auto& item : all_items) {
        if (raw_items.find(item) == raw_items.end()) {
            intermediate_items.insert(item);
        }
    }
    
    // --- 2. Phase 1: Find Max Feasible Target Rate ---
    // We create an LP to maximize the target rate, subject to all constraints.

    std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("GLOP")); // GLOP is good for LPs. Or use "DUAL_SIMPLEX"
    if (!solver) {
        // Fallback or error
        solver = MPSolver::CreateSolver("GLOP");
    }
    // Use DUAL_SIMPLEX for determinism and robustness
    MPSolverParameters params;
    params.SetIntegerParam(MPSolverParameters::LP_ALGORITHM, MPSolverParameters::DUAL_SIMPLEX);
    params.SetDoubleParam(MPSolverParameters::RANDOM_SEED, 1.0); // Deterministic seed

    const double infinity = solver->infinity();
    
    // Variables: x_r (crafts/min for recipe r)
    std::map<std::string, MPVariable*> x;
    for (auto const& [name, r] : recipes) {
        x[name] = solver->MakeNumVar(0.0, infinity, name);
    }
    // Variable: T (target rate)
    MPVariable* T = solver->MakeNumVar(0.0, infinity, "target_rate");

    // Constraints:
    // a) Item Conservation
    std::map<std::string, MPConstraint*> item_balance;
    for (const auto& item : all_items) {
        item_balance[item] = solver->MakeRowConstraint(0.0, 0.0, "balance_" + item);
        if (raw_items.count(item)) {
            // Raw items: consumption <= cap -> balance >= -cap
            item_balance[item]->SetBounds(-raw_caps[item], 0.0);
        }
    }
    // Target item balance: balance = T
    item_balance[target_item] = solver->MakeRowConstraint(0.0, 0.0, "balance_" + target_item);
    item_balance[target_item]->SetCoefficient(T, -1.0);

    // b) Machine Usage
    std::map<std::string, MPConstraint*> machine_usage;
    for (auto const& [name, cap] : machine_caps) {
        machine_usage[name] = solver->MakeRowConstraint(0.0, cap, "cap_" + name);
    }

    // Populate constraints
    for (auto const& [name, r] : recipes) {
        // Add to machine usage
        if (machine_usage.count(r.machine)) {
            machine_usage[r.machine]->SetCoefficient(x[name], r.machine_cost_per_craft);
        }

        // Add to item balances
        for (auto const& [item, qty] : r.in) {
            if (item_balance.count(item)) {
                item_balance[item]->SetCoefficient(x[name], -qty);
            }
        }
        for (auto const& [item, qty] : r.out) {
            if (item_balance.count(item)) {
                item_balance[item]->SetCoefficient(x[name], qty * r.prod_mult);
            }
        }
    }

    // Objective: Maximize T
    MPObjective* const objective1 = solver->MutableObjective();
    objective1->SetCoefficient(T, 1.0);
    objective1->SetMaximization();

    // Solve Phase 1
    const MPSolver::ResultStatus status1 = solver->Solve(params);

    if (status1 != MPSolver::OPTIMAL && status1 != MPSolver::FEASIBLE) {
        // This should not happen if 0 is a feasible solution
        json result = {
            {"status", "infeasible"},
            {"max_feasible_target_per_min", 0.0},
            {"bottleneck_hint", {"Initial solver failure"}}
        };
        std::cout << result.dump(2) << std::endl;
        return;
    }

    double max_feasible_rate = T->solution_value();
    const double tolerance = 1e-9;

    // --- 3. Check Feasibility and Handle Infeasible Case ---
    if (max_feasible_rate < requested_target_rate - tolerance) {
        json result;
        result["status"] = "infeasible";
        result["max_feasible_target_per_min"] = max_feasible_rate;
        
        std::set<std::string> hints;
        // Check binding constraints (bottlenecks)
        for (auto const& [name, constraint] : machine_usage) {
            if (constraint->dual_value() > tolerance) {
                hints.insert(name + " cap");
            }
        }
        for (auto const& [item, constraint] : item_balance) {
            if (raw_items.count(item) && constraint->dual_value() > tolerance) {
                hints.insert(item + " supply");
            }
        }
        if (hints.empty() && max_feasible_rate > tolerance) {
             hints.insert("Target rate conflicts with other constraints");
        } else if (hints.empty()) {
             hints.insert("Unknown bottleneck, possibly no production path");
        }

        result["bottleneck_hint"] = hints;
        std::cout << result.dump(2) << std::endl;
        return;
    }

    // --- 4. Phase 2: Minimize Machines for Feasible Target Rate ---
    // The requested rate is feasible. Now we solve a new LP to minimize machines.
    
    std::unique_ptr<MPSolver> solver2(MPSolver::CreateSolver("GLOP"));
    if (!solver2) {
        solver2 = MPSolver::CreateSolver("GLOP");
    }

    // Variables: x_r (crafts/min for recipe r)
    std::map<std::string, MPVariable*> x2;
    for (auto const& [name, r] : recipes) {
        x2[name] = solver2->MakeNumVar(0.0, infinity, name);
    }

    // Constraints:
    // a) Item Conservation
    std::map<std::string, MPConstraint*> item_balance2;
    for (const auto& item : intermediate_items) {
        item_balance2[item] = solver2->MakeRowConstraint(0.0, 0.0, "balance_" + item);
    }
    for (const auto& item : raw_items) {
        item_balance2[item] = solver2->MakeRowConstraint(-raw_caps[item], 0.0, "balance_" + item);
    }
    // Target item: Must be *exactly* the requested rate
    item_balance2[target_item] = solver2->MakeRowConstraint(requested_target_rate, requested_target_rate, "balance_" + target_item);

    // b) Machine Usage
    std::map<std::string, MPConstraint*> machine_usage2;
    for (auto const& [name, cap] : machine_caps) {
        machine_usage2[name] = solver2->MakeRowConstraint(0.0, cap, "cap_" + name);
    }
    
    // Objective: Minimize total machines
    MPObjective* const objective2 = solver2->MutableObjective();
    objective2->SetMinimization();

    // Populate constraints and objective
    for (auto const& [name, r] : recipes) {
        // Add to machine usage constraint
        if (machine_usage2.count(r.machine)) {
            machine_usage2[r.machine]->SetCoefficient(x2[name], r.machine_cost_per_craft);
        }
        // Add to objective function
        objective2->SetCoefficient(x2[name], r.machine_cost_per_craft);

        // Add to item balances
        for (auto const& [item, qty] : r.in) {
            if (item_balance2.count(item)) {
                item_balance2[item]->SetCoefficient(x2[name], -qty);
            }
        }
        for (auto const& [item, qty] : r.out) {
            if (item_balance2.count(item)) {
                item_balance2[item]->SetCoefficient(x2[name], qty * r.prod_mult);
            }
        }
    }

    // Solve Phase 2
    const MPSolver::ResultStatus status2 = solver2->Solve(params);

    if (status2 != MPSolver::OPTIMAL && status2 != MPSolver::FEASIBLE) {
        // This should not happen if Phase 1 succeeded
        json result = {
            {"status", "infeasible"},
            {"max_feasible_target_per_min", max_feasible_rate},
            {"bottleneck_hint", {"Phase 2 solver failure"}}
        };
        std::cout << result.dump(2) << std::endl;
        return;
    }

    // --- 5. Format Success Output ---
    json result;
    result["status"] = "ok";
    
    std::map<std::string, double> per_recipe_crafts;
    std::map<std::string, double> per_machine_counts;
    std::map<std::string, double> raw_consumption;

    // Initialize all machine counts to 0
    for (auto const& [name, cap] : machine_caps) {
        per_machine_counts[name] = 0.0;
    }
    // Initialize all recipe crafts to 0
    for (auto const& [name, r] : recipes) {
        per_recipe_crafts[name] = 0.0;
    }

    for (auto const& [name, r] : recipes) {
        double crafts = x2[name]->solution_value();
        if (crafts > tolerance) {
            per_recipe_crafts[name] = crafts;
            per_machine_counts[r.machine] += crafts * r.machine_cost_per_craft;
        }
    }
    
    for (const auto& item : raw_items) {
        // Activity is the value of the constraint's linear expression
        double balance = item_balance2[item]->activity();
        // Consumption is the negative of the balance
        double consumption = -balance;
        if (consumption > tolerance) {
            raw_consumption[item] = consumption;
        } else {
             // Explicitly add 0 if it was in the caps list
            if (raw_caps.count(item)) {
                raw_consumption[item] = 0.0;
            }
        }
    }

    result["per_recipe_crafts_per_min"] = per_recipe_crafts;
    result["per_machine_counts"] = per_machine_counts;
    result["raw_consumption_per_min"] = raw_consumption;

    std::cout << std::setprecision(10) << result.dump(2) << std::endl;
}

} // namespace operations_research

int main() {
    try {
        operations_research::solve_factory();
    } catch (const std::exception& e) {
        std::cerr << "Unhandled exception: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}