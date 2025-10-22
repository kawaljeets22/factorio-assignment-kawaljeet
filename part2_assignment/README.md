# Factorio Assignment - Design Notes

This repository contains C++ solutions for the "Factory Steady State" and "Belts with Bounds" problems.

The tools are built with CMake and are designed to read JSON from `stdin` and write a single JSON object to `stdout`.

## 1. Prerequisites

To build this project, you only need:
* **CMake** (version 3.21 or higher)
* A modern **C++17 compiler** (e.g., `g++` 9+ or `clang` 10+)
* **Git** (for downloading dependencies)

On Ubuntu, all prerequisites can be installed with:
`sudo apt update && sudo apt install build-essential cmake git`

All C++ library dependencies (`Google OR-Tools` and `nlohmann-json`) are handled automatically by CMake's `FetchContent` module. They will be downloaded and linked during the build process.

## 2. Build & Run

See **`RUN.md`** for detailed, step-by-step build and execution commands.

---

## 3. Factory Steady State (`factory`) Design

The `factory` problem is modeled as a two-phase **Linear Programming (LP)** problem, solved using Google's **OR-Tools** library.

### 3.1. Modeling Choices

* **Core Variables**: The primary variables are `x_r >= 0` for each recipe `r`, representing the number of "crafts per minute" for that recipe.
* **Item Conservation**: A conservation constraint is created for *every* item (raw, intermediate, and target).
    * `Balance(i) = [Total Production of i] - [Total Consumption of i]`
    * `Total Production(i) = Σ_r [x_r * out_r[i] * (1 + prod_r)]`
    * `Total Consumption(i) = Σ_r [x_r * in_r[i]]`
* **Constraint Bounds**:
    * **Intermediates**: `Balance(i) = 0`. The system is in steady-state, so all intermediates must be perfectly balanced.
    * **Raw Items**: `-raw_cap[i] <= Balance(i) <= 0`. The balance must be non-positive (net consumption) and the absolute value of the balance (the consumption) must not exceed the supply cap.
* **Machine Capacity**: A constraint is created for *each machine type* `m`.
    * `Total Machines(m) = Σ_{r uses m} [Machines Used by r] <= max_machines[m]`
    * `Machines Used by r = x_r / eff_crafts_per_min(r)`
    * The term `1.0 / eff_crafts_per_min(r)` is pre-calculated as `machine_cost_per_craft`.
* **Module Application**: Modules are applied *per-machine-type*. This logic is handled during data pre-processing.
* **Cycles & Byproducts**: This model handles cycles and byproducts automatically.
    * **Cycles (A→B→A)**: By enforcing `Balance(i) = 0` for all intermediates, the solver finds a valid flow `x_r` that satisfies the cycle.
    * **Byproducts**: By enforcing `Balance(i) = 0` for all intermediates, any byproduct `i` must either be fully consumed by another recipe or the recipe producing it will have its `x_r` set to 0.

### 3.2. Infeasibility and Tie-Breaking

A two-phase LP approach is used to meet the requirements.

* **Phase 1: Find Max Feasible Target Rate**
    1.  We introduce a new variable, `T >= 0`, representing the target item's production rate.
    2.  The target item's constraint becomes `Balance(target) - T = 0`.
    3.  The **objective is to MAXIMIZE `T`**.
    4.  We solve this LP. The result is `max_T`.

* **Phase 2: Minimize Machines**
    1.  We check if `max_T >= requested_target_rate - tolerance`.
    2.  **If No (Infeasible)**: We report `status: infeasible` with `max_feasible_target_per_min: max_T`. Bottlenecks are identified by checking the *dual values* (shadow prices) of the constraints from Phase 1.
    3.  **If Yes (Feasible)**: We formulate a *new* LP.
        * The target item constraint is now fixed: `Balance(target) = requested_target_rate`.
        * The **objective is to MINIMIZE `Σ_r [x_r * machine_cost_per_craft]`** (total machines).
        * We solve this second LP. The resulting `x_r` values are the final solution.

---

## 4. Belts with Bounds (`belts`) Design

The `belts` problem is modeled as a **Feasible Flow** problem (specifically, a *circulation with demands*) and solved using a **Max-Flow** algorithm (from **OR-Tools**) on a transformed graph.

### 4.1. Modeling Choices

1.  **Node Capacity Splitting**:
    * All nodes are given an `_in` index and an `_out` index.
    * For any uncapped node `v`, `v_in` and `v_out` point to the *same* index.
    * For a capped intermediate node `v` (not source/sink), `v_in` and `v_out` are *different* indices. An arc `(v_in -> v_out)` is added with capacity `cap(v)`.
    * All original edges `(u -> v)` are remapped to `(u_out -> v_in)`.

2.  **Lower Bounds Transformation (Circulation with Demands)**:
    * The problem is transformed into a standard max-flow problem on a *new* graph with a super-source `s*` and super-sink `t*`.
    * For each original edge `e = (u -> v)` with bounds `[lo, hi]`:
        * An arc is added to the max-flow graph `(u_out -> v_in)` with capacity `hi - lo`.
        * The lower bound `lo` is "pre-paid" by creating a *demand* of `lo` at `v` and a *supply* of `lo` at `u`.
    * A **balance** `B[i]` is tracked for each node `i`.
        * `B[i] = [Original Demand] - [Original Supply] + [Σ lo_in] - [Σ lo_out]`
        * For a source `s`, `B[s_out] -= supply(s)`.
        * For the sink `t`, `B[t_in] += total_supply`.

3.  **Feasibility Check**:
    * After calculating all balances, for every node `i`:
        * If `B[i] > 0` (net demand), add arc `(s* -> i)` with capacity `B_i`.
        * If `B[i] < 0` (net supply), add arc `(i -> t*)` with capacity `-B_i`.
    * We calculate `Total_Demand = Σ_{B[i] > 0} B[i]`.
    * We run max-flow from `s*` to `t*`.
    * The flow is **feasible if and only if `max_flow == Total_Demand`** (within `1e-9` tolerance).

### 4.2. Infeasibility Certificate

* If `max_flow < Total_Demand`, the flow is infeasible.
* The **min-cut** (which equals the max-flow) provides the certificate.
* We call `GetSourceSideMinCut()` to get all nodes reachable from `s*` in the residual graph. These are mapped back to their original string names and reported in `cut_reachable`.
* The `deficit` is `Total_Demand - max_flow`.

---

## 5. Numeric Approach & Determinism

* **Numeric Tolerance**: A standard tolerance of `1e-9` is used for all floating-point comparisons.
* **Solver**:
    * `factory`: **OR-Tools LP Solver (GLOP)**. We explicitly select the `DUAL_SIMPLEX` algorithm, which is deterministic and robust.
    * `belts`: **OR-Tools Max-Flow Solver**. This is a deterministic push-relabel implementation.
* **Tie-Breaking**:
    * `factory`: Handled by the `minimize total machines` objective in Phase 2. Any remaining ties are broken deterministically by the `DUAL_SIMPLEX` algorithm.
    * `belts`: The max-flow algorithm is deterministic.