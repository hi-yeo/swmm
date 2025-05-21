# Python SWMM Algorithms: Testing and Debugging Strategy

This document outlines the strategy for testing and debugging the Python functions translated from specific SWMM (Stormwater Management Model) C code algorithms. These Python functions are intended to replicate the core computational logic of their C counterparts for educational and analytical purposes.

## 1. Overall Goal of Testing

The primary goal of this testing effort is to **verify that the Python functions correctly implement the algorithmic logic and mathematical formulas** of their corresponding SWMM C code counterparts. The key Python functions under scrutiny include:

*   `calculate_conduit_trial_flow_momentum_eq`
*   Preissmann Slot geometry functions:
    *   `get_preissmann_slot_width_py`
    *   `get_effective_top_width_py`
    *   `get_effective_area_py`
    *   `get_effective_hydraulic_radius_py`
*   `calculate_conduit_trial_flow_with_slot_py`
*   `perform_picard_iteration_for_link_system`
*   Δt (time step) adjustment functions:
    *   `get_link_step_py`
    *   `get_node_step_py`
    *   `get_swmm_variable_routing_step_py`

It is crucial to emphasize that this testing is **not about replicating full SWMM simulation model outputs**. Instead, the focus is on validating these specific, isolated algorithms as faithful Python translations of the C code's computational steps and formulas.

## 2. General Testing Principles

The following principles will guide the testing process for each Python function or logical group of functions:

*   **Test Cases:** For each function, a diverse set of test cases will be defined. These will include:
    *   **Valid Inputs:** Typical values that the algorithm is expected to handle during a simulation.
    *   **Edge Case Inputs:** Values at the boundaries of valid ranges, zero values, very small (FUDGE-level) values, and values near thresholds that trigger different logical paths (e.g., around `y_full_conduit` for Preissmann slot functions).
    *   **Typical Inputs:** Representative values based on common hydraulic scenarios.
*   **Expected Outputs:** For each defined test case, the expected output will be determined through one or more of the following methods:
    *   **Manual Calculation:** For functions involving direct mathematical formulas (e.g., individual `dq` terms in the momentum equation, Preissmann slot geometry), outputs will be calculated manually using a calculator, spreadsheet, or symbolic math tools.
    *   **Comparison with SWMM C Code Values:** If feasible for very simple, isolated scenarios, intermediate values from the SWMM C code (obtained via debugging or temporary print statements in a controlled C environment) might serve as a reference. This is complex and will be used sparingly.
    *   **Logical Reasoning:** Based on the algorithm's purpose and the input values, the expected behavior or range of outputs will be determined (e.g., a very high friction factor should significantly reduce flow).
*   **Korean Comments in Test Code:** Any example test invocation code, particularly if included in the `if __name__ == '__main__':` blocks of the Python files, should ideally also feature Korean comments to maintain consistency with the function documentation.

## 3. Specific Strategy for Each Algorithm Group

### A. Conduit Flow Momentum Equation & Preissmann Slot Geometry

**Functions:**
*   `calculate_conduit_trial_flow_momentum_eq`
*   `get_preissmann_slot_width_py`
*   `get_effective_top_width_py`
*   `get_effective_area_py`
*   `get_effective_hydraulic_radius_py`
*   `calculate_conduit_trial_flow_with_slot_py`

**Testing `get_preissmann_slot_width_py`:**
*   **Inputs:**
    *   Vary `y_depth` relative to `y_full_conduit`:
        *   `y_depth < slot_crown_cutoff * y_full_conduit`
        *   `slot_crown_cutoff * y_full_conduit <= y_depth <= PREISSMANN_SLOT_HIGHER_Y_NORM_LIMIT * y_full_conduit` (Sjoberg range)
        *   `y_depth > PREISSMANN_SLOT_HIGHER_Y_NORM_LIMIT * y_full_conduit`
    *   Test `surcharge_method_option` with "SLOT" and "EXTRAN".
    *   Test `is_conduit_closed` with `True` and `False`.
*   **Expected:** Correct `w_slot` calculated: 0 for non-SLOT/open conduit/sub-cutoff depth; Sjoberg's formula value within its range; or the minimum slot width factor times `w_max_conduit` above the Sjoberg upper limit.

**Testing `get_effective_top_width_py`, `get_effective_area_py`, `get_effective_hydraulic_radius_py`:**
*   **Inputs:**
    *   Vary `y_depth` across sub-full, just-full, and surcharged conditions.
    *   Use `w_slot` values obtained from testing `get_preissmann_slot_width_py` (0 and positive values).
    *   For `physical_..._at_y_func` callbacks, use simplified functions for a rectangular conduit (e.g., `lambda y, y_full, w_max: w_max`).
*   **Expected:**
    *   `get_effective_top_width_py`: Should return `w_slot` if `w_slot > 0`, otherwise the physical top width (potentially adjusted near crown).
    *   `get_effective_area_py`: Should return `a_full_conduit + (y_depth - y_full_conduit) * w_slot` if `w_slot > 0` and `y_depth >= y_full_conduit`, otherwise the physical area.
    *   `get_effective_hydraulic_radius_py`: Should return `r_full_conduit` if `w_slot > 0` and `y_depth >= y_full_conduit`, otherwise the physical hydraulic radius.

**Testing `calculate_conduit_trial_flow_momentum_eq` (and by extension `calculate_conduit_trial_flow_with_slot_py` which uses it after geometry adjustment):**
*   **Inputs:**
    *   Vary hydraulic conditions: `q_old_barrel`, `h1_node`, `h2_node`.
    *   Vary geometric inputs (these would be the *effective* geometries if testing via `calculate_conduit_trial_flow_with_slot_py`): `a1_current_iter`, `a2_current_iter`, `a_mid_current_iter`, `r1_current_iter`, `r_mid_current_iter`.
    *   Test specific scenarios:
        *   Initial condition: `q_old_barrel = 0`, `q_last_iter_barrel = 0`.
        *   Flow primarily driven by head difference (`h1_node > h2_node`).
        *   Flow opposed by adverse head difference (`h2_node > h1_node`).
        *   Dominant friction: High `conduit_roughness_factor`, significant `v_current_iter_barrel`.
        *   Dominant local losses: `has_local_losses = True` with non-zero coefficients.
        *   Dominant inertial terms: `sigma_inertial_damping > 0` with `a_mid_current_iter != a_old_barrel` (for `dq3`) or `a1_current_iter != a2_current_iter` (for `dq4`).
        *   Test behavior when `rWtd` is very small (near `FUDGE`), expecting `dq1` to become very large if velocity is present.
        *   Test behavior when `length_modified` is very small.
*   **Expected:** Manually calculate each `dq` term (`dq1` to `dq5`, with `dq6=0`) and then the final `q_trial_momentum_barrel` using the formula. Compare the Python function's output with the manual calculation.
*   **Focus:** Verify the correct calculation of individual formula components (`dq` terms) and the final combined result. For `calculate_conduit_trial_flow_with_slot_py`, ensure the Preissmann Slot logic correctly modifies the geometric inputs fed into the momentum equation solver.

### B. Picard Iteration (`perform_picard_iteration_for_link_system`)

*   **Inputs:**
    *   Define a simple system: one conduit, two nodes. Provide fixed physical conduit properties (width, `y_full`, `a_full`, `r_full`, `z1_conduit`, `z2_conduit`) and SWMM-specific parameters (`length_modified`, `length_actual`, `conduit_roughness_factor`, `num_barrels`). Provide fixed node surface areas.
    *   Vary initial conditions:
        *   `h_node1_initial`, `h_node2_initial` (creating different initial gradients).
        *   `q_old_conduit_total` (previous time step's flow).
        *   `node1_inflow_external`, `node2_inflow_external`.
    *   Test different iteration control parameters: `dt`, `omega_relaxation`, `max_trials`, `head_tolerance`.
    *   Use fixed `sigma_inertial_damping_fixed` and `rho_factor_fixed` for simplicity in these tests, as their dynamic calculation is outside this function's scope.
*   **Expected Outputs/Behavior:**
    *   **Convergence:** Verify if the `converged` flag is `True` when changes in `h_node1_iter` and `h_node2_iter` between iterations fall below `head_tolerance` within `max_trials`.
    *   **Iteration Count:** Check if the returned `steps` is reasonable for the scenario. For simple, stable cases, it should be low.
    *   **Final Values:** Assess if `final_q_conduit_total`, `final_h_node1`, `final_h_node2` are physically plausible and stable for the given inputs.
        *   Example Case 1 (Easy Convergence): Steady external inflows, mild positive slope, reasonable `dt` and `omega_relaxation`. Expect quick convergence.
        *   Example Case 2 (Challenging): Rapidly changing or oscillating external inflows (if simulated over multiple calls, though this function is for one `dt`), very small `dt` (may take more iterations if `head_tolerance` implies very small absolute head changes), or a very strict `head_tolerance`.
        *   Example Case 3 (Non-Convergence): Set `max_trials` very low for a case that needs more iterations. Expect `converged = False`.
*   **Focus:** Correct implementation of the iteration loop, proper application of under-relaxation for both flow (within the `calculate_conduit_trial_flow_momentum_eq` call or directly after) and (simplified) node heads, iterative refinement of flow and head values, and accurate convergence checking against `head_tolerance`. The simplified geometry calculations within the loop should be consistent.

### C. Δt Adjustment Logic (`get_link_step_py`, `get_node_step_py`, `get_swmm_variable_routing_step_py`)

*   **Inputs for `get_link_step_py`:**
    *   Create `links_data` (list of dictionaries) with diverse conduit states:
        *   Varying `q_barrel`, `volume_new_barrel` (implies area), `froude_conduit`.
        *   Varying `length_mod_conduit` and `length_actual_conduit`.
        *   Include cases that should be skipped (e.g., `q_barrel` near FUDGE, `froude_conduit` very low).
*   **Expected for `get_link_step_py`:**
    *   For a few representative links in `links_data`, manually calculate the Courant time step `t` using the formula:
      `t = (volume_new_barrel / q_barrel) * (length_mod_conduit / length_actual_conduit) * (froude_conduit / (1.0 + froude_conduit)) * courant_factor`
    *   Verify that the Python function calculates these individual `t` values correctly.
    *   Verify that the function correctly returns the minimum `t` among all eligible links and the ID of that critical link.
*   **Inputs for `get_node_step_py`:**
    *   Create `nodes_data` (list of dictionaries) with diverse node states:
        *   Varying `depth_new`, `depth_rate_of_change` (`dYdT`).
        *   Varying `crown_elevation` and `invert_elevation`.
        *   Include cases that should be skipped (outfalls, very low depth, very low `dYdT`, depth near crown).
*   **Expected for `get_node_step_py`:**
    *   For a few representative nodes, manually calculate `allowable_depth_change = (crown_elevation - invert_elevation) * 0.25` and then `t1_node = allowable_depth_change / abs(depth_rate_of_change)`.
    *   Verify the Python function's individual `t1_node` calculations and its selection of the minimum.
*   **Inputs for `get_swmm_variable_routing_step_py`:**
    *   Use the `links_data` and `nodes_data` created above.
    *   Test various combinations of:
        *   `courant_factor_option` (0.0 and a typical value like 0.75).
        *   `min_route_step_option` (e.g., 0.1s, 0.5s).
        *   `fixed_routing_step` (e.g., 1s, 5s).
        *   `current_variable_step_global` (0.0 for initial step, and a non-zero previous step).
*   **Expected for `get_swmm_variable_routing_step_py`:**
    *   **Variable step disabled:** If `courant_factor_option = 0.0`, output should be `fixed_routing_step`.
    *   **Start of simulation:** If `current_variable_step_global = 0.0`, output should be `max(min_route_step_option, MINTIMESTEP_CONST_PY)`.
    *   **Link step governs:** Craft inputs where the calculated link step is the smallest positive controlling factor.
    *   **Node step governs:** Craft inputs where the calculated node step is the smallest positive controlling factor.
    *   **`min_route_step_option` governs:** Craft inputs where both link and node calculated steps are smaller than `min_route_step_option`.
    *   **`fixed_routing_step` governs:** Craft inputs where link and node steps are larger than `fixed_routing_step`.
    *   Verify the final flooring of the result to the nearest millisecond.
*   **Focus:** Correct implementation of the specific formulas for link and node stability criteria, and accurate application of the selection logic (min of link/node steps) and bounding logic (`min_route_step_option`, `fixed_routing_step`) in the main `get_swmm_variable_routing_step_py` function.

## 4. Debugging Approach

*   **Print Statements:** Liberal use of `print()` statements within the Python functions during development and testing. This is invaluable for tracing the execution flow, inspecting intermediate variable values (e.g., individual `dq` terms, effective geometry values, per-link/per-node time steps), and understanding how these contribute to the final result.
*   **Python Debugger (`pdb`):** For more complex scenarios, especially within iterative loops like in `perform_picard_iteration_for_link_system`, utilize `pdb` or IDE-specific debuggers. This allows setting breakpoints, executing code step-by-step, examining the call stack, and inspecting variable values at specific points in execution.
*   **Modular Testing:** Test functions in isolation before testing functions that depend on them. For instance, ensure `get_preissmann_slot_width_py` is correct before thoroughly testing `get_effective_area_py`, and ensure all geometry and momentum equation functions are reliable before focusing heavily on `perform_picard_iteration_for_link_system`.
*   **Comparison with Conceptual Logic:** Continuously refer to the `SWMM_DynamicWave_Conduit_Details.md` document (which this testing strategy is part of) and, if necessary, the original SWMM C code's structure and comments. Discrepancies in Python output should first prompt a review of the Python code against these conceptual references to ensure the translated algorithm aligns with the intended logic.
*   **Simplify and Isolate:** If a complex function (e.g., Picard iteration) shows unexpected behavior, simplify its inputs significantly or temporarily hardcode values for some of its internal calculations (or for functions it calls). This helps to isolate the problematic part of the logic. For example, when testing Picard iteration, one might initially use fixed, hand-calculated geometry values instead of dynamically calling the slot geometry functions.

## 5. Limitations to State

It is important to reiterate the scope and limitations of this testing strategy:
*   This testing is designed to verify the translated Python algorithms **in isolation**, focusing on their internal logic and mathematical correctness relative to the source C algorithms.
*   It **does not verify the integration** of these Python functions into a larger, comprehensive hydraulic simulation system.
*   The testing **does not involve comparison against full SWMM model outputs**. Such a comparison would require a complete Python port of SWMM, which is beyond the current scope.
*   The Python code, particularly in helper functions or simplified test setups (like the rectangular conduit geometry callbacks), uses **simplified representations** for many aspects that are more complex in the full SWMM application (e.g., detailed cross-section geometry for various shapes, a full range of node types, all possible link types and their specific behaviors). The core translated algorithms, however, are intended to be faithful to the C logic.

By following this strategy, we aim to build confidence in the correctness of the Python translations of these key SWMM dynamic wave algorithms.
