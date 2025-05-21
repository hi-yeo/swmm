## Discretization of St. Venant Equations for Conduit Flow in SWMM

This section details how the EPA Stormwater Management Model (SWMM), specifically within the `dwflow_findConduitFlow` function (located in `dwflow.c`), discretizes the St. Venant equations to model dynamic wave flow in conduits. The explanation is aimed at engineers familiar with basic hydraulics but not necessarily numerical methods experts, focusing on the implementation as seen in the code.

### 1. Starting Point: St. Venant Equations

Dynamic wave routing in SWMM is the most sophisticated method for simulating flow, as it solves the full one-dimensional St. Venant equations for unsteady open channel flow. These equations consist of:
1.  **Continuity Equation:** Mathematically expresses the conservation of mass (water volume) along the conduit.
2.  **Momentum Equation:** Mathematically expresses the conservation of momentum for the water flowing in the conduit, accounting for various forces acting on it (gravity, pressure, friction, inertia).

While the overall dynamic wave solution (`dynwave.c`) handles the coupled system of these equations for the entire network of nodes and links, the `dwflow_findConduitFlow` function plays a crucial role in this process.

### 2. Focus on Momentum Equation for a Conduit

The function `dwflow_findConduitFlow` is called iteratively within the larger dynamic wave solution framework for each conduit (link `j`) at each time step `dt`. Its primary responsibility is to assemble the terms of a finite difference approximation of the **momentum equation** for that single conduit and then calculate an updated trial flow `q` for it. This trial flow is then used by the broader algorithm in `dynwave.c` to update water levels (heads) at the connected nodes. This iterative process (node head updates and link flow updates) continues until the solution converges for the current time step.

### 3. Key Variables

The function utilizes several key variables to represent the hydraulic conditions at the upstream (node 1, often denoted with suffix `1`) and downstream (node 2, often denoted with suffix `2`) ends of the conduit, as well as average or mid-conduit conditions. These include:

*   **Node-related (at conduit ends):**
    *   `h1`, `h2`: Hydraulic heads at the upstream and downstream nodes connected to the conduit (ft).
    *   `y1`, `y2`: Flow depths within the conduit at its upstream and downstream ends, respectively (ft). These are calculated as `h1 - z1` and `h2 - z2`, where `z1` and `z2` are the conduit's invert elevations at each end.
    *   `a1`, `a2`: Cross-sectional flow areas at the upstream and downstream ends of the conduit (ft²), corresponding to depths `y1` and `y2`.

*   **Conduit-averaged/representative:**
    *   `aMid`: Average or mid-point cross-sectional flow area in the conduit (ft²).
    *   `rMid`: Average or mid-point hydraulic radius in the conduit (ft).
    *   `r1`: Hydraulic radius at the upstream end of the conduit (ft).
    *   `v`: Representative velocity in the conduit, typically calculated as `qLast / aMid` (ft/s).

*   **Flow variables:**
    *   `q`: The new trial flow rate being calculated for the current iteration (cfs).
    *   `qOld`: The flow rate in the conduit from the *previous time step* (cfs). This is `Link[j].oldFlow / barrels` (where `barrels` is the number of identical parallel conduit barrels).
    *   `qLast`: The flow rate in the conduit from the *previous iteration* within the current time step (cfs). This is `Conduit[k].q1` (where `k` is the conduit's sub-index).

*   **Other important variables:**
    *   `aOld`: Cross-sectional area from the *previous time step* (ft²), specifically `Conduit[k].a2`.
    *   `dt`: The current simulation time step duration (s).
    *   `length`: The effective length of the conduit used in calculations (`Conduit[k].modLength`) (ft).
    *   `barrels`: Number of identical parallel barrels for the conduit.

### 4. Momentum Equation Terms

The core of `dwflow_findConduitFlow` involves calculating several terms, labeled `dq1` through `dq6` in the code. These terms represent the individual components of the discretized momentum equation, reflecting the various forces and effects influencing the flow:

*   **`dq1`: Friction Slope Term**
    *   **Code (for non-force main):** `dt * Conduit[k].roughFactor / pow(rWtd, 1.33333) * fabs(v)`
    *   **Physical Meaning:** This term quantifies the momentum loss due to friction between the flowing water and the wetted perimeter of the conduit. It's essentially derived from Manning's equation for open channels. For force mains (pressurized pipes), a different formulation based on Hazen-Williams or Darcy-Weisbach is used via `forcemain_getFricSlope()`.
    *   `Conduit[k].roughFactor` incorporates Manning's roughness coefficient 'n'.
    *   `rWtd` is an upstream-weighted hydraulic radius, and `fabs(v)` is the absolute value of the flow velocity. The hydraulic radius is raised to the power of 4/3 (1.33333) as per Manning's formula structure when solving for friction slope.
    *   The term is multiplied by `dt` as it contributes to the change in momentum over the current time step.

*   **`dq2`: Pressure Gradient (or Water Surface Slope) Term**
    *   **Code:** `dt * GRAVITY * aWtd * (h2 - h1) / length`
    *   **Physical Meaning:** This term represents the net force exerted on the water body in the conduit due to the difference in hydrostatic pressure between its downstream (`h2`) and upstream (`h1`) ends. This pressure difference is directly related to the difference in water surface elevations.
    *   `aWtd` is the upstream-weighted flow area. `GRAVITY` is the gravitational acceleration.
    *   If `h2 > h1`, the downstream head is higher, creating an adverse pressure gradient that opposes forward flow (or drives reverse flow), hence its subtractive effect in the flow update formula (as `-dq2`).

*   **`dq3`: Inertial Term (Local Acceleration)**
    *   **Code:** `2.0 * v * (aMid - aOld) * sigma`
    *   **Physical Meaning:** This is the local acceleration component of inertia. It accounts for the change in momentum due to the change in flow velocity that occurs if the average flow area (`aMid`) changes over time (from `aOld` at the previous time step).
    *   `v` is the current velocity. `sigma` is an inertial damping factor.
    *   This term is only included if `sigma > 0.0`.

*   **`dq4`: Inertial Term (Convective Acceleration)**
    *   **Code:** `dt * v * v * (a2 - a1) / length * sigma`
    *   **Physical Meaning:** This is the convective acceleration component of inertia. It accounts for the change in momentum due to the change in flow velocity as water moves from a section with area `a1` to a section with area `a2` along the conduit's `length`.
    *   `v*v` is the velocity squared. `sigma` is an inertial damping factor.
    *   This term is also only included if `sigma > 0.0`.

*   **`dq5`: Local Losses Term**
    *   **Code:** `findLocalLosses(j, a1, a2, aMid, qLast) / 2.0 / length * dt`
    *   **Physical Meaning:** This term accounts for additional energy (momentum) losses that occur due to localized hydraulic disturbances, such as at conduit entrances (`Link[j].cLossInlet`), exits (`Link[j].cLossOutlet`), or other average losses along the conduit (`Link[j].cLossAvg`). These are typically quantified using empirical loss coefficients multiplied by a velocity head. The `findLocalLosses` function calculates the combined effect.
    *   This term is only computed if the conduit is defined as having such losses (`Conduit[k].hasLosses`).

*   **`dq6`: Term for Evaporation and Seepage Losses**
    *   **Code:** `link_getLossRate(j, DW, qLast, dt) * 2.5 * dt * v / link_getLength(j)`
    *   **Physical Meaning:** This term represents the effect of water volume lost per unit length of the conduit due to surface evaporation or seepage through the conduit walls. While these are mass losses, they are incorporated into the momentum equation as a sink term that can affect the flow rate.
    *   `link_getLossRate` fetches the user-defined evaporation and seepage rates. The additional factors (`2.5`, `v`, `link_getLength(j)`) are part of SWMM's empirical formulation to translate these volume loss rates into an equivalent momentum effect.

### 5. Flow Update Formula

Once all the `dq` terms (representing various components of force and momentum change per unit mass or volume) are calculated, the new trial flow rate `q` for the conduit in the current iteration is determined by the following explicit algebraic formula:

`q = (qOld - dq2 + dq3 + dq4 + dq6) / (1.0 + dq1 + dq5)`

**Explanation:**
*   `qOld` (flow from the *previous time step*) serves as the initial estimate of flow for the current time step.
*   The **numerator** adjusts this initial flow:
    *   `- dq2`: Accounts for the pressure gradient (if `h2 > h1`, `dq2` is positive, reducing `q`).
    *   `+ dq3 + dq4`: Adds the effects of local and convective acceleration (if areas/velocities change appropriately).
    *   `+ dq6`: Adjusts for flow lost due to evaporation/seepage (this term is often small or zero).
*   The **denominator** `(1.0 + dq1 + dq5)` incorporates terms that typically resist flow:
    *   `dq1`: Friction losses.
    *   `dq5`: Local losses.
    A larger friction or local loss term increases the denominator, thereby reducing the magnitude of the calculated flow `q`.

This calculated `q` is a *trial flow*. Within the `dwflow_findConduitFlow` function, this `q` is then further adjusted by an under-relaxation factor `omega` based on `qLast` (flow from the previous iteration) to improve the stability of the iterative solution:
`q = (1.0 - omega) * qLast + omega * q;` (if `steps > 0`)
This is part of a Picard iteration method, common for solving non-linear systems. The overall `dynwave.c` module uses these iteratively refined conduit flows to update node heads, repeating the process until convergence across the network for the time step `dt`.

### 6. Numerical Scheme Aspects (Conceptual)

Several numerical techniques are evident in `dwflow_findConduitFlow` to enhance the stability and robustness of the solution:

*   **Upstream Weighting (`rho`):**
    *   The calculation of `aWtd` (weighted area) and `rWtd` (weighted hydraulic radius) involves a factor `rho`:
        `rho = 1.0; if ( !isFull && qLast > 0.0 && h1 >= h2 ) rho = sigma;`
        `aWtd = a1 + (aMid - a1) * rho;`
        `rWtd = r1 + (rMid - r1) * rho;`
    *   `rho` is influenced by the inertial damping factor `sigma`. When flow is positive, the conduit isn't full, and the upstream head is higher (normal flow conditions), `rho` takes the value of `sigma`. Otherwise, `rho` is 1.0.
    *   This weighting scheme adjusts how the spatial averages of area and hydraulic radius (used primarily in the friction `dq1` and pressure gradient `dq2` terms) are computed. Giving more weight to upstream conditions (`a1`, `r1` when `rho` is closer to 1.0, or when `sigma` is 1.0) is a common technique to improve stability in finite difference methods, especially when advective forces are strong (e.g., supercritical flow).

*   **Inertial Damping (`sigma`):**
    *   `sigma` is a factor (0 to 1) that scales the inertial terms (`dq3` and `dq4`). It's calculated based on the Froude number: `sigma` is 1.0 (no damping) for Froude numbers <= 0.5, 0.0 (full damping) for Froude numbers >= 1.0, and varies linearly in between.
    *   Users can also set global inertial damping options (`InertDamping`: `NO_DAMPING`, `PARTIAL_DAMPING`, `FULL_DAMPING`). If `FULL_DAMPING` is chosen, or if a closed conduit is surcharged, `sigma` is set to 0.0, effectively removing the inertial terms.
    *   This damping is crucial for numerical stability, particularly in or near critical/supercritical flow regimes, or during rapid flow changes, where the inertial terms can otherwise cause oscillations and prevent convergence.

*   **Time Discretization:**
    *   The use of `qOld` (flow at the end of the previous time step, `t`) as the primary component in the numerator `(qOld - dq2 ...)` indicates that the scheme is essentially an explicit update for flow `q` at the new time `t + dt`.
    *   The terms `dq1` through `dq6` represent discretized approximations of the various terms in the momentum equation, averaged or evaluated over the time interval `dt`. The areas and depths used are from the current iteration's node head solution.

### 7. Conduit Length (`length`)

*   The variable `length` used in the denominator of the pressure gradient (`dq2`) and convective acceleration (`dq4`) terms is assigned `Conduit[k].modLength`.
*   The code comments state: "--- use Courant-modified length instead of conduit's actual length".
*   The Courant-Friedrichs-Lewy (CFL) condition is a fundamental stability requirement for explicit numerical solutions of hyperbolic partial differential equations, like the St. Venant equations. It links the time step (`dt`), the spatial discretization (conduit length `L`), and the wave celerity (`c`) such that `c * dt / L <= 1`.
*   By using a `modLength` (modified length), SWMM can internally adjust the effective spatial discretization used in the numerical scheme. This allows the model to maintain stability for a user-selected `dt` across conduits of varying physical lengths and under different flow conditions, effectively helping to satisfy the CFL condition without requiring the user to manually adjust `dt` to extremely small values for short conduits. This is a practical approach to make the model more robust for a wider range of input scenarios.

In essence, `dwflow_findConduitFlow` implements a carefully constructed finite difference approximation to the momentum equation. It linearizes non-linear terms (like friction) and uses iterative refinement along with numerical stabilization techniques (upstream weighting, inertial damping, Courant-modified length) to compute a stable and sufficiently accurate conduit flow within each time step of the dynamic wave simulation.

### 8. Picard Iteration for Solving the System

The dynamic wave solution in SWMM uses an iterative approach, specifically Picard iteration (also known as the method of successive approximations), to solve the non-linear system of equations that arise from discretizing the St. Venant equations over a time step `dt`.

#### 8.1 Purpose of Iteration

The St. Venant equations are non-linear, primarily due to:
*   The friction term in the momentum equation (which depends on velocity squared, or `Q|Q|/R^(4/3)`).
*   The convective acceleration terms (which involve velocity squared).
*   The fact that flow area `A` and hydraulic radius `R` are non-linear functions of depth `y`, and depth `y` at nodes is coupled with flows `Q` in links.

Because of these non-linearities and interdependencies, a direct analytical solution for the flows and depths at the end of a time step (`t + dt`) is generally not possible. Instead, an iterative numerical method is required to converge on a solution that satisfies the discretized equations for that time step. SWMM employs Picard iteration for this purpose.

#### 8.2 Overall Iteration Loop (managed by `dynwave_execute` in `dynwave.c`)

The function `dynwave_execute(double tStep)` in `dynwave.c` orchestrates the iterative solution process for each dynamic wave routing time step (`tStep`).
*   It initializes a loop that continues until either convergence is achieved across the network or a maximum number of iterations (`MaxTrials`, typically 8) is reached.
*   A variable `Steps` (static within `dynwave.c`) tracks the current iteration count within the time step.
*   The global variable `Omega` (static within `dynwave.c`, typically initialized to `OMEGA = 0.5`) is used for under-relaxation.

The basic structure of this loop is:
```c
// Simplified from dynwave_execute
Steps = 0;
converged = FALSE;
Omega = OMEGA; // Default under-relaxation parameter
initRoutingStep(); // Initializes conduit areas from previous time step

while ( Steps < MaxTrials )
{
    initNodeStates();      // Initialize node surface areas, inflows, outflows
    findLinkFlows(tStep);  // Calculate trial flows for all links
    converged = findNodeDepths(tStep); // Update node depths and check convergence
    Steps++;
    if ( Steps > 1 )
    {
        if ( converged ) break;
        findBypassedLinks(); // Optimization for converged links/nodes
    }
}
if ( !converged ) updateConvergenceStats();
findLimitedLinks();
```

#### 8.3 Calculating Trial Link Flows (`findLinkFlows` and `dwflow_findConduitFlow`)

Within each iteration of the `dynwave_execute` loop:
1.  `findLinkFlows(tStep)` is called. This function iterates through all links in the network.
2.  For each "true conduit" (a standard conduit, not a dummy link, pump, or regulator), `findLinkFlows` calls `dwflow_findConduitFlow(link_index, Steps, Omega, tStep)`.
3.  As detailed in Section 4 and 5, `dwflow_findConduitFlow` computes a new trial flow `q` for the conduit using the explicit formula based on conditions from the previous time step (`qOld`, `aOld`) and the current iteration's node heads (`h1`, `h2`).
4.  **Under-Relaxation:** If it's not the first iteration (`Steps > 0`), the calculated trial flow `q` is immediately modified by under-relaxation within `dwflow_findConduitFlow`:
    ```c
    // In dwflow_findConduitFlow:
    if ( steps > 0 ) // 'steps' here is the 'Steps' from dynwave_execute
    {
        q = (1.0 - omega) * qLast + omega * q;
        // 'omega' here is the 'Omega' passed from dynwave_execute
        // 'qLast' is Conduit[k].q1, the flow from the previous iteration
        if ( q * qLast < 0.0 ) q = 0.001 * SGN(q); // Prevents oscillation
    }
    ```
    This step is crucial for stability. The newly computed `q` (which might be very different from the previous iteration's `qLast`) is blended with `qLast`. With `omega = 0.5`, it's a simple average. This dampens oscillations and helps prevent the solution from diverging, especially in complex networks or rapidly changing conditions.

#### 8.4 Updating Node Depths (`findNodeDepths` and `setNodeDepth`)

After `findLinkFlows` has computed new trial flows (including under-relaxation) for all links:
1.  `dynwave_execute` calls `findNodeDepths(tStep)`.
2.  `findNodeDepths` iterates through all non-outfall nodes. For each node, it calls `setNodeDepth(node_index, tStep)`.
3.  `setNodeDepth` calculates the net inflow (`dQ = Node[i].inflow - Node[i].outflow`) to the node using the link flows computed in `findLinkFlows`. It then estimates the change in node volume `dV` over the time step.
4.  Based on `dV` and the node's surface area (`surfArea`), it calculates a change in depth `dy` and updates the node's `newDepth`.
    *   For non-surcharged nodes: `dy = dV / surfArea; yNew = yOld + dy;` followed by under-relaxation if `Steps > 0`: `yNew = (1.0 - Omega) * yLast + Omega * yNew;`
    *   For surcharged nodes (EXTRAN method), a different approach using `sumdqdh` (sum of flow-head derivatives from connecting links) is used to estimate `dy`, and under-relaxation is generally not applied directly in the same manner for depth.
5.  The `newDepth` at each node effectively represents the solution of the continuity equation at that node for the current iteration.

#### 8.5 Convergence Check (within `findNodeDepths`)

After updating all node depths:
1.  `findNodeDepths` checks if the solution has converged for the current iteration.
2.  For each non-outfall node, it compares the newly computed `Node[i].newDepth` with the depth from the previous iteration (`yOld`, which was stored before calling `setNodeDepth`).
3.  If `fabs(yOld - Node[i].newDepth) > HeadTol` for any non-outfall node, that node is marked as not converged (`Xnode[i].converged = FALSE`). `HeadTol` is a user-specifiable head tolerance (default is 0.005 ft).
4.  `findNodeDepths` returns `TRUE` if all non-outfall nodes have converged, and `FALSE` otherwise.
5.  If `converged` is `TRUE`, the `while` loop in `dynwave_execute` terminates, and the current time step is considered successfully solved. Otherwise, the loop continues to the next iteration (up to `MaxTrials`).

If `MaxTrials` is reached and convergence is still not achieved, `updateConvergenceStats()` is called to log information about non-converging nodes/links.

#### 8.6 Flowchart of the Iteration Process

```mermaid
graph TD
    DYN_START[Start dynwave_execute(tStep)] --> DYN_INIT_STEP[initRoutingStep(): Store prev. conduit areas (a2=a1)];
    DYN_INIT_STEP --> DYN_LOOP_COND{Steps < MaxTrials AND NOT Converged?};
    DYN_LOOP_COND -- Yes --> DYN_INIT_NODES[initNodeStates(): Init node SA, In/Outflow, sumdqdh];
    DYN_INIT_NODES --> DYN_LINK_FLOWS[findLinkFlows(tStep)];
    DYN_LINK_FLOWS --> DYN_DWFLOW[For each True Conduit (if not bypassed):\n call dwflow_findConduitFlow(Steps, Omega, tStep)\n  - Calculates trial q\n  - Applies under-relaxation: q = (1-Omega)*qLast + Omega*q];
    DYN_DWFLOW --> DYN_NONCONDUIT[For other links (Pumps, Regulators):\n call findNonConduitFlow()];
    DYN_NONCONDUIT --> DYN_UPDATE_NODEFLOWS[updateNodeFlows() for all links:\n Accumulate Node[n].inflow/outflow, Xnode[n].sumdqdh];
    DYN_UPDATE_NODEFLOWS --> DYN_NODE_DEPTHS[findNodeDepths(tStep)];
    DYN_NODE_DEPTHS --> DYN_SET_DEPTH[For each Non-Outfall Node:\n call setNodeDepth(dt)\n  - Updates Node[].newDepth\n  - Applies under-relaxation for non-surcharged node depths];
    DYN_SET_DEPTH --> DYN_CHECK_CONV[Check Convergence:\n fabs(yOld - newDepth) <= HeadTol for all nodes?];
    DYN_CHECK_CONV --> DYN_SET_CONV_FLAG[Set 'Converged' flag];
    DYN_SET_CONV_FLAG --> DYN_INC_STEPS[Steps++];
    DYN_INC_STEPS --> DYN_CHECK_BREAK{Steps > 1 AND Converged?};
    DYN_CHECK_BREAK -- Yes --> DYN_LOOP_END[End Iteration Loop];
    DYN_CHECK_BREAK -- No --> DYN_BYPASS[findBypassedLinks()];
    DYN_BYPASS --> DYN_LOOP_COND;
    DYN_LOOP_COND -- No --> DYN_LOOP_END;
    DYN_LOOP_END --> DYN_POST_CONV{NOT Converged?};
    DYN_POST_CONV -- Yes --> DYN_STATS[updateConvergenceStats()];
    DYN_POST_CONV -- No --> DYN_LIMITED_LINKS[findLimitedLinks()];
    DYN_STATS --> DYN_LIMITED_LINKS;
    DYN_LIMITED_LINKS --> DYN_END[End dynwave_execute];
```

This iterative process effectively decouples the solution of the momentum equation (primarily in `dwflow_findConduitFlow`) and the continuity equation (primarily in `setNodeDepth`) within each time step, using Picard iteration with under-relaxation to guide the entire system towards a converged solution for flows and depths across the network.

### 9. Preissmann Slot for Surcharged Flow

When a closed conduit becomes full and then surcharged (pressurized), the water surface disappears, and the standard St. Venant equations for open channel flow are no longer directly applicable. SWMM offers a numerical technique called the "Preissmann Slot" to handle this transition smoothly without switching to a completely different set of equations for pressurized flow.

#### 9.1 Purpose of the Preissmann Slot

*   **Unified Flow Equations:** The primary purpose is to allow the same St. Venant equations used for open channel flow to be applied to simulate full-flow and pressurized (surcharged) conditions in closed conduits. This avoids the complexity of implementing and switching between different mathematical models (e.g., open channel vs. rigid lid pressurized flow models).
*   **Hypothetical Slot:** The method introduces a narrow, hypothetical open slot at the crown (top) of the closed conduit. When the water level reaches the conduit crown, it can continue to rise into this imaginary slot. This ensures that a "free surface" always exists mathematically, even if the physical pipe is full. The pressure head in the surcharged pipe is then represented by the height of the water in this slot.
*   **Smooth Transition:** This approach allows for a continuous and numerically stable transition between open-channel, just-full, and surcharged flow conditions.

#### 9.2 Activation Conditions

*   **User Option:** The Preissmann Slot method is used if the user selects `SLOT` for the `SurchargeMethod` option in the SWMM input file's `[OPTIONS]` section.
*   **Depth Threshold:** For a specific closed conduit (e.g., circular, rectangular, but not open channels), the slot calculations become active when the flow depth `y` approaches or exceeds the conduit's full depth `yFull`.
    *   In `dwflow.c`, functions like `getSlotWidth()` and `getWidth()` use a global variable `CrownCutoff` to determine when to engage the slot. If `SurchargeMethod == SLOT`, `CrownCutoff` is initialized to `SLOT_CROWN_CUTOFF` (approximately 0.985257) in `dynwave_init()` (from `dynwave.c`).
    *   The condition `yNorm < CrownCutoff` (where `yNorm = y / xsect->yFull`) is checked. If `yNorm` is below this, the slot width is zero. If `yNorm` is at or above this cutoff, a non-zero slot width is calculated.
*   **Closed Conduits Only:** The Preissmann slot concept is only applicable to closed conduits that can surcharge. It is not used for open channels (e.g., `xsect_isOpen(xsect->type)` check in `getSlotWidth`).

#### 9.3 Calculating Slot Width (`getSlotWidth` in `dwflow.c`)

The function `getSlotWidth(TXsect* xsect, double y)` calculates the width of this hypothetical slot (`wSlot`).
*   **No Slot Conditions:**
    *   If `SurchargeMethod != SLOT`.
    *   If the conduit is an open type (`xsect_isOpen(xsect->type)` is true).
    *   If the normalized depth `yNorm = y / xsect->yFull` is less than `CrownCutoff`.
    In these cases, `wSlot` is 0.0.
*   **Sjoberg's Formula:** If the slot is active and `yNorm` is between `CrownCutoff` and a higher limit (1.78 in the code, i.e., `yNorm <= 1.78`), SWMM uses Sjoberg's formula to determine the slot width:
    `wSlot = xsect->wMax * 0.5423 * exp(-pow(yNorm, 2.4))`
    where `xsect->wMax` is the maximum physical top width of the conduit cross-section. This formula provides a slot width that decreases rapidly as the surcharge depth (and `yNorm`) increases, mimicking the rapid increase in wave speed (and thus decrease in effective "storage width") under pressurized conditions.
*   **Maximum Surcharge Depth:** If `yNorm` exceeds the higher limit (e.g., `yNorm > 1.78`), the slot width is capped at a very small percentage of the conduit's maximum physical width:
    `wSlot = 0.01 * xsect->wMax`
    This prevents the slot width from becoming excessively small or zero at very high surcharge depths, which could lead to numerical instability.

#### 9.4 Modifying Conduit Geometry for Flow Calculations

When `wSlot` is greater than zero, it modifies the effective geometric properties of the conduit used in the St. Venant equation terms:

*   **Top Width (`getWidth` in `dwflow.c`):**
    *   The function `getWidth(TXsect* xsect, double y)` first calls `getSlotWidth(xsect, y)`.
    *   If the returned `wSlot > 0.0`, then this `wSlot` is used as the effective top width of the flow.
    *   Otherwise (if `wSlot` is 0.0, meaning flow is below the crown or slot method is not used), the physical top width of the wetted cross-section at depth `y` is used (obtained via `xsect_getWofY(xsect, y)`). There's an additional check: if the depth `y` is very near the crown (`y / xsect->yFull >= CrownCutoff`) and the conduit is closed, the depth `y` used for `xsect_getWofY` is capped at `CrownCutoff * xsect->yFull`. This is likely to ensure a consistent top width for the final segment of open channel flow before the slot activates.

*   **Flow Area (`getArea` in `dwflow.c`):**
    *   The function `getArea(TXsect* xsect, double y, double wSlot)` (where `wSlot` is passed from `getSlotWidth`) calculates the effective flow area.
    *   If the depth `y` is greater than or equal to the conduit's full depth `yFull` (i.e., `y >= xsect->yFull`):
        `Area = xsect->aFull + (y - xsect->yFull) * wSlot`
        This calculation adds the area of a conceptual rectangle (the slot) on top of the full physical pipe area. The height of this rectangle is `(y - xsect->yFull)`, and its width is `wSlot`.
    *   If `y < xsect->yFull`, the standard physical cross-sectional area `xsect_getAofY(xsect, y)` is returned.

#### 9.5 Impact on Hydraulic Radius (`getHydRad` in `dwflow.c`)

*   The function `getHydRad(TXsect* xsect, double y)` determines the hydraulic radius.
*   If the depth `y` is greater than or equal to the full depth `yFull` (`y >= xsect->yFull`):
    `Hydraulic Radius = xsect->rFull`
    In this case, SWMM uses the hydraulic radius of the just-full pipe. It does *not* attempt to calculate a modified hydraulic radius that would include the wetted perimeter of the narrow slot walls. This is a common and practical simplification in Preissmann slot implementations, as the slot is primarily a numerical device to handle the transition to pressurized flow and maintain a free surface for the equations, rather than a physical feature whose frictional properties need detailed representation.
*   If `y < xsect->yFull`, the standard physical hydraulic radius `xsect_getRofY(xsect, y)` is returned.

#### 9.6 Conceptual Logic Flow for Applying Preissmann Slot

```mermaid
graph TD
    PS_START[Start Geometry Calculation for Conduit Link] --> PS_CHECKS{SurchargeMethod == SLOT AND Conduit is Closed Type?};
    PS_CHECKS -- No --> PS_NO_SLOT[Use Standard Physical Geometry (Area, Top Width, Hyd. Radius) based on depth y];
    PS_CHECKS -- Yes --> PS_CALC_YNORM[Calculate yNorm = y / yFull];
    PS_CALC_YNORM --> PS_CHECK_CROWNCUTOFF{yNorm >= CrownCutoff?};
    PS_CHECK_CROWNCUTOFF -- No --> PS_NO_SLOT;
    PS_CHECK_CROWNCUTOFF -- Yes --> PS_CALC_WSLOT[Calculate wSlot (Sjoberg or min value)];
    PS_CALC_WSLOT --> PS_MOD_GEOM[Apply Modified Geometry:];
    PS_MOD_GEOM --> PS_MOD_TOPWIDTH["Effective Top Width = wSlot"];
    PS_MOD_GEOM --> PS_MOD_AREA["Effective Area = aFull + (y - yFull) * wSlot"];
    PS_MOD_GEOM --> PS_MOD_HYDRAD["Hydraulic Radius = rFull"];
    PS_NO_SLOT --> PS_END[End Geometry Calculation];
    PS_MOD_HYDRAD --> PS_END;
```

**Textual Explanation of the Flowchart:**
1.  **Check Global Setting and Conduit Type:** First, verify if the `SurchargeMethod` is set to `SLOT` and if the conduit is a closed type (e.g., not an open channel). If not, standard physical geometry calculations are used.
2.  **Calculate Normalized Depth:** Compute `yNorm = y / yFull`.
3.  **Check Crown Cutoff:**
    *   If `yNorm` is less than `CrownCutoff` (e.g., ~0.985), the flow is considered open channel within the physical pipe. `wSlot` is 0. Standard physical geometry is used.
    *   If `yNorm` is greater than or equal to `CrownCutoff`, the Preissmann slot logic is engaged.
4.  **Calculate Slot Width (`wSlot`):**
    *   If `CrownCutoff <= yNorm <= 1.78` (approx.), `wSlot` is calculated using Sjoberg's formula.
    *   If `yNorm > 1.78` (approx.), `wSlot` is set to a minimum value (e.g., 1% of max physical width).
5.  **Apply Modified Geometry:**
    *   The **effective top width** for the St. Venant equations becomes `wSlot`.
    *   The **effective flow area** becomes `aFull + (y - yFull) * wSlot`.
    *   The **effective hydraulic radius** is taken as `rFull` (the hydraulic radius of the pipe when just full).
6.  These modified geometric properties are then used in the discretization of the St. Venant equations (specifically in calculating terms like `dq1`, `dq2`, `dq4` within `dwflow_findConduitFlow`).

By artificially maintaining a free surface through the Preissmann slot, SWMM can use a single set of equations (St. Venant) to model the complex transition from open channel to full, surcharged flow, and back again. This greatly simplifies the numerical model while providing a stable and generally effective way to handle pressurized conditions in pipe networks.

### 10. Dynamic Wave Time Step (Δt) Adjustment Logic

SWMM employs an adaptive time step adjustment mechanism when the "Variable Time Step" option (controlled by the `CourantFactor` in the `[OPTIONS]` section of the input file) is enabled for dynamic wave routing. This logic, primarily managed within `dynwave.c`, aims to ensure numerical stability while optimizing computational speed.

#### 10.1 Purpose of Variable Time Stepping

*   **Courant-Friedrichs-Lewy (CFL) Condition:** Explicit numerical solutions for hyperbolic Partial Differential Equations (PDEs), like the St. Venant equations, are subject to stability criteria. The most well-known is the CFL condition. It states that the numerical domain of dependence must contain the true domain of dependence. For the St. Venant equations, this means the time step `Δt` must be limited relative to the conduit length `L` and the speed `c` at which information (waves) can travel through the water.
*   **Simplified CFL Representation:** A common representation is `c * Δt / L <= 1.0`. SWMM uses a `CourantFactor` (a user-defined safety factor, typically between 0.5 and 1.0, found in the `[OPTIONS]` section as `COURANT_FACTOR`) in its criterion, effectively `(c / L) * Δt <= CourantFactor` or `Δt <= CourantFactor * L / c`.
*   **Stability and Efficiency:** If `Δt` is too large for a given `L` and `c`, the numerical solution can become unstable, leading to oscillations and unrealistic results. Conversely, an unnecessarily small `Δt` increases computation time without significantly improving accuracy.
*   **Adaptive Approach:** Since wave celerity `c` (which depends on flow depth and velocity) and flow conditions vary dynamically across the network and over time, a fixed `Δt` might be too large for some conditions/locations and too small for others. Variable time stepping allows SWMM to:
    *   Automatically reduce `Δt` when conditions become hydraulically "fast" (high velocities, rapid depth changes) to maintain stability.
    *   Increase `Δt` when conditions are "slow" to speed up the simulation, up to a user-defined maximum (`fixedStep` or `RouteStep`).

#### 10.2 Main Controlling Functions (`dynwave.c`)

*   **`dynwave_getRoutingStep(double fixedStep)`:** This is the primary function called by the main simulation control loop (in `swmm_step()` within `swmm5.c`) at the beginning of each potential routing period. `fixedStep` is the maximum routing time step specified by the user (e.g., the `ROUTING_STEP` from the `[OPTIONS]` section of the input file).
    *   **Fixed Time Step Check:** If `CourantFactor == 0.0` (meaning the variable time step option is disabled by the user), or if `fixedStep` is already very small (less than `MINTIMESTEP` which is 0.001s), this function simply returns `fixedStep`.
    *   **Variable Time Step Logic:** Otherwise, it proceeds to calculate a suitable variable time step. If it's the very start of the simulation (`VariableStep == 0.0`), it initializes `VariableStep` to `MinRouteStep` (user-defined minimum routing step, or a default). For subsequent calls, it invokes `getVariableStep(fixedStep)` to compute the new `VariableStep`.
    *   **Millisecond Adjustment:** The final `VariableStep` is adjusted to be a multiple of a millisecond: `floor(1000.0 * VariableStep) / 1000.0`.

*   **`getVariableStep(double maxStep)`:** This function is the core of the adaptive time step calculation.
    *   It aims to find a `tMin` (minimum stable time step) that satisfies stability criteria for both links (conduits) and nodes.
    *   It ensures `tMin` does not exceed `maxStep` (the user's maximum routing step) and does not fall below `MinRouteStep` (user-defined minimum, or `MINTIMESTEP`).
    *   The process involves:
        1.  Initializing `tMin = maxStep`.
        2.  Calling `getLinkStep(tMin, &minLink)` to get `tMinLink` based on Courant conditions in conduits.
        3.  Calling `getNodeStep(tMinLink, &minNode)` to get `tMinNode` based on rate of depth change at nodes, using `tMinLink` as an upper bound.
        4.  The final `tMin` before returning is the lesser of `tMinLink` and `tMinNode`, further constrained by `MinRouteStep`.

#### 10.3 Link Time Step Calculation (`getLinkStep`)

The function `getLinkStep(double tMin, int *minLink)` calculates the smallest time step (`tLink`) required across all "true" conduits to satisfy a form of the Courant condition.
*   It iterates through each conduit link `i`.
*   It skips conduits with negligible flow (`q <= FUDGE`), very small area (`Conduit[k].a1 <= FUDGE`), or very low Froude number (`Link[i].froude <= 0.01`), as these are unlikely to govern stability.
*   For eligible conduits, it calculates a trial time step `t` using the formula:
    `t = (Link[i].newVolume / Conduit[k].barrels / q) * (Conduit[k].modLength / link_getLength(i)) * (Link[i].froude / (1.0 + Link[i].froude)) * CourantFactor`

    *   **Term 1: `(Link[i].newVolume / Conduit[k].barrels / q)`:** This is `Volume / FlowRate`, which simplifies to `(Area * Length) / (Velocity * Area) = Length / Velocity = L/V`. This represents the travel time of water through the conduit based on its actual physical length (`link_getLength(i)` which is `Conduit[k].length`).
    *   **Term 2: `(Conduit[k].modLength / link_getLength(i))`:** This is the ratio of the "Courant-modified length" (`modLength`) to the actual physical length. As discussed in Section 7, `modLength` is an adjusted length used internally by SWMM to help satisfy the CFL condition over a range of conduit sizes for a given time step.
    *   **Term 3: `(Link[i].froude / (1.0 + Link[i].froude))`:** This term involves the Froude number `Fr = V / sqrt(gD)` (where D is hydraulic depth). For `Fr << 1` (subcritical flow), this term is approx. `Fr`. For `Fr >> 1` (supercritical flow), this term approaches 1. This factor attempts to relate the advective speed `V` to the gravity wave speed `c = V ± sqrt(gD)`. The term `V / (Fr / (1+Fr))` acts as an effective wave celerity in the Courant formulation.
    *   **Term 4: `CourantFactor`:** The user-defined safety factor (typically <= 1.0).

    The formula can be conceptually interpreted as `Δt_link = (L_actual / V_effective_wave) * (modLength / L_actual) * CourantFactor = (modLength / V_effective_wave) * CourantFactor`. The `modLength` helps normalize behavior across different physical lengths.
*   The function updates `tLink = MIN(tLink, t)` and stores the index of the link that resulted in this minimum in `*minLink`.

#### 10.4 Node Time Step Calculation (`getNodeStep`)

The function `getNodeStep(double tMin, int *minNode)` calculates the smallest time step (`tNode`) required to prevent excessive depth changes at nodes within a single step. This is a heuristic criterion rather than a strict CFL condition.
*   It iterates through all non-outfall nodes.
*   It skips nodes with very shallow depth, or nodes whose depth is already very close to their crown elevation (likely already stable or full).
*   For eligible nodes, it defines a maximum allowable depth change in one step:
    `maxDepthChange = (Node[i].crownElev - Node[i].invertElev) * 0.25`
    This is 25% of the available height from the node invert to its crown (top of highest connecting pipe).
*   It retrieves `dYdT = Xnode[i].dYdT`, which is the rate of change of depth at the node calculated in the previous iteration's `setNodeDepth` function (`fabs(yNew - yOld) / dt_previous_iteration`).
*   If `dYdT` is significant, it calculates a time step `t1` that would result in `maxDepthChange`:
    `t1 = maxDepthChange / dYdT`
*   The function updates `tNode = MIN(tNode, t1)` and stores the index of the node that resulted in this minimum in `*minNode`. The initial `tMin` passed to this function (which is `tMinLink`) acts as an upper bound for `tNode`.

#### 10.5 Final Time Step Selection and Adjustment

1.  **Minimum of Link and Node Steps:** In `getVariableStep`, `tMin` is first set to `maxStep` (the overall maximum routing step). It's then updated by `getLinkStep` to `tMinLink`. This `tMinLink` is then passed to `getNodeStep`, which may further reduce it to `tMinNode`. So, `tMin` emerging from these two calls is effectively `MIN(maxStep, tMinLink, tMinNode)`.
2.  **Enforce Absolute Minimum:** This `tMin` is then compared against `MinRouteStep` (a user-defined absolute minimum, defaulting to `MINTIMESTEP` if not set or if set lower).
    `tMin = MAX(MinRouteStep, tMin)`
    This ensures the time step doesn't become excessively small, which could lead to very long simulation times or potential floating-point precision issues.
3.  **Return from `getVariableStep`:** This potentially adjusted `tMin` is returned as `VariableStep`.
4.  **Millisecond Flooring (in `dynwave_getRoutingStep`):** The `VariableStep` obtained is then floored to the nearest millisecond:
    `VariableStep = floor(1000.0 * VariableStep) / 1000.0;`
    This is likely done for consistency in reporting or to avoid issues with very fine variations in time steps.

The global variable `MinRouteStep` (default 0.001s if not user-set or if user set lower) and `CourantFactor` (from `[OPTIONS]`) are key parameters governing this adaptive time-stepping behavior.

#### 10.6 Flowchart of Time Step Adjustment Logic

```mermaid
graph TD
    DGRS_START[Start dynwave_getRoutingStep(fixedStep)];
    DGRS_START --> DGRS_CHECK_CF{CourantFactor == 0.0 OR\nfixedStep < MINTIMESTEP?};
    DGRS_CHECK_CF -- Yes --> DGRS_RETURN_FIXED[Return fixedStep];
    DGRS_CHECK_CF -- No --> DGRS_INIT_VARSTEP{VariableStep == 0.0? (Start of Sim)};
    DGRS_INIT_VARSTEP -- Yes --> DGRS_SET_MINROUTE[VariableStep = MinRouteStep];
    DGRS_INIT_VARSTEP -- No --> DGRS_CALL_GETVAR[VariableStep = getVariableStep(fixedStep)];
    DGRS_SET_MINROUTE --> DGRS_CALL_GETVAR;
    DGRS_CALL_GETVAR --> DGRS_FLOOR[VariableStep = floor(1000.0 * VariableStep) / 1000.0];
    DGRS_FLOOR --> DGRS_RETURN_VAR[Return VariableStep];

subgraph getVariableStep(maxStep)
    direction LR
    GVS_START[Start] --> GVS_INIT_TMIN[tMin = maxStep];
    GVS_INIT_TMIN --> GVS_CALL_LINKSTEP[tMinLink = getLinkStep(tMin, &minLink)];
    GVS_CALL_LINKSTEP --> GVS_CALL_NODESTEP[tMinNode = getNodeStep(tMinLink, &minNode)];
    GVS_CALL_NODESTEP --> GVS_COMPARE_LINKNODE[tMin = MIN(tMinLink, tMinNode)];
    GVS_COMPARE_LINKNODE --> GVS_CHECK_MINROUTE[tMin = MAX(MinRouteStep, tMin)];
    GVS_CHECK_MINROUTE --> GVS_RETURN[Return tMin];
end

subgraph getLinkStep(tMin_in, &minLink)
    direction LR
    GLS_START[Start] --> GLS_INIT_TLINK[tLink = tMin_in];
    GLS_INIT_TLINK --> GLS_LOOP_LINKS{For each Conduit Link i};
    GLS_LOOP_LINKS -- Next Link --> GLS_CHECK_COND{Flow, Area, Froude > Min?};
    GLS_CHECK_COND -- Yes --> GLS_CALC_T[Calculate Courant t for Link i:\n t = (Vol/Q) * (modL/L) * (Fr/(1+Fr)) * CF];
    GLS_CALC_T --> GLS_UPDATE_TLINK{t < tLink?};
    GLS_UPDATE_TLINK -- Yes --> GLS_SET_TLINK[tLink = t, *minLink = i];
    GLS_UPDATE_TLINK -- No --> GLS_LOOP_LINKS;
    GLS_SET_TLINK --> GLS_LOOP_LINKS;
    GLS_CHECK_COND -- No --> GLS_LOOP_LINKS;
    GLS_LOOP_LINKS -- End Loop --> GLS_RETURN[Return tLink];
end

subgraph getNodeStep(tMin_in, &minNode)
    direction LR
    GNS_START[Start] --> GNS_INIT_TNODE[tNode = tMin_in];
    GNS_INIT_TNODE --> GNS_LOOP_NODES{For each Non-Outfall Node i};
    GNS_LOOP_NODES -- Next Node --> GNS_CHECK_COND{Depth & dYdT > Min?};
    GNS_CHECK_COND -- Yes --> GNS_CALC_MAXDY[maxDy = (CrownElev - InvertElev) * 0.25];
    GNS_CALC_MAXDY --> GNS_CALC_T1[t1 = maxDy / dYdT];
    GNS_CALC_T1 --> GNS_UPDATE_TNODE{t1 < tNode?};
    GNS_UPDATE_TNODE -- Yes --> GNS_SET_TNODE[tNode = t1, *minNode = i];
    GNS_UPDATE_TNODE -- No --> GNS_LOOP_NODES;
    GNS_SET_TNODE --> GNS_LOOP_NODES;
    GNS_CHECK_COND -- No --> GNS_LOOP_NODES;
    GNS_LOOP_NODES -- End Loop --> GNS_RETURN[Return tNode];
end

```

This adaptive time-stepping logic is critical for the stability and efficiency of SWMM's dynamic wave routing module, allowing it to handle a wide range of hydraulic conditions and network configurations.
