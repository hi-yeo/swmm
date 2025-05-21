## EPA SWMM: Core Concepts and Program Flow

This document provides a high-level overview of the EPA Stormwater Management Model (SWMM) codebase, its purpose, and its overall program flow. It's intended for engineers who are new to SWMM or its source code.

### 1. Purpose of SWMM

The Stormwater Management Model (SWMM) is a powerful simulation tool developed by the U.S. Environmental Protection Agency (EPA). As stated in its `README.md`:

> "SWMM is a dynamic hydrology-hydraulic water quality simulation model. It is used for single event or long-term (continuous) simulation of runoff quantity and quality from primarily urban areas."

In essence, SWMM helps engineers and planners predict how rainfall and snowmelt will affect urban drainage systems. It can model:
*   **Runoff:** How much water flows off surfaces like roads, roofs, and green spaces.
*   **Hydraulics:** How water moves through pipes, channels, storage units, and pumps.
*   **Water Quality:** How pollutants are transported and transformed within the urban water system.

This allows users to design and analyze stormwater infrastructure, assess the effectiveness of water quality improvement measures (like green infrastructure), and manage flood risks.

### 2. Main Inputs & Outputs

SWMM interacts with the user primarily through text-based files.

*   **Inputs:**
    *   **Main Input File (`.inp`):** This is the primary file that defines the entire simulation scenario. It's a plain text file containing various sections that describe:
        *   **Site Characteristics:** General information about the study area, simulation options, and overall parameters.
        *   **Subcatchment Data:** Details about individual land areas, including their size, slope, imperviousness (e.g., percentage of paved surfaces), soil properties, and how they are connected to the drainage network.
        *   **Drainage Network Layout:** Information about the manholes (junctions), pipes (conduits), channels, pumps, storage units, and outfalls that make up the stormwater system. This includes their physical dimensions, materials, and connectivity.
        *   **Rainfall Data References:** While some simple rainfall data can be directly entered, the `.inp` file typically references external rainfall data files. These files contain time series of precipitation intensity (e.g., inches per hour) at specific rain gages.
        *   **Pollutant Information:** Definitions of pollutants to be simulated and their buildup/washoff characteristics on land surfaces.
        *   **Treatment Information:** Descriptions of any treatment processes occurring within the drainage system elements.
    *   **External Data Files:** These can include:
        *   **Rainfall Files:** Time series of rainfall data.
        *   **Temperature Files:** For snowmelt calculations.
        *   **Inflow Files:** External hydrographs or pollutographs that enter the system at specific locations.
        *   **Control Rule Files:** Rules that define how elements like pumps or weirs operate based on conditions in the system.

*   **Outputs:**
    Based on the command-line arguments in `src/run/main.c` (`runswmm f1 f2 f3`), SWMM produces two main output files:
    *   **Report File (`.rpt`):** This is a human-readable text file that summarizes the simulation. It includes:
        *   A copy of the input data.
        *   Status messages, warnings, and errors generated during the simulation.
        *   Summary tables for key results (e.g., peak flows, total volumes, pollutant loads) for subcatchments, nodes, and links.
        *   Detailed time series results for specific elements if requested by the user in the input file.
    *   **Binary Output File (`.out`):** This is a binary (non-human-readable) file that stores detailed time-series results for *all* elements in the model at each reporting time step. This file is typically used for:
        *   Post-processing with other software tools.
        *   Visualization and detailed analysis using SWMM's graphical user interface (GUI) or other specialized programs.
        *   Providing initial conditions for subsequent simulations.
        This file is optional and will only be created if a filename is provided on the command line.

### 3. Flowchart Description for `main.c` (from `src/run/main.c`)

The `main.c` file provides the command-line interface for running SWMM simulations. Its program flow is as follows:

```mermaid
graph TD
    A[Start] --> B{Parse Command Line Arguments};
    B --> C{Number of Arguments?};
    C -- 1 --> D[Print "Not Enough Arguments"];
    D --> X[Exit];
    C -- 2 --> E{Argument Type?};
    E -- "--help" or "-h" --> F[Display Help Text];
    F --> X;
    E -- "--version" or "-v" --> G[Display Version];
    G --> X;
    E -- Other --> H[Print "Unknown Argument"];
    H --> X;
    C -- "3 or more" --> I[Extract Input, Report, Binary File Names];
    I --> J[Print SWMM Version];
    J --> K[Call swmm_run(inputFile, reportFile, binaryFile)];
    K --> L[Record End Time];
    L --> M[Calculate & Print Run Time];
    M --> N{Check for Errors/Warnings};
    N -- Errors/Warnings --> O[Print Status Message];
    O --> P[End];
    N -- No Errors/Warnings --> P;
```

**Explanation of `main.c` Flow:**

1.  **Start:** The program begins.
2.  **Parse Command Line Arguments:** It examines the arguments provided when `runswmm` is executed from the command line (e.g., `runswmm project.inp project.rpt project.out`).
3.  **Check Number of Arguments:**
    *   If only one argument (just `runswmm`) is given, it prints "Not Enough Arguments" and exits.
    *   If two arguments are given:
        *   It checks if the argument is `--help` (or `-h`) or `--version` (or `-v`).
        *   If `--help`, it displays usage instructions and exits.
        *   If `--version`, it displays the SWMM build version and exits.
        *   Otherwise, it's an "Unknown Argument", and the program prints an error and exits.
    *   If three or more arguments are provided:
        *   The first argument is taken as the input file path (`.inp`).
        *   The second is the report file path (`.rpt`).
        *   The third (if present) is the binary output file path (`.out`). If not provided, no binary output is generated.
        *   It prints the SWMM version being used.
        *   Crucially, it calls `swmm_run()`. This function is the entry point into the main SWMM simulation engine (which is often compiled as a separate library, like `swmm5.dll`). This is where all the hydrological and hydraulic calculations occur.
        *   After `swmm_run()` completes, it records the end time.
        *   It calculates the total simulation run time and prints it to the console.
        *   It then checks if `swmm_run()` reported any errors (via `swmm_getError()`) or warnings (via `swmm_getWarnings()`) and prints an appropriate status message (e.g., "There are errors.", "There are warnings.", or a clean completion message).
4.  **End:** The program terminates.

### 4. Flowchart Description for `swmm_run()` (Conceptual)

The `swmm_run()` function is the heart of the SWMM engine. It orchestrates the entire simulation process. While the actual implementation involves many C functions and modules, the conceptual flow can be visualized as follows:

```mermaid
graph TD
    SR_A[Start swmm_run(inputFile, reportFile, binaryFile)] --> SR_B[Call swmm_open()];
    SR_B --"Internally"--> SR_B1[project_open(): Setup project data structures];
    SR_B --"Internally"--> SR_B2[report_open(): Prepare report file];
    SR_B --"Internally"--> SR_B3[output_open(): Prepare binary output file (if specified)];
    
    SR_B --> SR_C[Call project_readInput(inputFile)];
    SR_C --"Tasks"--> SR_C1[Read & Parse .inp file];
    SR_C --"Tasks"--> SR_C2[Populate internal data structures (subcatchments, nodes, links, etc.)];
    SR_C --"Tasks"--> SR_C3[Validate input data];

    SR_C --> SR_D[Initialize Simulation State];
    SR_D --"Tasks"--> SR_D1[Perform initial calculations];
    SR_D --"Tasks"--> SR_D2[Set initial conditions (e.g., water levels)];
    SR_D --"Tasks"--> SR_D3[stats_open(): Initialize statistics collection];

    SR_D --> SR_E[Call swmm_start() - Initiate Simulation Loop];
    SR_E --> SR_F{Loop Through Simulation Time Steps};
    SR_F -- Each Time Step --> SR_G[Determine Current Time Step Size];
    SR_G --> SR_H[Execute Runoff Calculations (runoff_execute)];
    SR_H --"For each subcatchment"--> SR_H1[Calculate precipitation, snowmelt, infiltration, runoff, pollutant washoff];
    
    SR_H --> SR_I[Execute Routing Calculations (routing_execute)];
    SR_I --"Through drainage network"--> SR_I1[Route flows and pollutants (nodes & links)];
    SR_I --"Involves"--> SR_I2[Hydraulic modeling (e.g., kinematic/dynamic wave)];
    SR_I --"Involves"--> SR_I3[Water quality routing];
    SR_I --"Involves"--> SR_I4[Control actions (pumps, weirs)];

    SR_I --> SR_J[Update Overall System State];
    SR_J --> SR_K[Collect Time Series Results & Statistics];
    SR_K --> SR_F;
    SR_F -- End of Loop --> SR_L[Call swmm_end()];
    SR_L --"Tasks"--> SR_L1[Finalize simulation (e.g., stats_close() to compute summary stats)];

    SR_L --> SR_M[Call swmm_report()];
    SR_M --"Tasks"--> SR_M1[Generate text-based report file (.rpt)];
    SR_M1 --"Includes"--> SR_M1a[Input data summaries];
    SR_M1a --"Includes"--> SR_M1b[Simulation statistics];
    SR_M1b --"Includes"--> SR_M1c[Detailed results for selected objects];

    SR_M --> SR_N[Call swmm_close()];
    SR_N --"Tasks"--> SR_N1[Clean up: Close files, free memory (project_close, output_close)];
    
    SR_N --> SR_Z[End swmm_run];
```

**Explanation of `swmm_run()` Conceptual Flow:**

1.  **Start `swmm_run(inputFile, reportFile, binaryFile)`:** The function is called with the file paths obtained from the command line.
2.  **Call `swmm_open()`:** This is the first major step.
    *   Internally, this involves functions like `project_open()` which set up the main data structures to hold all the information about the simulation (subcatchments, network, etc.).
    *   It also prepares the report file (`report_open()`) and the binary output file (`output_open()`) for writing, if a binary file was specified.
3.  **Call `project_readInput(inputFile)`:** (Corresponds to `input.c` module)
    *   This function is responsible for reading the entire `.inp` file section by section.
    *   It parses the data and populates the internal data structures created in the previous step (e.g., arrays of subcatchment objects, node objects, link objects).
    *   It also performs validation checks on the input data to ensure consistency and correctness.
4.  **Initialize Simulation State:**
    *   Before the simulation loop starts, various initial calculations are performed, and initial conditions are set (e.g., initial water levels in storage units, initial pollutant concentrations).
    *   Statistics collection is also initialized (e.g., via `stats_open()` from `stats.c`).
5.  **Call `swmm_start()` - Initiate Simulation Loop:** This function kicks off the core of the simulation.
    *   **Loop through simulation time steps:** The simulation proceeds step-by-step from a defined start time to an end time. The duration of each time step can be fixed or variable (adaptive time stepping based on simulation conditions).
        *   **Determine current time step size.**
        *   **Execute Runoff Calculations** (e.g., using `runoff_execute()` from `runoff.c`): For each subcatchment defined in the input:
            *   It calculates how much rain falls, how much snow melts (if applicable), how much water infiltrates into the ground, how much groundwater flows, how much water runs off the surface, and how pollutants are washed off by the runoff.
        *   **Execute Routing Calculations** (e.g., using `routing_execute()` from `routing.c`): This part handles the movement of water and pollutants through the drainage network (pipes, channels, etc.).
            *   It routes the runoff generated from subcatchments, along with any other external inflows, through the system of nodes (junctions, storage units, outfalls) and links (conduits, pumps, weirs).
            *   This involves complex hydraulic calculations which can use different physics models (like kinematic wave or the more complex dynamic wave for pressurized flow and backwater effects).
            *   It also routes pollutants through the network and can simulate control actions (e.g., turning a pump on or off based on water level).
        *   **Update overall system state:** After runoff and routing for the current time step, the overall state of the system (water levels, flow rates, pollutant concentrations) is updated.
        *   **Collect time series results and statistics:** Data for the current time step is written to the binary output file (if active) and aggregated for summary statistics (using functions from `output.c` and `stats.c`).
    *   **End of time step loop:** The loop continues until the simulation end time is reached.
6.  **Call `swmm_end()`:**
    *   Once the simulation loop is complete, this function finalizes the simulation. For example, `stats_close()` might be called to compute final summary statistics from the data collected during the simulation.
7.  **Call `swmm_report()`:**
    *   This function generates the human-readable report file (`.rpt`). It uses various modules (`report.c`, `statsrpt.c`, `inputrpt.c`) to write:
        *   Summaries of the input data.
        *   Overall simulation statistics (e.g., continuity errors, peak values).
        *   Detailed results (tables, time series) for specific objects as requested by the user in the `.inp` file.
8.  **Call `swmm_close()`:**
    *   This is the final cleanup step. It closes all opened files (input, report, binary output) and frees any dynamically allocated memory (e.g., `project_close()`, `output_close()`).
9.  **End `swmm_run`:** The main simulation processing is complete, and control returns to `main()`.

This overview should provide a foundational understanding of what SWMM does and how its command-line version processes a simulation. The actual codebase is extensive and involves many more specialized modules and functions for each step described.

### 5. Key Module Categories and Representative Modules

This section details the key categories of modules within the SWMM engine and provides examples of specific C files within those categories. This categorization is based on a conceptual understanding of SWMM's architecture, as a `Roadmap.txt` file detailing this explicitly was not available. The roles and responsibilities are inferred from common practices in hydrological simulation software and the functions generally associated with these module names.

### 5.1 Core Modules

*   **Overall Purpose:** These modules form the backbone of the SWMM engine. They are responsible for managing the overall simulation process, handling the project's data structures, managing input and output operations, and defining the global objects and enumerations used throughout the program.
*   **Representative Module 1: `input.c`**
    *   **Role:** Its primary responsibility is to read and parse the user-provided `.inp` (input) file, which contains all the definitions and parameters for a simulation scenario.
    *   **Key Functions/Responsibilities (Conceptual):**
        *   Reading data section by section (e.g., `[OPTIONS]`, `[SUBCATCHMENTS]`, `[CONDUITS]`).
        *   Parsing keywords, identifiers (names of objects), and numerical data for each type of object or setting.
        *   Storing the parsed data into the appropriate internal data structures. These structures are typically defined in header files like `objects.h` and their memory is managed by functions likely found in a module like `project.c`.
        *   Performing basic validation of input values, such as checking for required data fields, ensuring values are within acceptable ranges, and verifying that object names are unique.
        *   Establishing connections and relationships between different model components (e.g., linking subcatchments to their outlet nodes, or connecting nodes together via links like conduits or pumps).
*   **Representative Module 2: `report.c`**
    *   **Role:** Its primary responsibility is to generate the human-readable text report file (`.rpt`) at the end of a successful simulation run. This file provides a summary of the simulation's setup, its performance, and key results.
    *   **Key Functions/Responsibilities (Conceptual):**
        *   Writing summary tables for various simulation results. This can include subcatchment runoff summaries (total precipitation, infiltration, runoff, peak flows), node summaries (average and peak depths, overflow occurrences), link flow summaries (peak flows, capacities), and overall mass balance/continuity errors.
        *   Formatting and writing detailed time-series output for specific objects if the user has requested this in the input file.
        *   Interacting with other modules to gather information for the report. For example, it might call functions from `statsrpt.c` (for detailed simulation statistics) and `inputrpt.c` (for echoing and summarizing the input data within the report).

### 5.2 Runoff Calculation Modules

*   **Overall Purpose:** These modules are responsible for simulating all the hydrologic processes that occur on the land surface (subcatchments). Their main goal is to calculate the quantity (flow) and quality (pollutant concentrations) of stormwater runoff generated from these areas, which then becomes input to the conveyance system.
*   **Representative Module 1: `subcatch.c`**
    *   **Role:** This is a central module for subcatchment-level runoff calculations. It orchestrates the computation of runoff generated from individual subcatchments, along with associated pollutant buildup during dry weather and washoff during storm events.
    *   **Key Functions/Responsibilities (Conceptual):**
        *   Iterating through each subcatchment defined in the model at each runoff calculation time step.
        *   Calculating the total precipitation input to the subcatchment, considering data from assigned rain gages and accounting for any snow accumulation or melt.
        *   Calling functions from other specialized modules to determine hydrological abstractions (losses):
            *   `infil.c`: To calculate infiltration of water into pervious land surfaces (e.g., using Horton's, Green-Ampt, or Curve Number methods).
            *   `snow.c`: To manage snowpack accumulation, redistribution, and melting processes if snow is part of the simulation.
        *   Calculating the amount of water stored on the surface (depression storage) and the resulting surface runoff. This often involves methods based on Manning's equation for overland flow or non-linear reservoir routing.
        *   Simulating pollutant buildup on the land surface during periods without rain.
        *   Simulating pollutant washoff from the land surface by the generated runoff, using various functions (e.g., exponential, rating curve).
        *   Interacting with `lid.c` / `lidproc.c` if Low Impact Development (LID) controls (like rain gardens or permeable pavement) are present within the subcatchment, as these will modify the runoff and pollutant generation processes.
        *   Aggregating the final runoff flow and pollutant concentrations from the subcatchment and passing these as lateral inflows to specific nodes (junctions or outfalls) in the drainage system model.

### 5.3 Flow and Water Quality Routing Modules

*   **Overall Purpose:** These modules handle the movement (routing) of water and pollutants through the conveyance network of the drainage system. This network consists of nodes (like manholes, storage units, and outfalls) and links (like pipes, channels, pumps, and weirs).
*   **Representative Module 1: `flowrout.c`**
    *   **Role:** This module provides the top-level control and logic for the flow routing process within the conveyance system. It determines which hydraulic routing method to use (as specified by the user, typically kinematic wave or dynamic wave) and manages the sequence of calculations required to simulate flow through the network.
    *   **Key Functions/Responsibilities (Conceptual):**
        *   Initializing the hydraulic state (e.g., depths at nodes, flows in links) at the beginning of each routing time step.
        *   Accumulating all inflows to each node in the system. These inflows can come from subcatchment runoff, direct user-defined inflows, dry weather sanitary flows, and rainfall-dependent infiltration/inflow (RDII) into the sewer system.
        *   Calling the appropriate core routing function based on the user's choice:
            *   `kinwave.c`: If kinematic wave routing is selected. This is a simpler method suitable for many systems where flow is primarily gravity-driven and downstream conditions don't significantly affect upstream flow.
            *   `dynwave.c`: If dynamic wave routing is selected. This is a more complex method that solves the full St. Venant equations.
        *   Updating node depths (or hydraulic heads) and link flows based on the results of the chosen routing calculation for the current time step.
        *   For dynamic wave routing, it plays a role in managing the iterative solution process, checking for convergence, and adjusting time steps if necessary to maintain stability.
*   **Representative Module 2: `dynwave.c`**
    *   **Role:** Implements the dynamic wave hydraulic routing method. This method solves the full one-dimensional St. Venant equations of unsteady flow, which consist of the continuity and momentum equations. It is computationally more intensive but allows for modeling complex hydraulic phenomena such as backwater effects, flow reversals, pressurized (surcharged) flow, and looped networks.
    *   **Key Functions/Responsibilities (Conceptual):**
        *   Discretizing the St. Venant equations (partial differential equations) into a system of algebraic equations for each link (conduit) in the network, typically using a finite difference scheme.
        *   Setting up a large system of non-linear equations that represent the flow and head (water level) conditions across the entire network for the current time step.
        *   Iteratively solving this system of equations. This often involves matrix solution techniques (e.g., variants of Newton's method) to find the unknown node heads and link flows that satisfy the equations.
        *   Handling various boundary conditions, such as water levels at outfall nodes, operational status of pumps, and settings of weirs or orifices.
        *   May interact with lower-level modules (e.g., `dwflow.c` or similar) that contain the specific mathematical formulations for solving the flow equations for individual conduits or handling the matrix operations.

### 5.4 Support Modules

*   **Overall Purpose:** These modules provide various utility functions, data management tools, and specialized calculations that are used by many other parts of the SWMM engine. They don't typically represent a core simulation process on their own but are essential for the functioning of the entire system.
*   **Representative Module 1: `datetime.c`**
    *   **Role:** Provides a suite of functions for handling dates, times, and time steps. These are fundamental to a time-based simulation model like SWMM, which tracks processes occurring over a simulation period.
    *   **Key Functions/Responsibilities (Conceptual):**
        *   Converting dates and times between different string representations (e.g., "MM/DD/YYYY HH:MM:SS") and internal numerical formats (like Julian date/time or seconds elapsed).
        *   Performing arithmetic operations on dates and times, such as adding or subtracting time durations (e.g., adding a time step to the current simulation time).
        *   Comparing different dates and times to determine sequence or duration.
        *   Managing the simulation clock, including advancing it by the current time step, and ensuring it aligns with reporting intervals and the overall simulation duration.
*   **Representative Module 2: `mathexpr.c`**
    *   **Role:** Parses and evaluates mathematical expressions that can be supplied by the user in certain sections of the input file. This allows for more flexible and customized model behavior, for example, to define complex pump operating curves, control rules for hydraulic structures, or time-varying pollutant concentrations.
    *   **Key Functions/Responsibilities (Conceptual):**
        *   Tokenizing the input mathematical expression string: breaking it down into individual components like numbers, variables, operators, and function calls.
        *   Parsing the sequence of tokens to build an internal representation of the expression, often as an expression tree or a postfix (Reverse Polish Notation) stack.
        *   Evaluating the expression at runtime during the simulation. This involves substituting the current values of simulation variables (like node depth, link flow, or current time) into the expression and computing the result.
        *   Handling a library of built-in mathematical functions (e.g., sin, cos, log, sqrt).

### 6. Flowchart Descriptions for Core Computational Processes

This section details the conceptual logic for the two main computational engines within SWMM: the runoff process (how water and pollutants are generated from land surfaces) and the routing process (how water and pollutants move through the drainage network).

### 6.1 Runoff Process (`runoff_execute`)

*   **Overall Purpose:** This process simulates all land-surface hydrology and water quality for each subcatchment in the model. It determines how much rainfall becomes runoff, how much infiltrates into the ground, how snow accumulates and melts, and how pollutants are picked up by the flowing water. The output from this process becomes the input to the hydraulic routing module.
*   **Inputs (for each subcatchment, at each runoff time step):**
    *   Precipitation data (intensity from the assigned rain gage, considering time series in rainfall files).
    *   Subcatchment properties (physical characteristics like area, width, slope; hydrologic parameters like imperviousness, depression storage depths, infiltration model parameters – all originally read by `input.c` and stored in data structures defined in `objects.h`).
    *   Current snowpack conditions (depth, water content, etc., if snowmelt is being modeled, managed by `snow.c`).
    *   Groundwater conditions (water table elevation, soil moisture, if groundwater interaction is modeled, managed by `gwater.c`).
    *   LID (Low Impact Development) control status and parameters (type of LID, dimensions, material properties, if LIDs are present in the subcatchment, managed by `lid.c`).
*   **Key Calculation Steps (Conceptual Flowchart):**
    ```mermaid
    graph TD
        R_Start[Start runoff_execute for a subcatchment] --> R_A[Get Precipitation];
        R_A --> R_B{Snow Present/Possible?};
        R_B -- Yes --> R_C[Calculate Snow Accumulation/Melt (snow.c)];
        R_C --> R_D[Net Rainfall/Snowmelt at Surface];
        R_B -- No --> R_D;
        R_D --> R_E[Calculate Depression Storage Losses (fill/spill)];
        R_E --> R_F[Water available for Infiltration & Runoff (Excess Precip)];
        R_F --> R_G{Pervious Area Component?};
        R_G -- Yes --> R_H[Calculate Infiltration into Soil (infil.c)];
        R_H --> R_I[Surface Runoff (Pervious Portion)];
        R_G -- No (Fully Impervious) --> R_I;
        R_F --> R_J[Surface Runoff (Impervious Portion, after its depression storage)];
        R_I --> R_K[Combine Pervious & Impervious Runoff Components];
        R_J --> R_K;
        R_K --> R_L{LID Controls Present in Subcatchment?};
        R_L -- Yes --> R_M[Process Combined Runoff through LIDs (lidproc.c)];
        R_M --> R_N[Final Subcatchment Runoff (Post-LID)];
        R_L -- No --> R_N[Final Subcatchment Runoff (Pre-LID, if no LIDs)];
        R_N --> R_O[Calculate Pollutant Buildup/Washoff (subcatch.c / landuse.c)];
        R_O --> R_P[Total Runoff Flow & Quality for Subcatchment];
        R_P --> R_End[End for subcatchment: Output to routing module];
    ```
*   **Textual Explanation of the Flowchart:**
    1.  **Start `runoff_execute` for a subcatchment:** The process begins for a specific subcatchment at the current runoff time step.
    2.  **Get Precipitation:** The current precipitation rate (rain or snow) is obtained from the relevant rain gage data.
    3.  **Snow Present/Possible?:** Check if conditions are suitable for snow processes (e.g., temperature, existing snowpack).
        *   **Yes:** Call `snow.c` functions to calculate snow accumulation on the subcatchment's snowpack or melt from the existing snowpack. This will adjust the amount of water available at the surface.
        *   **No:** Proceed with the precipitation as liquid rain.
    4.  **Net Rainfall/Snowmelt at Surface:** This is the actual liquid water available at the land surface after accounting for immediate snow processes.
    5.  **Calculate Depression Storage Losses:** Water first fills up small surface depressions. Only after this storage is filled does runoff begin.
    6.  **Water available for Infiltration & Runoff (Excess Precip):** This is the water remaining after depression storage is satisfied.
    7.  **Pervious Area Component?:** Subcatchments are often modeled as a combination of pervious (e.g., grass, soil) and impervious (e.g., pavement, roofs) areas. This step checks for the pervious component.
        *   **Yes:** For the pervious portion of the subcatchment, call `infil.c` functions to calculate how much water infiltrates into the soil based on the chosen infiltration model (e.g., Horton, Green-Ampt, Curve Number) and current soil moisture conditions.
        *   **Surface Runoff (Pervious Portion):** The water that doesn't infiltrate becomes surface runoff from the pervious area.
        *   **No (Fully Impervious):** If the subcatchment (or the portion being analyzed) is entirely impervious, infiltration is typically zero or very small, and most water becomes runoff directly.
    8.  **Surface Runoff (Impervious Portion, after its depression storage):** For the impervious portion of the subcatchment, water fills its depression storage, and the excess becomes surface runoff.
    9.  **Combine Pervious & Impervious Runoff Components:** The runoff generated from the pervious and impervious parts of the subcatchment are combined to get the total initial surface runoff.
    10. **LID Controls Present in Subcatchment?:** Check if any Low Impact Development controls (like rain gardens, permeable pavement, etc.) are defined for this subcatchment.
        *   **Yes:** Call `lidproc.c` functions. These functions will take the combined runoff and simulate its passage through the LID units. LIDs can reduce runoff volume (through infiltration and evapotranspiration within the LID) and alter the timing of the runoff.
        *   **No:** The combined runoff bypasses LID processing.
    11. **Final Subcatchment Runoff (Post-LID or Pre-LID):** This is the runoff hydrograph (flow rate over time) leaving the subcatchment after all land surface processes, including LIDs, have been accounted for.
    12. **Calculate Pollutant Buildup/Washoff:** Based on the land uses defined for the subcatchment and the generated runoff, calculate the buildup of pollutants during dry weather and their washoff during wet weather. This involves functions within `subcatch.c` and potentially `landuse.c`.
    13. **Total Runoff Flow & Quality for Subcatchment:** The final runoff flow rate and the associated pollutant concentrations (pollutograph) are now determined.
    14. **End for subcatchment: Output to routing module:** These runoff hydrographs and pollutographs are passed as lateral inflows to the specified nodes (manholes or other inlets) in the conveyance network model.
*   **Outputs (for each subcatchment, sent to the routing module):**
    *   Runoff flow rate (as a time series, or hydrograph).
    *   Pollutant concentrations or loads (as a time series, or pollutograph).

### 6.2 Routing Process (`routing_execute`)

*   **Overall Purpose:** This process simulates the movement of water and pollutants through the conveyance network, which consists of nodes (such as manholes, storage units, and outfalls) and links (such as pipes, channels, pumps, weirs, and orifices). It determines water depths at nodes and flow rates in links throughout the system over time.
*   **Inputs (at each routing time step):**
    *   Runoff hydrographs/pollutographs from `runoff_execute` (these become lateral inflows at specified nodes).
    *   Other direct/external inflows (user-defined time series, dry weather flows from `inflow.c`, etc.).
    *   Drainage network topology and properties (characteristics of nodes like invert elevation, maximum depth; and links like length, roughness, shape, size – all originally from `input.c` and stored in `objects.h`).
    *   Initial hydraulic conditions (water depths at nodes, flows in links) carried over from the end of the previous routing time step.
    *   Control rules (operational rules for pumps, weirs, etc., defined by the user and processed by `controls.c`).
    *   User-selected routing method (typically Kinematic Wave or Dynamic Wave, specified in the options section of the `.inp` file).
*   **Key Calculation Steps (Conceptual Flowchart):**
    ```mermaid
    graph TD
        RT_Start[Start routing_execute for the network] --> RT_A[Initialize Node/Link States for current time step];
        RT_A --> RT_B[Gather All Inflows to Nodes];
        RT_B -- Runoff from Subcatchments (runoff_execute) --> RT_B1;
        RT_B -- Direct/External Inflows (inflow.c) --> RT_B1;
        RT_B1[Total Inflow at each Node for current time step] --> RT_C{User-Selected Routing Method?};
        RT_C -- Kinematic Wave --> RT_D_KW[Calculate Link Flows (kinwave.c) based on upstream Node conditions, Link properties, & normal flow equations];
        RT_D_KW --> RT_E[Update Node Depths based on Continuity (Volume Change = Inflow - Outflow)];
        RT_C -- Dynamic Wave --> RT_D_DW[Iteratively Solve Full St. Venant Equations for all Links & Nodes (dynwave.c)];
        RT_D_DW -- Solution Converged --> RT_E;
        RT_E --> RT_F[Apply Control Rules (controls.c) - may adjust pump status, weir settings, etc.];
        RT_F --> RT_G[Calculate Pollutant Transport & Treatment through Network (qualrout.c, treatmnt.c)];
        RT_G --> RT_H[Update Final Node Depths, Link Flows, Water Quality for current time step];
        RT_H --> RT_End[End for network: Store results, prepare for next time step];
    ```
*   **Textual Explanation of the Flowchart:**
    1.  **Start `routing_execute` for the network:** The process begins for the entire conveyance system at the current routing time step.
    2.  **Initialize Node/Link States:** Set initial conditions for nodes (e.g., depth) and links (e.g., flow) based on the results from the end of the previous time step.
    3.  **Gather All Inflows to Nodes:** For each node in the network:
        *   Sum the lateral inflows from subcatchment runoff (output of `runoff_execute`).
        *   Add any other direct inflows like user-defined time series or dry weather flows (managed by `inflow.c`).
    4.  **Total Inflow at each Node for current time step:** This represents the total volume of water entering each node during this time step before any routing calculations.
    5.  **User-Selected Routing Method?:** SWMM checks which hydraulic routing method the user has chosen.
        *   **Kinematic Wave (`kinwave.c`):**
            *   **Calculate Link Flows:** For each link, the flow is calculated based on the hydraulic conditions at the upstream node (e.g., water depth), the link's geometric and hydraulic properties (shape, roughness, slope), and simplified normal flow equations (like Manning's equation). This method assumes flow is primarily driven by gravity and that downstream conditions do not affect upstream flow (no backwater).
            *   **Update Node Depths:** After link flows are determined, the change in water volume (and thus depth) at each node is calculated based on the principle of continuity: change in storage = (total inflow to node) - (total outflow from node via links).
        *   **Dynamic Wave (`dynwave.c`):**
            *   **Iteratively Solve Full St. Venant Equations:** This method solves the complete one-dimensional St. Venant equations (governing continuity and momentum) for all links and nodes in the network simultaneously. This is typically done using an implicit numerical method that forms a large system of non-linear equations. These equations are solved iteratively until a stable solution (convergence) for node depths and link flows is found for the current time step. Dynamic wave routing can simulate complex phenomena like backwater effects, flow reversals, and pressurized flow.
            *   **Solution Converged:** Once the iterative solver finds a solution, it proceeds to the next step.
    6.  **Apply Control Rules (`controls.c`):** After the primary hydraulic calculations, any user-defined control rules are evaluated. These rules might change the status of pumps (turn on/off), adjust weir opening heights, or modify orifice settings based on current simulated conditions (e.g., water depth in a node, flow in a conduit). If controls change operational parameters, it might necessitate a re-calculation or adjustment within the dynamic wave solver.
    7.  **Calculate Pollutant Transport & Treatment through Network (`qualrout.c`, `treatmnt.c`):** With the hydraulic conditions (flows and depths) established for the time step, the movement and transformation of pollutants are simulated. This involves:
        *   Advection and dispersion of pollutants within links.
        *   Mixing of pollutants at nodes.
        *   Applying any defined treatment processes (e.g., decay or removal) within links or nodes (using functions from `treatmnt.c`).
    8.  **Update Final Node Depths, Link Flows, Water Quality for current time step:** The final values for water depth/head at nodes, flow rate/velocity in links, and pollutant concentrations throughout the network are established for this time step.
    9.  **End for network: Store results, prepare for next time step:** These results are stored for reporting (to the `.rpt` file), for detailed output (to the binary `.out` file), and serve as the initial conditions for the next routing time step.
*   **Outputs (at each routing time step):**
    *   Updated hydraulic states for all nodes (e.g., water depth, hydraulic head, stored volume).
    *   Updated flow rates and other hydraulic parameters for all links (e.g., flow rate, flow velocity, water depth in link).
    *   Updated water quality concentrations in all nodes and links.
    *   Results are written to the binary output file and aggregated for the summary report file.

### 7. Guidance for a Novice Engineer

This section offers some practical advice and starting points for a novice engineer aiming to understand the intricacies of the SWMM C codebase. Diving into a large, established codebase can be daunting, but a structured approach can make it manageable and rewarding.

### 7.1 Recommended Learning Approach

*   **Start with the Big Picture:**
    *   Read the `README.md` for a general overview of SWMM's purpose.
    *   Familiarize yourself with this document (`SWMM_Explanation.md`) to understand the high-level structure, core concepts, and main program flow.
    *   Locate and review `src/solver/Roadmap.txt` (if available, or rely on the module descriptions in this document) for a module-by-module breakdown. *(Note: `Roadmap.txt` was not found in this repository during the generation of this document, so rely on Section 5 for module descriptions).*
*   **Understand the Entry Point and Core API:**
    *   Examine `src/run/main.c` to see how the command-line executable initializes and runs a simulation. Pay attention to how it calls `swmm_run()`.
    *   Look at `src/solver/swmm5.c` (and its corresponding header `swmm5.h`). This file defines the main API functions (`swmm_run`, `swmm_open`, `swmm_start`, `swmm_step`, `swmm_end`, `swmm_report`, `swmm_close`, etc.) that control the simulation. Understanding their sequence, as outlined in Section 4 of this document, is key.
*   **Explore Core Modules:**
    *   Dive into `input.c` (as described in Section 5.1) to see how data is read from the `.inp` file and stored.
    *   Look at `project.c` (conceptually related to `swmm_open` and data structure management) to understand how project data (subcatchments, nodes, links, etc., defined in `objects.h`) is organized and accessed.
    *   Review `report.c` (as described in Section 5.1) to see how the summary report is generated.
*   **Branch into Specific Processes based on Interest:**
    *   **Hydrology (Runoff):** If interested in how runoff is generated (Section 6.1), explore `runoff.c`, `subcatch.c` (Section 5.2), `infil.c`, `snow.c`, `gwater.c`, and `lid.c`/`lidproc.c`.
    *   **Hydraulics (Routing):** If interested in how water moves through the network (Section 6.2), explore `routing.c`, `flowrout.c` (Section 5.3), and then delve into `kinwave.c` or `dynwave.c` (Section 5.3) (and potentially `dwflow.c` for dynamic wave details) depending on the routing method of interest.
    *   **Water Quality:** Examine `qualrout.c` and `treatmnt.c`.
*   **Use the Code as a Reference:** When you see a function call, try to find its definition. Use code searching tools (like `grep` or IDE features) to navigate. The function prototypes in `funcs.h` and `swmm5.h` can also guide you.

### 7.2 Key Header Files

These header files are crucial for understanding data structures, global variables, constants, and function prototypes used throughout the SWMM engine. Regularly refer to them:

*   `src/solver/objects.h`: Defines the C structures for all major simulation objects (subcatchments, nodes, links, pollutants, etc.). Understanding these structures is fundamental.
*   `src/solver/enums.h`: Defines various enumerated types (symbolic constants) used to represent different options, object types, parameter codes, etc. This helps make the code more readable.
*   `src/solver/globals.h`: Declares global variables that are used across multiple modules. While extensive use of globals can make code harder to follow, in SWMM, they hold important simulation-wide data and state.
*   `src/solver/funcs.h`: Contains function prototypes for many general-purpose functions that are called from various modules.
*   `src/solver/consts.h`: Defines various numerical constants used in calculations.
*   `src/solver/text.h`: Defines text strings, often used for keywords in the input file or messages in the report file.
*   `src/solver/macros.h`: Contains useful macros that might simplify common coding patterns.
*   `src/solver/include/swmm5.h`: The public API header for the SWMM engine library, defining functions like `swmm_run`, `swmm_open`, etc.

### 7.3 Using the Command-Line Version for Tracing

*   **Build the Executable:** Follow the instructions in `Build.md` or `Readme.txt` (if available in the repository) to compile `runswmm` (the command-line executable, typically found in a directory like `src/run`).
*   **Find or Create Sample Input Files:** SWMM often comes with example `.inp` files (check for an "examples" directory or similar, or download from the official EPA SWMM website). Start with a small, simple example to minimize complexity.
*   **Run with a Debugger (Advanced):**
    *   If you are familiar with a C debugger (like GDB on Linux or the Visual Studio Debugger on Windows), running `runswmm` with a simple input file under the debugger is an invaluable way to learn.
    *   You can set breakpoints in key functions (e.g., in `src/run/main.c`, `src/solver/swmm5.c` functions like `swmm_run` or `swmm_step`, `runoff_execute` in `runoff.c`, `routing_execute` in `routing.c`, or specific functions in `input.c`) and step through the code line by line.
    *   This allows you to see the call stack (how functions call each other), inspect variable values (like those defined in `objects.h` or `globals.h`), and understand the program's state at different points in the simulation.
*   **Add Print Statements (Simpler):**
    *   As a simpler alternative to a debugger, you can temporarily add `printf` statements (e.g., `#include <stdio.h>` then `printf("Now in function_x with var_y = %f\n", var_y);`) at the beginning and end of functions, or to print specific variable values to the console. Recompile and run to see the execution flow or data changes. (Remember to remove these diagnostic prints before making any formal changes or contributions).

By combining these strategies—understanding the high-level concepts first, then diving into specific code modules, and using tools like debuggers or print statements to trace execution—a novice engineer can systematically explore and gain a solid understanding of the SWMM codebase. Be patient, take notes, and don't be afraid to jump between different files to follow the logic. Good luck!

### 8. Discretization of St. Venant Equations in `dwflow_findConduitFlow`

This section delves into how the SWMM engine, specifically within the `dwflow_findConduitFlow` function (located in `dwflow.c`), discretizes the St. Venant equations to model dynamic wave flow in conduits.

### 8.1 Starting Point: St. Venant Equations

Dynamic wave routing in SWMM solves the one-dimensional St. Venant equations for unsteady open channel flow. These equations consist of:
1.  **Continuity Equation:** Describes conservation of mass (water volume).
2.  **Momentum Equation:** Describes conservation of momentum (forces acting on the water).

While the overall dynamic wave solution (`dynwave.c`) handles the coupled system of equations for the entire network, the `dwflow_findConduitFlow` function focuses on calculating an updated flow for a *single conduit* within an iterative solution process. It primarily assembles and solves a finite difference approximation of the **momentum equation** for that conduit.

### 8.2 Focus on Momentum Equation for a Conduit

The function `dwflow_findConduitFlow` takes the current state of the system (heads at nodes, previous flow in the conduit) and calculates a new trial flow `q` for the conduit. This new flow will then be used in the broader iterative solution of the `dynwave.c` module to update node heads, and the process repeats until convergence for the current time step.

### 8.3 Key Variables

The function uses several variables representing conditions at the upstream (node 1, denoted with suffix `1`) and downstream (node 2, denoted with suffix `2`) ends of the conduit `j`, as well as average or mid-conduit conditions:

*   `h1`, `h2`: Heads at the upstream and downstream nodes (ft).
*   `y1`, `y2`: Flow depths in the conduit at the upstream and downstream ends (ft), calculated as `h1 - z1` and `h2 - z2` where `z1, z2` are conduit invert elevations.
*   `a1`, `a2`: Cross-sectional flow areas at the upstream and downstream ends (ft²).
*   `r1`: Hydraulic radius at the upstream end (ft).
*   `yMid`, `aMid`, `rMid`: Flow depth, area, and hydraulic radius at the mid-point or as an average for the conduit (ft, ft², ft).
*   `qOld`: Flow rate in the conduit from the *previous time step* (cfs). This is `Link[j].oldFlow / barrels`.
*   `qLast`: Flow rate in the conduit from the *previous iteration* of the current time step (cfs). This is `Conduit[k].q1`.
*   `aOld`: Cross-sectional area from the *previous time step* (ft²), specifically `Conduit[k].a2`.
*   `dt`: Current time step (s).
*   `length`: Effective conduit length (`Conduit[k].modLength`) (ft).
*   `v`: Velocity in the conduit, typically `qLast / aMid` (ft/s).

### 8.4 Momentum Equation Terms

The core of `dwflow_findConduitFlow` is calculating various terms (`dq1` through `dq6`) that represent the components of a finite difference form of the momentum equation. These terms are then used to solve for the new flow `q`.

The physical meaning of these terms:

*   **`dq1`: Friction Slope Term**
    *   `dq1 = dt * Conduit[k].roughFactor / pow(rWtd, 1.33333) * fabs(v)` (for non-force mains)
    *   This term represents the energy loss due to friction along the wetted perimeter of the conduit. It's derived from Manning's equation (or Darcy-Weisbach for force mains, handled by `forcemain_getFricSlope`).
    *   `Conduit[k].roughFactor` incorporates Manning's n, and `rWtd` is the upstream-weighted hydraulic radius. `fabs(v)` is the absolute velocity.
    *   It's proportional to `dt` because it influences the change in momentum over the time step.

*   **`dq2`: Pressure Gradient (or Energy Slope) Term**
    *   `dq2 = dt * GRAVITY * aWtd * (h2 - h1) / length`
    *   This term accounts for the net force due to the difference in water surface elevation (and thus pressure) between the downstream (`h2`) and upstream (`h1`) ends of theconduit.
    *   `aWtd` is the upstream-weighted flow area. `GRAVITY` is the acceleration due to gravity.
    *   A positive `h2 - h1` (downstream head higher) would mean an adverse pressure gradient, tending to reduce forward flow or cause reverse flow.

*   **`dq3`: Inertial Term (Temporal Acceleration)**
    *   `dq3 = 2.0 * v * (aMid - aOld) * sigma`
    *   This term represents the local inertia, i.e., the force required to change the flow velocity due to a change in the average flow area (`aMid`) over time (compared to `aOld` from the previous time step).
    *   `v` is the current velocity. `sigma` is an inertial damping factor.
    *   This term is only added if `sigma > 0.0`.

*   **`dq4`: Inertial Term (Convective Acceleration)**
    *   `dq4 = dt * v * v * (a2 - a1) / length * sigma`
    *   This term represents the convective inertia, i.e., the force required to change the flow velocity as it moves along the conduit due to a change in flow area from upstream (`a1`) to downstream (`a2`).
    *   `v*v` is velocity squared. `sigma` is an inertial damping factor.
    *   This term is only added if `sigma > 0.0`.

*   **`dq5`: Local Losses Term**
    *   `dq5 = findLocalLosses(j, a1, a2, aMid, qLast) / 2.0 / length * dt`
    *   This term accounts for energy losses due to abrupt changes in geometry or flow conditions, such as at the entrance (`Link[j].cLossInlet`), exit (`Link[j].cLossOutlet`), or along the conduit (`Link[j].cLossAvg`).
    *   The `findLocalLosses` function calculates these based on coefficients and velocities at different sections.
    *   This term is only added if `Conduit[k].hasLosses` is true.

*   **`dq6`: Evaporation and Seepage Losses Term**
    *   `dq6 = link_getLossRate(j, DW, qLast, dt) * 2.5 * dt * v / link_getLength(j)`
    *   This term accounts for water volume lost per unit length due to evaporation from the surface or seepage through the conduit walls.
    *   `link_getLossRate` retrieves these loss rates (which could be zero if not specified). The factor `2.5` and scaling by `v / link_getLength(j)` seem to be specific empirical adjustments or unit conversions within the SWMM formulation.

### 8.5 Flow Update Formula

After calculating all the `dq` terms, the new flow `q` for the conduit in the current iteration is computed using the following algebraic formula:

`q = (qOld - dq2 + dq3 + dq4 + dq6) / (1.0 + dq1 + dq5)`

Here:
*   `qOld` is the flow from the *previous time step*. This acts as the starting point for the current time step's calculation.
*   The numerator sums the flow from the previous time step and adjusts it based on the pressure gradient (`-dq2` because `dq2` is defined with `h2-h1`), inertial terms (`+dq3`, `+dq4`), and evaporation/seepage losses (`+dq6`).
*   The denominator (`1.0 + dq1 + dq5`) acts as a divisor that incorporates the effects of friction (`dq1`) and local losses (`dq5`). These terms typically resist flow, so a larger `dq1` or `dq5` will reduce the resulting `q`.

This equation is an explicit solution for `q` based on conditions from the previous time step and the previous iteration. The overall `dynwave.c` module will iterate, using this new `q` to update node heads, and then `dwflow_findConduitFlow` will be called again with updated heads until `q` converges for the current time step `dt`.

### 8.6 Numerical Scheme Aspects (Conceptual)

Several numerical techniques are employed to enhance stability and accuracy:

*   **Upstream Weighting (`rho` and `aWtd`, `rWtd`):**
    *   `rho = 1.0; if ( !isFull && qLast > 0.0 && h1 >= h2 ) rho = sigma;`
    *   `aWtd = a1 + (aMid - a1) * rho;`
    *   `rWtd = r1 + (rMid - r1) * rho;`
    *   The variables `aWtd` (weighted area) and `rWtd` (weighted hydraulic radius) are calculated using an upstream weighting factor `rho`. `rho` itself is influenced by the inertial damping factor `sigma`.
    *   When flow is positive (`qLast > 0.0`), not full, and the upstream head is greater than or equal to the downstream head (`h1 >= h2`), `rho` takes the value of `sigma`. Otherwise, `rho` is 1.0 (implying full upstream weighting for `aWtd` and `rWtd` if `sigma` is 1, or more central weighting if `sigma` is less than 1).
    *   This weighting is a common technique in finite difference schemes to ensure stability, particularly when dealing with advection-dominated flows or sharp gradients, by giving more influence to the upstream conditions which are "driving" the flow.

*   **Inertial Damping (`sigma`):**
    *   `if ( Link[j].froude <= 0.5 ) sigma = 1.0; else if ( Link[j].froude >= 1.0 ) sigma = 0.0; else sigma = 2.0 * (1.0 - Link[j].froude);`
    *   `sigma` is an inertial damping factor, ranging from 0.0 to 1.0. It reduces the influence of the inertial terms (`dq3` and `dq4`) under certain conditions.
    *   If the Froude number (`Link[j].froude`) is less than or equal to 0.5 (clearly subcritical), `sigma` is 1.0 (no damping).
    *   If the Froude number is greater than or equal to 1.0 (critical or supercritical), `sigma` is 0.0 (full damping of inertial terms).
    *   Between these Froude numbers, `sigma` transitions linearly.
    *   The user can also globally control inertial damping via the `InertDamping` option (`NO_DAMPING`, `PARTIAL_DAMPING` (which uses the Froude-based sigma), or `FULL_DAMPING` (sigma = 0)).
    *   Furthermore, if a closed conduit is flowing full (surcharged), `sigma` is forced to 0.0.
    *   Damping the inertial terms helps to stabilize the numerical solution, especially in rapidly changing flow conditions or near critical flow, where inertial effects can lead to oscillations or instability.

*   **Time Weighting and Iteration:**
    *   The use of `qOld` (flow from the previous time step) as the starting point in the numerator indicates that the scheme is stepping forward in time over the interval `dt`.
    *   The function `dwflow_findConduitFlow` is called iteratively within a larger loop in `dynwave.c` for each time step. In each iteration `steps`, the calculated `q` is refined using an under-relaxation parameter `omega`:
        `q = (1.0 - omega) * qLast + omega * q;`
        where `qLast` is the flow from the previous iteration. This is a common Picard iteration technique to improve convergence.
    *   The check `if ( q * qLast < 0.0 ) q = 0.001 * SGN(q);` prevents the flow from rapidly oscillating in direction during iterations, promoting stability.

### 8.7 Conduit Length (`length`)

*   The `length` variable used in the calculations is `Conduit[k].modLength`.
*   The comment in the code states: "--- use Courant-modified length instead of conduit's actual length".
*   The Courant-Friedrichs-Lewy (CFL) condition is a critical stability criterion for explicit numerical solutions of hyperbolic partial differential equations (like the St. Venant equations). It relates the time step, spatial discretization (length), and wave speed.
*   By adjusting the "effective length" (`modLength`), SWMM can internally manage stability for the chosen time step `dt` across a range of conduit lengths and flow conditions, helping to satisfy the CFL condition implicitly. This is a common technique in implementations of the St. Venant equations to allow for more flexible time stepping by the user.

In summary, `dwflow_findConduitFlow` implements a finite difference solution to the momentum equation for a single conduit, incorporating terms for friction, pressure gradient, inertia, local losses, and other physical losses. It employs numerical techniques like upstream weighting, inertial damping, and iterative refinement within a time step to achieve a stable and reasonably accurate solution for conduit flow as part of the overall dynamic wave routing algorithm.
