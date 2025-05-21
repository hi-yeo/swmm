# Python SWMM Algorithm Test Results: `calculate_conduit_trial_flow_momentum_eq`

This document presents the test results for the Python function `calculate_conduit_trial_flow_momentum_eq`, 
which translates core momentum equation logic from SWMM's `dwflow.c`.

## Test Execution Summary

### Test Case 1: Zero Initial Flow, Head Difference Driving Flow

#### Input Parameters:
| Parameter                        | Value          |
|----------------------------------|----------------|
| q_old_barrel                     | 0.0000         |
| q_last_iter_barrel               | 0.0000         |
| h1_node                          | 10.0000        |
| h2_node                          | 9.8000         |
| z1_conduit                       | 8.0000         |
| z2_conduit                       | 8.0000         |
| a_old_barrel                     | 1.9000         |
| length_modified                  | 100.0000       |
| length_actual                    | 100.0000       |
| conduit_roughness_factor         | 0.0025         |
| dt                               | 1.0000         |
| sigma_inertial_damping           | 0.5000         |
| rho_factor_upstream_weighting    | 0.7500         |
| has_local_losses                 | False          |
| coeff_loss_inlet                 | 0.0000         |
| coeff_loss_outlet                | 0.0000         |
| coeff_loss_avg                   | 0.0000         |
| GRAVITY                          | 32.2000        |
| FUDGE                            | 0.000001       |
| y1_current_iter                  | 2.0000         |
| y2_current_iter                  | 1.8000         |
| a1_current_iter                  | 2.0000         |
| a2_current_iter                  | 1.8000         |
| a_mid_current_iter               | 1.9000         |
| r1_current_iter                  | 0.4000         |
| r_mid_current_iter               | 0.3958         |

#### Comparison of Intermediate and Final Values:
| Variable                         | Expected Value  | Actual Value    | Abs Diff   | Rel Diff   | Status |
|----------------------------------|-----------------|-----------------|------------|------------|--------|
| v_current_iter_barrel            | 0.000000        | 0.000000        | 0.000000e+00 | 0.000000   | Pass   |
| aWtd                             | 1.925000        | 1.925000        | 0.000000e+00 | 0.000000   | Pass   |
| rWtd                             | 0.396875        | 0.396875        | 0.000000e+00 | 0.000000   | Pass   |
| dq1                              | 0.000000        | 0.000000        | 0.000000e+00 | 0.000000   | Pass   |
| dq2                              | -0.124030       | -0.124030       | 0.000000e+00 | 0.000000   | Pass   |
| dq3                              | 0.000000        | 0.000000        | 0.000000e+00 | 0.000000   | Pass   |
| dq4                              | 0.000000        | 0.000000        | 0.000000e+00 | 0.000000   | Pass   |
| dq5                              | 0.000000        | 0.000000        | 0.000000e+00 | 0.000000   | Pass   |
| dq6                              | 0.000000        | 0.000000        | 0.000000e+00 | 0.000000   | Pass   |
| denominator                      | 1.000000        | 1.000000        | 0.000000e+00 | 0.000000   | Pass   |
| q_trial_momentum_barrel          | 0.124030        | 0.124030        | 0.000000e+00 | 0.000000   | Pass   |

### Test Case 2: Existing Flow, Adverse Head, Friction Dominant

#### Input Parameters:
| Parameter                        | Value          |
|----------------------------------|----------------|
| q_old_barrel                     | 5.0000         |
| q_last_iter_barrel               | 5.0000         |
| h1_node                          | 10.0000        |
| h2_node                          | 10.1000        |
| z1_conduit                       | 8.0000         |
| z2_conduit                       | 8.0000         |
| a_old_barrel                     | 2.0500         |
| length_modified                  | 100.0000       |
| length_actual                    | 100.0000       |
| conduit_roughness_factor         | 0.0100         |
| dt                               | 1.0000         |
| sigma_inertial_damping           | 0.5000         |
| rho_factor_upstream_weighting    | 0.7500         |
| has_local_losses                 | False          |
| coeff_loss_inlet                 | 0.0000         |
| coeff_loss_outlet                | 0.0000         |
| coeff_loss_avg                   | 0.0000         |
| GRAVITY                          | 32.2000        |
| FUDGE                            | 0.000001       |
| y1_current_iter                  | 2.0000         |
| y2_current_iter                  | 2.1000         |
| a1_current_iter                  | 2.0000         |
| a2_current_iter                  | 2.1000         |
| a_mid_current_iter               | 2.0500         |
| r1_current_iter                  | 0.4000         |
| r_mid_current_iter               | 0.4020         |

#### Comparison of Intermediate and Final Values:
| Variable                         | Expected Value  | Actual Value    | Abs Diff   | Rel Diff   | Status |
|----------------------------------|-----------------|-----------------|------------|------------|--------|
| v_current_iter_barrel            | 2.439024        | 2.439024        | 1.929874e-07 | 0.000000   | Pass   |
| aWtd                             | 2.037500        | 2.037500        | 0.000000e+00 | 0.000000   | Pass   |
| rWtd                             | 0.401470        | 0.401471        | 5.882353e-07 | 0.000001   | Pass   |
| dq1                              | 0.084150        | 0.084152        | 1.674903e-06 | 0.000020   | Pass   |
| dq2                              | 0.065600        | 0.065604        | 3.750000e-06 | 0.000057   | Pass   |
| dq3                              | 0.000000        | 0.000000        | 0.000000e+00 | 0.000000   | Pass   |
| dq4                              | 0.002974        | 0.002974        | 2.302161e-08 | 0.000008   | Pass   |
| dq5                              | 0.000000        | 0.000000        | 0.000000e+00 | 0.000000   | Pass   |
| dq6                              | 0.000000        | 0.000000        | 0.000000e+00 | 0.000000   | Pass   |
| denominator                      | 1.084150        | 1.084152        | 1.674903e-06 | 0.000002   | Pass   |
| q_trial_momentum_barrel          | 4.554100        | 4.554093        | 7.110610e-06 | 0.000002   | Pass   |

### Test Case 3: Inertial Terms Significant, No Friction/Losses

#### Input Parameters:
| Parameter                        | Value          |
|----------------------------------|----------------|
| q_old_barrel                     | 2.0000         |
| q_last_iter_barrel               | 2.0000         |
| h1_node                          | 10.0000        |
| h2_node                          | 10.0000        |
| z1_conduit                       | 8.0000         |
| z2_conduit                       | 8.0000         |
| a_old_barrel                     | 1.0000         |
| length_modified                  | 100.0000       |
| length_actual                    | 100.0000       |
| conduit_roughness_factor         | 0.0000         |
| dt                               | 1.0000         |
| sigma_inertial_damping           | 1.0000         |
| rho_factor_upstream_weighting    | 0.7500         |
| has_local_losses                 | False          |
| coeff_loss_inlet                 | 0.0000         |
| coeff_loss_outlet                | 0.0000         |
| coeff_loss_avg                   | 0.0000         |
| GRAVITY                          | 32.2000        |
| FUDGE                            | 0.000001       |
| y1_current_iter                  | 2.0000         |
| y2_current_iter                  | 1.0000         |
| a1_current_iter                  | 2.0000         |
| a2_current_iter                  | 1.0000         |
| a_mid_current_iter               | 1.5000         |
| r1_current_iter                  | 0.4000         |
| r_mid_current_iter               | 0.3750         |

#### Comparison of Intermediate and Final Values:
| Variable                         | Expected Value  | Actual Value    | Abs Diff   | Rel Diff   | Status |
|----------------------------------|-----------------|-----------------|------------|------------|--------|
| v_current_iter_barrel            | 1.333333        | 1.333333        | 0.000000e+00 | 0.000000   | Pass   |
| aWtd                             | 1.625000        | 1.625000        | 0.000000e+00 | 0.000000   | Pass   |
| rWtd                             | 0.381250        | 0.381250        | 0.000000e+00 | 0.000000   | Pass   |
| dq1                              | 0.000000        | 0.000000        | 0.000000e+00 | 0.000000   | Pass   |
| dq2                              | 0.000000        | 0.000000        | 0.000000e+00 | 0.000000   | Pass   |
| dq3                              | 1.333333        | 1.333333        | 0.000000e+00 | 0.000000   | Pass   |
| dq4                              | -0.017778       | -0.017778       | 1.234568e-08 | 0.000001   | Pass   |
| dq5                              | 0.000000        | 0.000000        | 0.000000e+00 | 0.000000   | Pass   |
| dq6                              | 0.000000        | 0.000000        | 0.000000e+00 | 0.000000   | Pass   |
| denominator                      | 1.000000        | 1.000000        | 0.000000e+00 | 0.000000   | Pass   |
| q_trial_momentum_barrel          | 3.315550        | 3.315556        | 4.440892e-07 | 0.000000   | Pass   |
