# Stepped-Wedge Power Analysis App

This repository contains an interactive Shiny application for conducting simulation-based power analyses for stepped-wedge randomized trials.

The app was designed to provide a flexible and user-friendly interface for exploring how different study design decisions influence statistical power in longitudinal stepped-wedge studies.

Users can modify core design assumptions, visualize the rollout structure of the intervention, and estimate statistical power under a range of implementation scenarios.

---

# Features

## Interactive stepped-wedge study design

The app allows users to customize major features of a stepped-wedge design, including:

- Number of rollout steps
- Number of sites assigned to each step
- Timing of rollout steps
- Length of intervention periods
- Study start and end dates
- Repeated follow-up measurement intervals

Both traditional stepped-wedge structures and clustered rollout structures are supported.


# Statistical approach

The application uses:

- Logistic regression (`glm`)
- Cluster-robust standard errors (`sandwich::vcovCL`)
- Site-level clustering adjustments
- Repeated longitudinal observations
- Simulation-based inference

Importantly, the app intentionally uses a standard generalized linear model (`glm`) with cluster-robust standard errors rather than a full mixed-effects model (`glmer`). 
This design choice substantially improves computational speed, making the application responsive enough for interactive use while still accounting for within-site dependence in the estimation of standard errors and statistical significance.

This approach provides a practical balance between:

- computational efficiency,
- realistic handling of clustered observations,
- and accessibility for rapid study planning and sensitivity analyses.
