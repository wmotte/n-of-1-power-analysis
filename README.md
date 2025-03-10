# N-of-1 Trial Power Simulator

## Overview

This R package simulates statistical power for multicycle N-of-1 trials using Bayesian linear models. It is specifically designed to help researchers and clinicians evaluate the feasibility of N-of-1 trials for rare disorders, with a particular focus on rare epilepsies and antiseizure medication (ASM) evaluation.

Unlike traditional group-based clinical trials, N-of-1 trials focus on individual treatment responses, making them ideal for rare disorders with small, heterogeneous patient populations. This simulator helps predict the likelihood of detecting clinically meaningful treatment effects under various conditions.

## Background

Clinical trials for rare disorders face significant challenges due to small sample sizes and heterogeneous patient populations. N-of-1 trials offer an alternative by focusing on individual treatment responses rather than group averages. 

For patients with rare epilepsies, the critical question is often not about group-level efficacy but whether a specific treatment reduces seizures to a clinically meaningful degree for an individual patient. This package helps answer this question by simulating trial outcomes across multiple scenarios.

## Features

- Simulate multicycle N-of-1 trials with alternating verum and placebo periods
- Model seizure count data using Bayesian linear models
- Assess model certainty in detecting treatment effects
- Evaluate the impact of various parameters:
  - Minimal clinically important difference (MCID)
  - Baseline seizure frequency
  - Treatment effect size
  - Placebo effect size
  - Number of cycles
  - Cycle duration
- Generate probability estimates for trial success
- Visualize certainty progression across cycles

## Installation

```r
# Install from GitHub
devtools::install_github("wmotte/n-of-1-power-simulator")
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `mcid` | Minimal clinically important difference (as proportion) | 0.5 |
| `baseline_freq` | Average seizure frequency during baseline | 7 |
| `treatment_effect` | Proportional reduction in seizures with treatment | 0.6 |
| `placebo_effect` | Proportional reduction in seizures with placebo | 0.2 |
| `cycles` | Number of cycles to simulate | 3 |
| `periods_per_cycle` | Number of periods per cycle (must be even) | 2 |
| `period_duration` | Duration of each period in days | 14 |

## Model Details

The package uses Bayesian linear models to analyze simulated seizure count data. For each simulation:

1. Seizure counts are generated based on specified parameters
2. A Bayesian linear model is fitted to the data after each cycle
3. The posterior distribution of the treatment effect is analyzed
4. The probability that the treatment effect exceeds the MCID is calculated

The model accounts for time trends, placebo effects, and random variation in seizure counts.

## Clinical Relevance

This simulator is particularly relevant for rare epilepsies where:

- Standard clinical trials are challenging due to small patient numbers
- Individual treatment responses are highly variable
- The goal is to determine if a treatment provides clinically meaningful seizure reduction for a specific patient
- Multiple treatment cycles may be needed to reach adequate certainty

For refractory epilepsy, where complete seizure freedom is often unlikely, detecting partial reductions (e.g., 30-50%) is clinically valuable. This tool helps researchers determine the feasibility of detecting such effects using an N-of-1 approach.

## Dependencies

- R (>= 4.0.0)
- brms
- dplyr
- tidyr
- ggplot2
- purrr
- rstan


