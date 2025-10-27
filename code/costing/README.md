The Costing Model
================

- [1 Overview](#1-overview)
- [2 Parameters](#2-parameters)
- [3 Model details](#3-model-details)
- [4 Results](#4-results)
- [5 Attributions / Authors](#5-attributions--authors)

<!-- # Figures (temporary) {.unlisted .unnumbered} -->

This document describes the costing model that is used in the CEPI
application.

# 1 Overview

# 2 Parameters

| Math notation | Description | Distribution | Parameter 1 | Parameter 2 | Parameter 3 |
|:---|:--:|:--:|----|----|----|
| $D_{0; 365}$ | Preclinical trial duration (365); weeks | Constant | 14 |  |  |
| $D_{0; 200}$ | Preclinical trial duration (200DM); weeks | Constant | 5 |  |  |
| $D_{0; 100}$ | Preclinical trial duration (100DM); weeks | Constant | 5 |  |  |
| $D_{1; 365}$ | Phase I trial duration (365); weeks | Constant | 19 |  |  |
| $D_{1; 200}$ | Phase I trial duration (200DM); weeks | Constant | 7 |  |  |
| $D_{1; 100}$ | Phase I trial duration (100DM); weeks | Constant | 0 |  |  |
| $D_{2; 365}$ | Phase II trial duration (365); weeks | Constant | 19 |  |  |
| $D_{2; 200}$ | Phase II trial duration (200DM); weeks | Constant | 0 |  |  |
| $D_{2; 100}$ | Phase II trial duration (100DM); weeks | Constant | 0 |  |  |
| $D_{3; 365}$ | Phase III trial duration (365); weeks | Constant | 16 |  |  |
| $D_{3; 200}$ | Phase III trial duration (200DM); weeks | Constant | 15 |  |  |
| $D_{3; 100}$ | Phase III trial duration (100DM); weeks | Constant | 8 |  |  |
| $V_{L; 0}$ | Cost of vaccine delivery at start up (0–10%) in LIC; USD per dose | Triangular | 1 | 1.5 | 2 |
| $V_{L; 11}$ | Cost of vaccine delivery during ramp up (11–30%) in LIC; USD per dose | Triangular | 0.75 | 1 | 1.5 |
| $V_{L; 31}$ | Cost of vaccine delivery getting to scale (31–80%) in LIC; USD per dose | Triangular | 1 | 2 | 4 |
| $V_{LM; 0}$ | Cost of vaccine delivery at start up (0–10%) in LMIC; USD per dose | Triangular | 3 | 4.5 | 6 |
| $V_{LM; 11}$ | Cost of vaccine delivery during ramp up (11–30%) in LMIC; USD per dose | Triangular | 2.25 | 3 | 4.5 |
| $V_{LM; 31}$ | Cost of vaccine delivery getting to scale (31–80%) in LMIC; USD per dose | Triangular | 1.5 | 2 | 2.5 |
| $V_{UM; 0}$ | Cost of vaccine delivery at start up (0–10%) in UMIC; USD per dose | Triangular | 6 | 9 | 12 |
| $V_{UM; 11}$ | Cost of vaccine delivery during ramp up (11–30%) in UMIC; USD per dose | Triangular | 4.5 | 6 | 9 |
| $V_{UM; 31}$ | Cost of vaccine delivery getting to scale (31–80%) in UMIC; USD per dose | Triangular | 3 | 4 | 5 |
| $V_{H; 0}$ | Cost of vaccine delivery at start up (0–10%) in HIC; USD per dose | Triangular | 30 | 40 | 75 |
| $V_{H; 11}$ | Cost of vaccine delivery during ramp up (11–30%) in HIC; USD per dose | Triangular | 30 | 40 | 75 |
| $V_{H; 31}$ | Cost of vaccine delivery getting to scale (31–80%) in HIC; USD per dose | Triangular | 30 | 40 | 75 |
| $M_G$ | Global annual manufacturing volume; billion doses | Constant | 15 |  |  |
| $M_C$ | Current annual manufacturing volume; billion doses | Constant | 6 |  |  |
| $F$ | Facility transition start; weeks before vaccine approval | Constant | 7 |  |  |
| $I_R$ | Weeks to initial manufacturing, reserved infrastructure | Constant | 12 |  |  |
| $I_E$ | Weeks to initial manufacturing, existing and unreserved infrastructure | Constant | 30 |  |  |
| $I_B$ | Weeks to initial manufacturing, built and unreserved infrastructure | Constant | 48 |  |  |
| $C_R$ | Weeks to scale up to full capacity, reserved infrastructure | Constant | 10 |  |  |
| $C_E$ | Weeks to scale up to full capacity, existing and unreserved infrastructure | Constant | 16 |  |  |
| $C_B$ | Weeks to scale up to full capacity, built and unreserved infrastructure | Constant | 16 |  |  |
| $P_0$ | Probability of success; preclinical | Constant | 0.5 |  |  |
| $P_1$ | Probability of success; Phase I | Constant | 0.73 |  |  |
| $P_2$ | Probability of success; Phase II | Constant | 0.55 |  |  |
| $P_3$ | Probability of success; Phase III | Constant | 0.84 |  |  |
| $T_0$ | Trial cost, preclinical; million USD | Uniform | 1.7 | 140 |  |
| $T_1$ | Trial cost, Phase I; million USD | Uniform | 1.9 | 70 |  |
| $T_2$ | Trial cost, Phase II; million USD | Uniform | 3.8 | 140 |  |
| $T_3$ | Trial cost, Phase III; million USD | Uniform | 15 | 910 |  |
| $L$ | Licensure; USD | Constant | 287750 |  |  |
| $T_0^{(D)}$ | Trial duration, preclinical; years | Uniform | 1 | 2 |  |
| $T_1^{(D)}$ | Trial duration, Phase I; years | Uniform | 1 | 2 |  |
| $T_2^{(D)}$ | Trial duration, Phase II; years | Constant | 2 |  |  |
| $T_3^{(D)}$ | Trial duration, Phase III; years | Uniform | 2 | 4 |  |
| $L^{(D)}$ | Licensure duration; years | Constant | 2 |  |  |
| $B$ | BPSV cost of goods supplied; USD per dose | Constant | 4.68 |  |  |
| $A$ | Advanced capacity reservation fee; USD per dose per year | Constant | 0.53 |  |  |
| $S_R$ | SSV procurement price, reserved capacity; USD per dose | Constant | 6.29 |  |  |
| $S_U$ | SSV procurement price, reactive capacity; USD per dose | Constant | 18.94 |  |  |
| $E$ | Enabling activities; million USD per year | Constant | 700 |  |  |

Notation and parametric assumptions for inputs to the costing model

# 3 Model details

# 4 Results

# 5 Attributions / Authors
