The Costing Model
================

- [1 Parameters](#1-parameters)
- [2 Preparedness cost equation](#2-preparedness-cost-equation)
  - [2.1 BPSV advanced R&D](#21-bpsv-advanced-rd)
  - [2.2 BPSV investigational reserve](#22-bpsv-investigational-reserve)
  - [2.3 SSV capacity reservation](#23-ssv-capacity-reservation)
  - [2.4 Enabling activities](#24-enabling-activities)
- [3 Response cost equation](#3-response-cost-equation)
  - [3.1 Risk-adjusted R&D cost per candidate
    calculation](#31-risk-adjusted-rd-cost-per-candidate-calculation)
    - [3.1.1 SSV](#311-ssv)
    - [3.1.2 BPSV](#312-bpsv)
  - [3.2 Procurement cost calculation](#32-procurement-cost-calculation)
    - [3.2.1 SSV](#321-ssv)
    - [3.2.2 BPSV](#322-bpsv)
  - [3.3 Delivery Cost Equation](#33-delivery-cost-equation)
    - [3.3.1 SSV](#331-ssv)
    - [3.3.2 BPSV](#332-bpsv)
- [4 SSV delivery](#4-ssv-delivery)
  - [4.1 Timing](#41-timing)
  - [4.2 Production](#42-production)
  - [4.3 Allocation](#43-allocation)
  - [4.4 Delivery](#44-delivery)
- [5 BPSV delivery](#5-bpsv-delivery)
  - [5.1 Timing](#51-timing)
  - [5.2 Production](#52-production)
  - [5.3 Allocation](#53-allocation)
  - [5.4 Delivery](#54-delivery)
- [6 Parameter samples](#6-parameter-samples)
- [7 Contributors](#7-contributors)
- [8 References](#8-references)

<!-- # Figures (temporary) {.unlisted .unnumbered} -->

This document describes the costing model that is used in the CEPI
application.

# 1 Parameters

| Math notation | Code | Description | Distribution | Parameters | Source |
|:--:|:--:|:--:|:--:|:--:|:--:|
| $N^{\text{(BPSV)}}$ | n_bpsv_candidates | Number of BPSV candidates | Constant | 14 |  |
| $N^{\text{(BPSV-1)}}$ | n_bpsv_p1 | Number of BPSV candidates starting at phase 1 | Constant | 1 |  |
| $P_0^{\text{(BPSV)}}$ | pos_0 | Probability of success; preclinical | Multinomial | 0.40, 0.41, 0.41, 0.42, 0.48, 0.57 | Gouglas et al. (2018) |
| $P_1^{\text{(BPSV)}}$ | pos_1 | Probability of success; Phase I | Multinomial | 0.33, 0.40, 0.50, 0.68, 0.70, 0.72, 0.74, 0.77, 0.81, 0.90 | Gouglas et al. (2018) |
| $P_2^{\text{(BPSV)}}$ | pos_2 | Probability of success; Phase II | Multinomial | 0.22, 0.31, 0.33, 0.43, 0.46, 0.54, 0.58, 0.58, 0.74, 0.79 | Gouglas et al. (2018) |
| $P_3^{\text{(BPSV)}}$ | pos_3 | Probability of success; Phase III | Uniform | 0.4, 0.8 | Wong, Siah, and Lo (2019) |
| $Q^{\text{(SSV)}}$ | q_ssv_successes | Probability of N or more SSV successes | Constant | 0.9 | Model choice |
| $n^{\text{(SSV)}}$ | n_ssv_successes | Number of SSV successes | Constant | 5 | Model choice |
| $X_0$ | covid_ssv_0 | COVID-19 candidates failed at preclinical | Constant | 33 | Linksbridge SPC (2025) |
| $X_1$ | covid_ssv_1 | COVID-19 candidates failed at Phase 1 | Constant | 20 | Linksbridge SPC (2025) |
| $X_2$ | covid_ssv_2 | COVID-19 candidates failed at Phase 2 | Constant | 8 | Linksbridge SPC (2025) |
| $X_3$ | covid_ssv_3 | COVID-19 candidates failed at Phase 3 | Constant | 8 | Linksbridge SPC (2025) |
| $X_4$ | covid_ssv_4 | COVID-19 candidates successful | Constant | 27 | Linksbridge SPC (2025) |
| $Y_0^{(B)}$ | duration_0 | BPSV preclinical duration; years | Multinomial | 1, 2 | CEPI (2022) |
| $Y_1^{(B)}$ | duration_1 | BPSV Phase I duration; years | Multinomial | 1, 2 | CEPI (2022) |
| $Y_2^{(B)}$ | duration_2 | BPSV Phase II duration; years | Constant | 2 | CEPI (2022) |
| $Y_3^{(B)}$ | duration_3 | BPSV Phase III duration; years | Multinomial | 2, 3, 4 | CEPI (2022) |
| $W_3^{(B)}$ | duration_3_resp | BPSV response Phase III duration; weeks | Constant | 18 |  |
| $Y^{(200)}$ | years_200 | Years of R&D to 200-day readiness | Constant | 5 |  |
| $Y^{(100)}$ | years_100 | Years of R&D to 100-day readiness | Constant | 15 |  |
| $E$ | cost_enab | Enabling activities; million USD per year | Constant | 700 | CEPI (2021) |
| $W_{0; 365}^{(S)}$ | weeks_P0_365 | SSV preclinical duration (365); weeks | Constant | 14 |  |
| $W_{0; 200}^{(S)}$ | weeks_P0_200 | SSV preclinical duration (200 days); weeks | Constant | 5 |  |
| $W_{0; 100}^{(S)}$ | weeks_P0_100 | SSV preclinical duration (100 days); weeks | Constant | 5 |  |
| $W_{1; 365}^{(S)}$ | weeks_P1_365 | SSV phase I duration (365); weeks | Constant | 0 |  |
| $W_{1; 200}^{(S)}$ | weeks_P1_200 | SSV phase I duration (200 days); weeks | Constant | 0 |  |
| $W_{1; 100}^{(S)}$ | weeks_P1_100 | SSV phase I duration (100 days); weeks | Constant | 0 |  |
| $W_{2; 365}^{(S)}$ | weeks_P2_365 | SSV phase II duration (365); weeks | Constant | 19 |  |
| $W_{2; 200}^{(S)}$ | weeks_P2_200 | SSV phase II duration (200 days); weeks | Constant | 7 |  |
| $W_{2; 100}^{(S)}$ | weeks_P2_100 | SSV phase II duration (100 days); weeks | Constant | 0 |  |
| $W_{3; 365}^{(S)}$ | weeks_P3_365 | SSV phase III duration (365); weeks | Constant | 16 |  |
| $W_{3; 200}^{(S)}$ | weeks_P3_200 | SSV phase III duration (200 days); weeks | Constant | 15 |  |
| $W_{3; 100}^{(S)}$ | weeks_P3_100 | SSV phase III duration (100 days); weeks | Constant | 8 |  |
| $T_0^{(e)}$ | cost_0_ex | Cost, preclinical, experienced manufacturer; USD | Exponential | 24213683 | Gouglas et al. (2018) |
| $T_0^{(n)}$ | cost_0_inex | Cost, preclinical, inexperienced manufacturer; USD | Inverse Gaussian | 7882792, 13455907 | Gouglas et al. (2018) |
| $T_1^{(e)}$ | cost_1_ex | Cost, Phase I, experienced manufacturer; USD | Inverse Gaussian | 15339198, 8076755 | Gouglas et al. (2018) |
| $T_1^{(n)}$ | cost_1_inex | Cost, Phase I, inexperienced manufacturer; USD | PearsonV | 2.2774, 9799081 | Gouglas et al. (2018) |
| $T_2^{(e)}$ | cost_2_ex | Cost, Phase II, experienced manufacturer; USD | Log normal | 28297339, 24061641 | Gouglas et al. (2018) |
| $T_2^{(n)}$ | cost_2_inex | Cost, Phase II, inexperienced manufacturer; USD | Inverse Gaussian | 17124622, 35918793 | Gouglas et al. (2018) |
| $T_3^{(e)}$ | cost_3_ex | Cost, Phase III, experienced manufacturer; USD | PearsonV | 1.3147, 51397313 | Gouglas et al. (2018) |
| $T_3^{(n)}$ | cost_3_inex | Cost, Phase III, inexperienced manufacturer; USD | PearsonVI | 4.8928, 1.6933, 11400026 | Gouglas et al. (2018) |
| $T_4$ | cost_lic | Licensure cost, 2018; USD | Constant | 287750 | Gouglas et al. (2018) |
| $I$ | inflation | Inflation (2018 t0 2025) | Constant | 0.28 | U.S. BLS (2025) |
| $\omega$ | inex_weight | Share of manufacturers that are inexperienced | Constant | 0.875 | See Table <a href="#tab:inex">2.1</a> |
| $A_4$ | bpsv_inv_res | Size of BPSV investigational reserve, doses | Constant | 100000 | Model choice |
| $A_1$ | cost_bpsv_res | Annual BPSV reservation cost, USD per dose | Constant | 0.0101213633333333 |  |
| $A_5$ | bpsv_res_upfront | BPSV reserve upfront cost, USD per dose | Constant | 0.115 |  |
| $Y_{rep}$ | bpsv_replenishment | Years after which BPSV doses are to be replaced | Constant | 3 |  |
| $A_2$ | cost_capres | Advanced capacity reservation fee; USD per dose per year | Constant | 0.531692307692308 | Pfizer (2023) |
| $A_3$ | hic_cap_res | Reserved capacity for HIC, billions | Constant | 0.5 |  |
| $S_U$ | cost_un | SSV procurement price, reactive capacity; USD per dose | Constant | 18.9392 | Linksbridge SPC (2025) |
| $G$ | cost_cogs | Drug substance cost; USD per dose | Constant | 4.68 | Kazaz (2021) |
| $M_p$ | profit | Profit margin | Constant | 0.2 | Kazaz (2021) |
| $M_f$ | cost_ff | Fill/finish cost | Constant | 0.1398 | Kazaz (2021) |
| $M_t$ | cost_travel | Cost to transport product | Constant | 0.12 | Kazaz (2021) |
| $M_G$ | man_glo | Global annual manufacturing volume; billion doses | Constant | 15 | Linksbridge SPC (2025) |
| $M_C$ | man_curr | Current annual manufacturing volume; billion doses | Constant | 9 | Linksbridge SPC (2025) |
| $\lambda$ | final_vaccine_coverage | Final vaccine coverage, proportion of population | Constant | 0.8 | Model choice |
| $\delta$ | vaccine_wastage | Fraction of BPSV expected to go to waste | Constant | 0.3142532 | Model choice |
| $N^{\text{(boost)}}$ | n_boosters | Number of boosters given, one per year | Constant | 2 | Model choice |
| $I_0$ | week_trans_start | Facility transition start; weeks before vaccine approval | Constant | 7 |  |
| $I_R$ | weeks_init_res | Weeks to initial manufacturing, reserved infrastructure | Constant | 12 | Vaccines Europe (2023) |
| $I_{E,0}$ | weeks_init_ex_nb | Weeks to initial manufacturing when there’s no BPSV, existing and unreserved infrastructure | Constant | 30 | Vaccines Europe (2023) |
| $I_{E,1}$ | weeks_init_ex_bp | Weeks to initial manufacturing when there’s BPSV, existing and unreserved infrastructure | Constant | 12 | Vaccines Europe (2023) |
| $I_B$ | weeks_init_bui | Weeks to initial manufacturing, built and unreserved infrastructure | Constant | 48 |  |
| $C_R$ | weeks_scale_res | Weeks to scale up to full capacity, reserved infrastructure | Constant | 10 | Vaccines Europe (2023) |
| $C_E$ | weeks_scale_ex | Weeks to scale up to full capacity, existing and unreserved infrastructure | Constant | 16 |  |
| $C_B$ | weeks_scale_bui | Weeks to scale up to full capacity, built and unreserved infrastructure | Constant | 16 |  |
| $V_{L; 0}$ | cost_lic_0 | Cost of vaccine delivery at start up (0–10%) in LIC; USD per dose | Triangular | 1, 1.5, 2 | See Table <a href="#tab:delcosts">3.1</a> |
| $V_{L; 11}$ | cost_lic_11 | Cost of vaccine delivery during ramp up (11–30%) in LIC; USD per dose | Triangular | 0.75, 1, 1.5 | See Table <a href="#tab:delcosts">3.1</a> |
| $V_{L; 31}$ | cost_lic_31 | Cost of vaccine delivery getting to scale (31% and over) in LIC; USD per dose | Triangular | 1, 2, 4 | See Table <a href="#tab:delcosts">3.1</a> |
| $V_{LM; 0}$ | cost_lmic_0 | Cost of vaccine delivery at start up (0–10%) in LMIC; USD per dose | Triangular | 3, 4.5, 6 | See Table <a href="#tab:delcosts">3.1</a> |
| $V_{LM; 11}$ | cost_lmic_11 | Cost of vaccine delivery during ramp up (11–30%) in LMIC; USD per dose | Triangular | 2.25, 3, 4.5 | See Table <a href="#tab:delcosts">3.1</a> |
| $V_{LM; 31}$ | cost_lmic_31 | Cost of vaccine delivery getting to scale (31% and over) in LMIC; USD per dose | Triangular | 1.5, 2, 2.5 | See Table <a href="#tab:delcosts">3.1</a> |
| $V_{UM; 0}$ | cost_umic_0 | Cost of vaccine delivery at start up (0–10%) in UMIC; USD per dose | Triangular | 6, 9, 12 | See Table <a href="#tab:delcosts">3.1</a> |
| $V_{UM; 11}$ | cost_umic_11 | Cost of vaccine delivery during ramp up (11–30%) in UMIC; USD per dose | Triangular | 4.5, 6, 9 | See Table <a href="#tab:delcosts">3.1</a> |
| $V_{UM; 31}$ | cost_umic_31 | Cost of vaccine delivery getting to scale (31% and over) in UMIC; USD per dose | Triangular | 3, 4, 5 | See Table <a href="#tab:delcosts">3.1</a> |
| $V_{H; 0}$ | cost_hic_0 | Cost of vaccine delivery at start up (0–10%) in HIC; USD per dose | Triangular | 30, 40, 75 | See Table <a href="#tab:delcosts">3.1</a> |
| $V_{H; 11}$ | cost_hic_11 | Cost of vaccine delivery during ramp up (11–30%) in HIC; USD per dose | Triangular | 30, 40, 75 | See Table <a href="#tab:delcosts">3.1</a> |
| $V_{H; 31}$ | cost_hic_31 | Cost of vaccine delivery getting to scale (31% and over) in HIC; USD per dose | Triangular | 30, 40, 75 | See Table <a href="#tab:delcosts">3.1</a> |
| $r$ | discount | Discount rate | Uniform | 0.02, 0.06 | Glennerster, Snyder, and Tan (2023) |
| $N_{HIC}^{(0)}$ | pop_hic_0 | Population, HIC | Constant | 1260028362 |  |
| $N_{UMIC}^{(0)}$ | pop_umic_0 | Population, UMIC | Constant | 2854556263.5 |  |
| $N_{LMIC}^{(0)}$ | pop_lmic_0 | Population, LMIC | Constant | 3314048516 |  |
| $N_{LIC}^{(0)}$ | pop_lic_0 | Population, LIC | Constant | 762656294\.5 |  |
| $N_{HIC}^{(15)}$ | pop_hic_15 | Population aged 15 and older, HIC | Constant | 1064531991.5 |  |
| $N_{UMIC}^{(15)}$ | pop_umic_15 | Population aged 15 and older, UMIC | Constant | 2308984518 |  |
| $N_{LMIC}^{(15)}$ | pop_lmic_15 | Population aged 15 and older, LMIC | Constant | 2363976954.5 |  |
| $N_{LIC}^{(15)}$ | pop_lic_15 | Population aged 15 and older, LIC | Constant | 450976596\.5 |  |
| $N_{HIC}^{(65)}$ | pop_hic_65 | Population aged 65 and older, HIC | Constant | 256715334 |  |
| $N_{UMIC}^{(65)}$ | pop_umic_65 | Population aged 65 and older, UMIC | Constant | 359824402\.5 |  |
| $N_{LMIC}^{(65)}$ | pop_lmic_65 | Population aged 65 and older, LMIC | Constant | 215830985\.5 |  |
| $N_{LIC}^{(65)}$ | pop_lic_65 | Population aged 65 and older, LIC | Constant | 24812768 |  |

Notation and parametric assumptions for inputs to the costing model.
Parameters are used as follows: uniform distributions go from Parameter
1 to Parameter 2. Triangular distributions go from Parameter 1 to
Parameter 3 with a peak at Parameter 2. Multinomial distributions have
equally probable values listed individually. Exponential distributions
have as a mean Parameter 1. Inverse Gaussian distributions have as a
mean Parameter 1, and as a shape Parameter 2. Log normal distributions
have as a mean Parameter 1, and as a standard deviation Parameter 2.
PearsonV distributions have shape Parameter 1, scale Parameter 2, and
location 0. PearsonVI distributions have shape Parameters 1 and 2, scale
Parameter 3, and location 0. Where given, distributions are truncated at
bounds.

# 2 Preparedness cost equation

<span style="color:red;">(BPSV R&D + BPSV Stockpile + SARS-X Reserved
capacity + Enabling activities) / (1 + discount rate) ^ (year –
2025)</span>

$$D_y^{\text{(prep)}} = \frac{1}{(1+r)^y}\left(D_s^{\text{(BP-adRD)}} + D_{s,y}^{\text{(BP-man)}} + D_{s,y}^{\text{(BP-inv)}} + D_s^{\text{(S-cap)}} + D_{s,y}^{\\text{(en)}}\right)$$

- $D_s^{\text{(BP-adRD)}}$ is the R&D cost of BPSV prior to an outbreak;
  see Equation (2.1)
- $D_{s,y}^{\text{(BP-man)}}$ is the upfront cost of maintaining an
  investigational reserve of 100,000 BPSV doses; see Equation (2.2)
- $D_{s,y}^{\text{(BP-inv)}}$ is the annual cost of maintaining an
  investigational reserve of 100,000 BPSV doses; see Equation (2.3)
- $D_s^{\text{(S-cap)}}$ is the cost of reserved capacity for SSV; see
  Equation (2.4)
- $D_{s,y}^{\\text{(en)}}$ is the annual cost of enabling activities;
  see Equation (2.5).

## 2.1 BPSV advanced R&D

**These values match the spreadsheet results**

| Developer | Licensure Experience |
|:----------|:---------------------|
| CalTech   | No                   |
| SK Bio    | Yes                  |
| Codiak    | No                   |
| Panacea   | No                   |
| NEC Onco  | No                   |
| Intravacc | No                   |
| VIDO      | No                   |
| IVI       | No                   |

<span id="tab:inex"></span>Table 2.1: Manufacturers working on BPSV and
whether or not they have licensure experience

Probabilities of success for preclinical, Phase I, Phase II, and Phase
III are $P_0$, $P_1$, $P_2$ and $P_3$. Then probabilities of occurrence
are:

``` math
\hat{P}_i^{(0)} = \left\{\begin{array}{lr}1 & i=0 \\ 
\prod_{j=0}^{i-1}P_j & i>0 
\end{array}\right.
```

for $i \in \{0,1,2,3,4\}$, with $i=4$ corresponding to licensure.

For $N^{\text{(BPSV-1)}} = 1$ candidate(s), which have already been
through the preclinical phase, we have

``` math
\hat{P}_i^{(1)} = \left\{\begin{array}{lr}1 & i=1 \\ 
\prod_{j=1}^{i-1}P_j & i>1 
\end{array}\right.
```

for $i \in \{1,2,3,4\}$.

The cost of each phase is $T_i$, a weighted average of experienced and
inexperienced manufacturers (with $\omega = 0.875$):

$$T_{i} = (1+I)(\omega T_i^{(n)} + (1-\omega)T_i^{(e)}).$$

where $I = 0.28$ is inflation from 2018 to 2025. Then the total weighted
cost for phases 0 through 2 for $N^{\text{(BPSV)}}$ candidates is

$$\begin{equation}
D_s^{\text{(BP-adRD)}} = \left\\{\begin{array}{lr}
 \left(N^{\text{(BPSV)}}-N^{\text{(BPSV-1)}}\right)\sum_{i=0}^2 \hat{P}_i^{(0)}T_{i} + N^{\text{(BPSV-1)}}\sum_{i=1}^2 \hat{P}_i^{(1)}T_{i} \\; & \\; s=1 \\\\
0  \\; & \\; s\neq 1
\end{array}\right.
\qquad(2.1)
\end{equation}$$

<div class="figure">

<img src="README_files/figure-gfm/posbpsv-1.png" alt="Risk-adjusted R&amp;D cost for 8 BPSV candidates"  />
<p class="caption">
<span id="fig:posbpsv"></span>Figure 2.1: Risk-adjusted R&D cost for 8
BPSV candidates
</p>

</div>

Min. 1st Qu. Median Mean 3rd Qu. Max. 22.24 58.29 69.77 68.45 85.04
105.38

Target: 146 (103 135 177)

## 2.2 BPSV investigational reserve

**The annual cost annual is correct, at around 162 thousand, but the
total cost is slightly too high**

The time taken to complete development of the BPSV up to the end of
phase II, from which point it is manufactured to be held in an
investigational reserve, is:

$$Y^{(B)} = Y_0^{(B)} + Y_1^{(B)} + Y_2^{(B)}.$$

The upfront cost of securing the investigational reserve is

$$\begin{equation}
D_{s,y}^{\text{(BP-man)}} = \left\\{\begin{array}{lr}
 A_4A_5
\\; & \\; s=1 \\;\\&\\;y=Y^{(B)}+1\\\\
0  \\; & \\; s\neq 1\\;\\|\\;y\neq Y^{(B)}+1
\end{array}\right.
\qquad(2.2)
\end{equation}$$

where $A_4 =100,000$ is the size of the reserve and $A_5 =0.115$ is the
cost per dose in USD.

The cost of goods supplied is $G = 4.68$. Then the cost of drug
substance is $G(1-M_f)(1+M_p) = 4.83$ USD per dose. The reserve is
replenished every $Y_{rep} = 3$ years. Then the annual cost to maintain
the reserve of $A_4 =100,000$ doses is

$$\begin{equation}
D_{s,y}^{\text{(BP-inv)}} = \left\\{\begin{array}{lr}
 \frac{A_4}{Y_{rep}}G  (1-M_f)(1+M_p) + A_1
\\; & \\; s=1 \\;\\&\\;y>Y^{(B)}\\\\
0  \\; & \\; s\neq 1\\;\\|\\;y\leq Y^{(B)}
\end{array}\right.
\qquad(2.3)
\end{equation}$$

where $A_1 = 0.01$ USD is the annual reservation cost per dose.

<div class="figure">

<img src="README_files/figure-gfm/bpsvinv-1.png" alt="BPSV investigational reserve costs accumulated from the completion of Phase II to year 15 with uniformly distributed discount rate."  />
<p class="caption">
<span id="fig:bpsvinv"></span>Figure 2.2: BPSV investigational reserve
costs accumulated from the completion of Phase II to year 15 with
uniformly distributed discount rate.
</p>

</div>

Min. 1st Qu. Median Mean 3rd Qu. Max. 1067 1171 1215 1240 1287 1458

Target: 1 (0.9 1 1.1)

## 2.3 SSV capacity reservation

**This matches the spreadsheet results.**

The cost per dose reservation per year is $A_2 = 0.53$ USD. Reservation
sizes, in billions, depend on scenarios, including the $A_3 = 0.5$
billion doses reserved for HIC, as follows:

$$\begin{equation}
M_{R,s} = \left\\{\begin{array}{lr}A_3 & s\in\\{0, 1, 4, 7, 10\\} \\\\ 
A_3+0.7 & s\in\\{2, 5, 8\\} \\\\ 
A_3+2 & s\in\\{3, 6, 9\\} \end{array}\right.
\end{equation}$$

Then the total cost per year is

$$\begin{equation}
D_s^{\text{(S-cap)}} =  M_{R,s} A_2
\qquad(2.4)
\end{equation}$$

The annual costs in billion USD are 0.27, 0.64, and 1.33, respectively.

<div class="figure">

<img src="README_files/figure-gfm/capres-1.png" alt="Capacity reservation costs accumulated over 15 years with uniformly distributed discount rate."  />
<p class="caption">
<span id="fig:capres"></span>Figure 2.3: Capacity reservation costs
accumulated over 15 years with uniformly distributed discount rate.
</p>

</div>

0 Min. 1st Qu. Median Mean 3rd Qu. Max. 3034 3233 3320 3301 3431 3472

0.7 Min. 1st Qu. Median Mean 3rd Qu. Max. 7283 7758 7969 7923 8235 8334

2 Min. 1st Qu. Median Mean 3rd Qu. Max. 15172 16163 16601 16506 17156
17362

Targets: 3,086 (2,897 3,074 3,269)

7,407 (6,954 7,378 7,845)

15,431 (14,487 15,370 16,344)

## 2.4 Enabling activities

**This matches the spreadsheet results**

Denote the “Days Mission” by $\zeta$, so that
$\zeta\in\lbrace 365, 200, 100 \rbrace$. Then annual costs, $E=700$
million, accumulate depending on the year and the mission:

$$\begin{equation}
D_{s,y}^{\text{(en)}} = \left\\{\begin{array}{lr}E & \zeta(s)=200 \\;\\&\\; y\leq 5 \\; |\\; \zeta(s)=100\\; \\& \\;y\leq 15 \\\\ 
0 & \zeta(s)=365 \\;|\\; y > 15 \\;|\\; \zeta(s)=200 \\;\\&\\; y \\;>\\; 5  \end{array}\right.
\qquad(2.5)
\end{equation}$$

For our scenarios, we have

$$\begin{equation}
\zeta(s) = \left\\{\begin{array}{lr} 365 & s\in\\{0, 1, 2, 3, 10\\} \\\\ 
200 & s\in\\{4, 5, 6\\} \\\\ 
100 & s\in\\{7, 8, 9\\} \end{array}\right.
\end{equation}$$

<div class="figure">

<img src="README_files/figure-gfm/en-1.png" alt="Enabling costs accumulated over 15 years with uniformly distributed discount rate."  />
<p class="caption">
<span id="fig:en"></span>Figure 2.4: Enabling costs accumulated over 15
years with uniformly distributed discount rate.
</p>

</div>

365 Min. 1st Qu. Median Mean 3rd Qu. Max. 0 0 0 0 0 0

200 Min. 1st Qu. Median Mean 3rd Qu. Max. 3.23 3.29 3.32 3.31 3.35 3.36

100 Min. 1st Qu. Median Mean 3rd Qu. Max. 7.99 8.51 8.74 8.69 9.03 9.14

Targets:

3,242 (3,182 3,241 3,302)

8,126 (7,629 8,094 8,607)

# 3 Response cost equation

<span style="color:red;">(BPSV R&D + SARS-X R&D + BPSV Procurement +
SARS-X Procurement + BPSV Delivery + SARS-X Delivery) / (1 + discount
rate) ^ (year – 2025)</span>

$$D_{s,y}^{\text{(res)}} = \frac{1}{(1+r)^y}\left(D_{s,y}^{\text{(BP-resRD)}} + D_{s,y}^{\text{(S-RD)}} + D_{s,y}^{\text{(BP-proc)}} + D_{s,y}^{\text{(S-proc)}} + D_{s,y}^{\text{(BP-del)}} + D_{s,y}^{\text{(S-del)}}\right)$$

- $D_{s,y}^{\text{(BP-resRD)}}$ is the R&D cost of BPSV after an
  outbreak; see Equation (3.2)
- $D_{s,y}^{\text{(S-RD)}}$ is the R&D cost for SSV; see Equation (3.1)
- $D_{s,y}^{\text{(BP-proc)}}$ is the cost of procuring BPSV; see
  Equation (3.4)
- $D_{s,y}^{\text{(S-proc)}}$ is the cost of procuring SSV; see Equation
  (3.3)
- $D_{s,y}^{\text{(BP-del)}}$ is the cost of delivering BPSV; see
  Equation (3.6)
- $D_{s,y}^{\text{(S-del)}}$ is the cost of delivering SSV; see Equation
  (3.5)

## 3.1 Risk-adjusted R&D cost per candidate calculation

### 3.1.1 SSV

**These don’t match the spreadsheet results. Values too low.**

Trial costs are adjusted for the duration of the trial, which depend on
the R&D investment, denoted $\zeta\in\lbrace 365, 200, 100\rbrace$:

$$T_{\zeta,i}^{(e)} = (1+I)\frac{W_{i;\zeta}^{(S)}}{52Y_{i}^{(B)}}T_i^{(e)}.$$

where $I = 0.28$ is inflation from 2018 to 2025.

The probability of success of each phase comes from COVID-19 data.

$$\begin{equation}
P_i^{\text{(SSV)}} \sim \text{Beta}\left(\sum_{j=i+1}^4X_j+1, X_i + 1 \right)
\end{equation}$$

for $i \in \lbrace 0,1,2,3 \rbrace$ where $X_i$ is the number of
candidates that failed in phase $i$ and $X_4$ the number that succeeded.

Then the total cost is

$$\begin{equation}
D_s^{\text{(S-RD)}} = N^{\text{(SSV)}}\left(\sum_{i=0}^3 \hat{P}_i^{\text{(SSV)}} \cdot T_{\zeta(s),i}^{(e)} +  \hat{P}_4^{\text{(SSV)}} T_4\right)
\qquad(3.1)
\end{equation}$$

We multiply by the number of candidates, $N^{\text{(SSV)}}$, to get the
total cost from the weighted average per candidate, where

$$\begin{equation}
N^{\text{(SSV)}} = n^{\text{(SSV)}} + F^{-1}_{NegBin}\left(Q^{\text{(SSV)}}; n^{\text{(SSV)}},  \hat{P}_3^{\text{(SSV)}} \right)
\end{equation}$$

is chosen to secure at least $n^{\text{(SSV)}} = 5$ successful
candidates with probability $Q^{\text{(SSV)}} = 90$%. Here,
$F^{-1}_{NegBin}\left(q; n,  p \right)$ is the cumulative density of a
negative binomial distribution with parameters $n$ and $p$ evaluated at
quantile $q$.

<div class="figure">

<img src="README_files/figure-gfm/posssv-1.png" alt="Risk-adjusted R&amp;D cost for 18 SSV candidates"  />
<p class="caption">
<span id="fig:posssv"></span>Figure 3.1: Risk-adjusted R&D cost for 18
SSV candidates
</p>

</div>

| DM  | Min. | 1st Qu. | Median | Mean | 3rd Qu. | Max. |
|:---:|:----:|:-------:|:------:|:----:|:-------:|:----:|
| 365 |  26  |   38    |   46   |  45  |   52    |  61  |
| 200 |  16  |   22    |   28   |  27  |   30    |  39  |
| 100 |  10  |   15    |   17   |  17  |   19    |  24  |

Targets:

284 (105 170 283)

195 (61 97 164)

118 (35 61 108)

### 3.1.2 BPSV

**This is a little higher than the spreadsheet results**

**I have basically assumed the same as SSV except for the numbers given
(8 candidates and 18 weeks)**

The BPSV has $N^{\text{(BPSV)}}$ candidates. Those that have passed
through Phases 0 to 2 prior to the outbreak go through Phase 3 during
the response. The duration is $W_3^{(B)}=18$ weeks. Thus we write the
BPSV R&D response cost

$$\begin{equation}
D_s^{\text{(BP-resRD)}} = \left\\{\begin{array}{lr}N^{\text{(BPSV)}}\hat{P}_3\left( (1+I)\frac{W_3^{(B)}}{52Y_3^{(B)}}T_3^{(e)} + P_3 T_4\right) \\; & \\; s=1 \\\\
0  \\; & \\; s\neq 1
\end{array}\right.
\qquad(3.2)
\end{equation}$$

<div class="figure">

<img src="README_files/figure-gfm/bpsvresrd-1.png" alt="Reactive R&amp;D cost for BPSV"  />
<p class="caption">
<span id="fig:bpsvresrd"></span>Figure 3.2: Reactive R&D cost for BPSV
</p>

</div>

Min. 1st Qu. Median Mean 3rd Qu. Max. 0.8 1.0 2.9 4.1 6.7 11.9

Target: 14 (3 5 10)

## 3.2 Procurement cost calculation

The cost per dose comes from the cost of goods supplied ($G = 4.68$)
adjusted for profits ($M_p = 0.2$) and the transportation cost
($M_t = 0.12$).

$S_R = G(1+M_p)(1+M_t)$ evaluates to 6.29.

This cost is used both for SSV doses manufactured using reserved
capacity, and all newly manufactured BPSV doses.

### 3.2.1 SSV

**These values are close, but not identical, to the spreadsheet results
if I adjust for the total demand**

We write billion doses procured from channel $x\in\lbrace R,E,B \rbrace$
in year $y$ and scenario $s$ as $A_{x,s,y}$ (see Equation (4.1)). Then
the total cost, in billion USD, is:

$$\begin{equation}
D_{s,y}^{\text{(S-proc)}} = A_{R,s,y} S_R  + \sum_{x\in\lbrace E,B \rbrace}  A_{x,s,y} S_U
\qquad(3.3)
\end{equation}$$

Here, $S_R = 6.29$ is the cost per reserved dose and $S_U = 18.94$ the
cost per unreserved dose in USD.

<div class="figure">

<img src="README_files/figure-gfm/costperyear-1.png" alt="SSV procurement cost"  />
<p class="caption">
<span id="fig:costperyear"></span>Figure 3.3: SSV procurement cost
</p>

</div>

| Scenario |  Min.  | 1st Qu. | Median |  Mean  | 3rd Qu. |  Max.  |
|:--------:|:------:|:-------:|:------:|:------:|:-------:|:------:|
|   BAU    | 177626 | 213237  | 230030 | 227100 | 252152  | 260578 |
|   S01    | 180975 | 216598  | 233358 | 230410 | 255409  | 263802 |
|   S02    | 164357 | 197208  | 212693 | 209988 | 233088  | 240856 |
|   S03    | 133031 | 159380  | 171786 | 169609 | 188114  | 194330 |
|   S04    | 183579 | 219158  | 235866 | 232906 | 257822  | 266173 |
|   S05    | 166819 | 199088  | 214238 | 211552 | 234145  | 241716 |
|   S06    | 139015 | 165783  | 178343 | 176111 | 194842  | 201115 |
|   S07    | 183695 | 218990  | 235550 | 232605 | 257298  | 265567 |
|   S08    | 165625 | 197361  | 212246 | 209596 | 231790  | 239221 |
|   S09    | 135390 | 161160  | 173238 | 171081 | 189089  | 195113 |
|   S10    | 175500 | 210673  | 227258 | 224365 | 249107  | 257429 |

Costs summed and discounted from year 16 to year 20, million USD

Targets:

184,127 ( 151,271 180,171 214,966 ) 187,255 ( 154,376 183,358 218,147 )
167,519 ( 137,713 163,938 195,495 ) 135,910 ( 111,925 133,050 158,444 )
189,820 ( 157,000 185,976 220,684 ) 169,549 ( 140,293 166,133 197,067 )
141,440 ( 117,134 138,613 164,309 ) 189,878 ( 157,295 186,091 220,526 )
168,378 ( 139,564 165,035 195,494 ) 137,984 ( 114,513 135,278 160,078 )
178,766 ( 146,883 174,927 208,686 )

### 3.2.2 BPSV

**This is pretty close**

$$\begin{equation}
D_s^{\text{(BP-proc)}} = \left\\{\begin{array}{lr}
A_{BPSV,s}\cdot S_R +  A_4(M_f+M_t)(1+M_p)G\\; & \\; s=1 \\\\
0  \\; & \\; s\neq 1
\end{array}\right.
\qquad(3.4)
\end{equation}$$

For a world population aged 65 and over of 0.9 billion, an uptake of 80%
(accounting for wastage of 31%), and a cost per dose of $S_R = 6.29$ USD
(the same as for SSV via reserved capacity), the procurement cost for
BPSV is 6.68 billion USD.

Although 1.0625 billion doses are manufactured, as manufacturing stops
once one billion doses have been made.

Min. 1st Qu. Median Mean 3rd Qu. Max. 3597 4180 4448 4397 4796 4927

Target: 3,628 (3,062 3,568 4,165)

## 3.3 Delivery Cost Equation

### 3.3.1 SSV

**These values are ballpark correct but too concentrated**

For populations aged 15 and above $N_i^{(15)}$ in income group
$i\in\lbrace\text{LIC, LMIC, UMIC, HIC}\rbrace$, we write

$$L_i = 2\cdot \lambda\cdot N_i^{(15)}$$

as the total demand for first-schedule doses for income group $i$,
representing two doses each for $\lambda = 80$% of the population.

We write the delivery cost for $h_{s,i,w}$ doses given in week $w$ and
country type $i$ as follows. There are three cost tiers, the first of
which is applied to the first 10% of $L_i$, the second to the subsequent
20%, and the third to all doses thereafter. The same costing schedule
applies both to the first-schedule plus booster SSV doses and the BPSV
rollout.

$$\begin{equation}
H_{s,i,w} = 
\left\\{\begin{array}{lr}
V_{i; 0}h_{s,i,w}  & \sum_{j=1}^{w-1}h_{s,i,w}\leq \frac{1}{10}L_i \\\\
V_{i; 10}h_{s,i,w}  & \frac{1}{10}L_i < \sum_{j=1}^{w-1}h_{s,i,w}\leq \frac{3}{10}L_i  \\\\
V_{i; 30}h_{s,i,w}  & \frac{3}{10}L_i < \sum_{j=1}^{w-1}h_{s,i,w}
\end{array}\right.
\end{equation}$$

Then the delivery cost in year $y$ and scenario $s$ is

$$\begin{equation}
D_{s,y}^{\text{(S-del)}} = 
\sum_{w\in y} \sum_i H_{s,i,w}
\qquad(3.5)
\end{equation}$$

<!-- We set  -->
<!-- $$V_{LLMIC; j} = \frac{1}{N_{LMIC}^{(15)} + N_{LIC}^{(15)}} \left(N_{LMIC}^{(15)}V_{LMIC; j} + N_{LIC}^{(15)}V_{LIC; j} \right)$$ -->

<div class="figure">

<img src="README_files/figure-gfm/deliverycost-1.png" alt="SSV delivery cost"  />
<p class="caption">
<span id="fig:deliverycost"></span>Figure 3.4: SSV delivery cost
</p>

</div>

     BAU Min.   :106597   1st Qu.:121357   Median :136500   Mean   :138771  
     S01 Min.   :106924   1st Qu.:121657   Median :136762   Mean   :139059  
     S02 Min.   :106669   1st Qu.:121430   Median :136558   Mean   :138837  
     S03 Min.   :106761   1st Qu.:121519   Median :136626   Mean   :138915  
     S04 Min.   :107261   1st Qu.:121996   Median :137045   Mean   :139367  
     S05 Min.   :107525   1st Qu.:122259   Median :137248   Mean   :139599  
     S06 Min.   :107702   1st Qu.:122439   Median :137397   Mean   :139764  
     S07 Min.   :109228   1st Qu.:124100   Median :138597   Mean   :141279  
     S08 Min.   :108483   1st Qu.:123260   Median :137997   Mean   :140468  
     S09 Min.   :108719   1st Qu.:123506   Median :138179   Mean   :140680  
     S10 Min.   :105737   1st Qu.:120443   Median :135796   Mean   :137847  
                                          
     BAU 3rd Qu.:155601   Max.   :175655  
     S01 3rd Qu.:155878   Max.   :175903  
     S02 3rd Qu.:155668   Max.   :175707  
     S03 3rd Qu.:155748   Max.   :175769  
     S04 3rd Qu.:156184   Max.   :176159  
     S05 3rd Qu.:156416   Max.   :176346  
     S06 3rd Qu.:156578   Max.   :176487  
     S07 3rd Qu.:158275   Max.   :178027  
     S08 3rd Qu.:157277   Max.   :177099  
     S09 3rd Qu.:157484   Max.   :177280  
     S10 3rd Qu.:154493   Max.   :174588  

Targets:

114,526 ( 90,654 110,005 134,444 ) 114,771 ( 91,321 111,130 134,341 )
114,769 ( 91,620 110,604 133,752 ) 114,811 ( 91,647 110,815 133,856 )
114,527 ( 91,170 110,664 133,720 ) 114,615 ( 91,074 110,653 133,836 )
115,095 ( 91,858 111,205 134,355 ) 115,634 ( 92,514 111,639 134,375 )
116,385 ( 93,116 112,183 135,664 ) 117,196 ( 93,427 113,114 136,861 )
116,913 ( 93,536 112,957 136,414 ) 118,141 ( 94,682 114,649 137,100 )
113,540 ( 89,745 109,012 132,595 )

### 3.3.2 BPSV

**These values match the spreadsheet results. (NB: more doses are
purchased and delivered than there are eligible people in the
population)**

For the BPSV, which goes only to people aged 65 or older, with
populations $N_i^{(65)}$, coverage is reached earlier in the process, so
the cost is weighted more heavily towards start up and ramp up:

$$\begin{equation}
D_s^{\text{(BP-del)}} = 
\left\\{\begin{array}{lr}
\sum_{i}D_{\text{BPSV},i}
\\; & \\; s=1 \\\\
0  \\; & \\; s\neq 1
\end{array}\right.
\qquad(3.6)
\end{equation}$$

$$\begin{equation}
D_{\text{BPSV},i} = 
\left\\{\begin{array}{lr}
N_i^{(65)}V_{i; 0}  & N_i^{(65)}\leq \frac{1}{10}N_i^{(15)} \\\\
\frac{N_i^{(15)}}{10} V_{i; 0} + \left(N_i^{(65)}-\frac{N_i^{(15)}}{10} \right)V_{i; 11}  & \frac{1}{10}N_i^{(15)} < N_i^{(65)}\leq \frac{3}{10}N_i^{(15)} \\\\
\frac{N_i^{(15)}}{10} V_{i; 0} + \frac{2}{10}N_i^{(15)} V_{i; 11} + \left(N_i^{(65)}-\frac{3}{10}N_i^{(15)} \right)V_{i; 31} & N_i^{(65)}> \frac{3}{10} N_i^{(15)}
\end{array}\right.
\end{equation}$$

The logic of this is as follows:

- The increments in cost correspond to numbers of eligible people in the
  whole population, namely those aged 15 and above.
- If the number of people eligible for the BPSV is less than 10% of the
  population aged 15 and over, then all doses cost the “start up”
  amount.
- If the number of people eligible for the BPSV is more than 10% and
  less than 30% of the 15+ population, then cost of the first doses, a
  number equal to 10% of the 15+ population, is the “start up” amount.
  All remaining doses cost the “ramp up” amount.
- If the number of people eligible for the BPSV is more than 30% of the
  15+ population, then the cost of the first doses, a number equal to
  10% of the 15+ population, is the “start up” amount. The cost of the
  second tranche of doses, a number equal to 20% of the 15+ population,
  is the “ramp up” amount. All remaining doses cost the “getting to
  scale” amount.

<div class="figure">

<img src="README_files/figure-gfm/bpsvdeliverycost-1.png" alt="BPSV delivery cost"  />
<p class="caption">
<span id="fig:bpsvdeliverycost"></span>Figure 3.5: BPSV delivery cost
</p>

</div>

Min. 1st Qu. Median Mean 3rd Qu. Max. 11993 12858 13409 14118 15611
17557

Target: 11,206 (9,037 10,865 13,054)

| Country | Country status | Study type | Financial Cost per dose (USD) | Source |
|:---|----|:---|----|:---|
| WHO, Gavi, and UNICEF AMC Estimate | AMC | Top down | 1.66 | Griffiths et al. (2021) |
| UNICEF Global Estimate | All | Model | 0.73 | Oyatoye (2023) |
| DRC | LIC | Bottom up | 1.91 | Moi et al. (2024) |
| Malawi | LIC | Bottom up | 4.55 | Ruisch et al. (2025) |
| Mozambique | LIC | Bottom up | 0.5 | Namalela et al. (2025) |
| Uganda | LIC | Bottom up | 0.79 | Tumusiime et al. (2024) |
| Bangladesh | LMIC | Bottom up | 0.29 | Yesmin et al. (2024) |
| Cote d’Ivoire | LMIC | Bottom up | 0.67 | K. Vaughan et al. (2023) |
| Nigeria | LMIC | Bottom up | 0.84 | Noh et al. (2024) |
| Philippines | LMIC | Bottom up | 2.16 | Banks et al. (2023) |
| Vietnam | LMIC | Bottom up | 1.73 | Nguyen et al. (2024) |
| Ghana | LMIC | CVIC tool | 2.2–2.3 | Nonvignon et al. (2022) |
| Lao PDR | LMIC | CVIC tool | 0.79–0.81 | Yeung et al. (2023) |
| Kenya | LMIC | Top down | 3.29–4.28 | Orangi et al. (2022) |
| Botswana | UMIC | Mixed | 19 | Kelsey Vaughan et al. (2025) |
| South Africa | UMIC | Top down | 3.84 | Edoka et al. (2024) |

<span id="tab:delcosts"></span>Table 3.1: Literature review of global
and country-specific delivery costs

# 4 SSV delivery

| Category | Reserved capacity | Private response (existing capacity) | Private response (built capacity) |
|:---|----|:---|:---|
| Annual manufacturing volume | By scenario (0.5–2.5B) | 9B minus reserved volume | 6B |
| Facility transition start | 7 weeks before vaccine approval | 7 weeks before vaccine approval | 7 weeks before vaccine approval |
| Weeks to initial manufacturing | 12 | 12 (BPSV) or 30 (no BPSV) | 48 |
| Scale-up weeks to full capacity | 10 | 16 | 16 |

Manufacturing response timeline assumptions

<!-- | Weeks from transition start | 0-11 | 12-21 | 22-29 | 30-45  | 46-47 | 48-63 | 64+ | -->
<!-- |---|---|---|---|---|---|---|---| -->
<!-- | Reserved Capacity (%)  || Scaling from 0-100 | 100 | 100 | 100 | 100 | 100 |  -->
<!-- | Private Capacity (Existing; %)  || | | Scaling from 0-100 | 100 | 100 | 100 |  -->
<!-- | Private Capacity (Response; %)  | | | |  | | Scaling from 0-100 | 100 | -->

| Weeks from transition start | Reserved Capacity (%) | Existing Private Capacity (%) | Response Private Capacity (%) |
|:---|:---|:---|:---|
| 0–11 |  |  |  |
| 12–21 | Scales from 0 to 100 |  |  |
| 22–29 | 100 |  |  |
| 30–45 | 100 | Scales from 0 to 100 |  |
| 46–47 | 100 | 100 |  |
| 48–63 | 100 | 100 | Scales from 0 to 100 |
| 64+ | 100 | 100 | 100 |

Vaccine Production Timeline when there is no BPSV. When BPSV is also
modelled, Existing Private Capacity scales from 0 to 100 in weeks 12–21.

## 4.1 Timing

Facility transition occurs $I_0=7$ weeks before vaccine approval, which
in turn depends on R&D investments. We have three levels in our
scenarios, corresponding to SSVs available in 100 days, 200 days, and
365 days. The total weeks taken for vaccine approval can be written as
follows:

$$W_{j}^{(S)} = \sum_{i=0}^3 W_{i;j}^{(S)}$$

for $j\in\lbrace 365, 200, 100\rbrace$. These work out as 52, 29, and 14
weeks, respectively. Thus “week 0” for manufacturing occurs 45, 22, and
7 weeks, respectively, after the new pathogen has been sequenced. We
denote this variable $w_s^{(0)}$.

## 4.2 Production

The total global manufacturing volume is $M_G=15$ billion doses. The
amount that is reserved, in billion doses, including the HIC-specific
reservation of $A_3=0.5$ billion doses, depends on the scenarios as
follows:

$$\begin{equation}
M_{R,s} = \left\\{\begin{array}{lr}A_3 & s\in\\{0, 1, 4, 7, 10\\} \\\\ 
A_3 + 0.7 & s\in\\{2, 5, 8\\} \\\\ 
A_3 + 2 & s\in\\{3, 6, 9\\} \end{array}\right.
\end{equation}$$

where $s=0$ denotes the BAU scenario. By definition,
$M_{E,s} = M_C - M_{R,s}$, and $M_B=M_G-M_C$.

Then the number of doses, in billions, that are made from capacity
$x\in \lbrace R, E, B\rbrace$ in week $w$ of scenario $s$ is:

$$\begin{equation}
Z_{x,s,w} = \left\\{\begin{array}{lr}0 & w-w_s^{(0)} \leq I_x \\\\ 
\frac{1}{52}\frac{w-w_s^{(0)}-I_x}{C_x}M_{x,s} & w-w_s^{(0)}\in(I_x, I_x+C_x] \\\\ 
\frac{1}{52}M_{x,s}  & w-w_s^{(0)}> I_x+C_x
\end{array}\right.
\end{equation}$$

<!-- \frac{1}{52}M_{R,s}  & w\in[I_R+C_R, I_E)\\\\  -->
<!-- \frac{1}{52}\left(M_{R,s} + \frac{w-I_E+1}{C_E}M_{E,s}\right) & w\in[I_E, I_E+C_E)\\\\  -->
<!-- \frac{1}{52}\left(M_{R,s} + M_{E,s}\right)  & w\in[I_E+C_E, I_B)\\\\  -->
<!-- \frac{1}{52}\left(M_{R,s} + M_{E,s} + \frac{w-I_B+1}{C_B}M_{B}\right) & w\in[I_B, I_B+C_B)\\\\  -->

where $I_R = 12$ is the number of weeks to initial manufacturing for
reserved capacity, $C_R = 10$ is the number of weeks to scale up to full
capacity; $I_B = 48$ is the number of weeks to initial manufacturing for
built and unreserved capacity, $C_B = 16$ is the number of weeks to
scale up to full capacity, and

$$\begin{equation}
I_E = \left\\{\begin{array}{lr}
 I_{E,1} \\; & \\; s=1 \\\\
I_{E,0}  \\; & \\; s\neq 1
\end{array}\right.
\end{equation}$$

where $I_{E,0} = 30$ and $I_{E,1} = 12$ are the number of weeks to
initial manufacturing for existing and unreserved capacity, $C_E = 16$
is the number of weeks to scale up to full capacity.

<!-- Then the total number of doses produced in week $w$ is  -->
<!-- \begin{equation} -->
<!-- Z_{T,s,w} = Z_{R,s,w}+Z_{E,s,w}+Z_{B,s,w}. -->
<!-- \end{equation} -->

<div class="figure">

<img src="README_files/figure-gfm/supply-1.png" alt="Doses made available from manufacturing per scenario. Weeks are in reference to the sequencing of the pathogen."  />
<p class="caption">
<span id="fig:supply"></span>Figure 4.1: Doses made available from
manufacturing per scenario. Weeks are in reference to the sequencing of
the pathogen.
</p>

</div>

In Figure <a href="#fig:supply">4.1</a>, the following scenarios have
identical supply (because they have the same capacity reservations and
R&D investments): BAU & S10.

## 4.3 Allocation

Denote the weekly allocated doses at week $w$ from capacity $x$ to
income level $i$ $k_{s,x,i,w}$, and the cumulative number $K_{s,i,w}$,
such that

$$K_{s,i,w} = \sum_{x\in\lbrace R,E,B\rbrace}\sum_{j=0}^w k_{s,x,i,j}.$$

<!-- assuming vaccine wastage of $\delta = 0.31$. -->

$$\begin{equation}
k_{s,R,i,w} = \left\\{ \begin{array}{lr}
\frac{A_3}{M_{R,s}}Z_{R,s,w}             & K_{s,\text{HIC},w} <  L_{\text{HIC}}\\;\\&\\; i=\text{HIC} \\\\
\frac{M_{R,s}-A_3}{M_{R,s}}\frac{N_{i}}{N_{HIC}+N_{UMIC}+N_{LMIC}+N_{LIC}}Z_{R,s,w}                    & K_{s,\text{HIC},w} < L_{\text{HIC}}  \\\\
0 & K_{s,\text{HIC},w} \geq L_{\text{HIC}} \\;\\&\\; i=\text{HIC}\\\\
\frac{N_{i}}{N_{UMIC}+N_{LMIC}+N_{LIC}}Z_{R,s,w} & K_{s,\text{HIC},w} \geq L_{\text{HIC}} \\;\\&\\;  K_{s,\text{UMIC},w} < L_{\text{UMIC}} \\;\\&\\; i\neq\text{HIC}\\\\
0 & K_{s,\text{UMIC},w} \geq L_{\text{UMIC}} \\;\\&\\; i=\text{UMIC}\\\\
\frac{N_{i}}{N_{LMIC}+N_{LIC}}Z_{R,s,w} & K_{s,\text{UMIC},w} \geq L_{\text{UMIC}} \\;\\&\\;  K_{s,\text{LMIC},w} < L_{\text{LMIC}} \\;\\&\\; i\notin  \lbrace\text{HIC},\text{UMIC}\rbrace\\\\
0 & K_{s,\text{LMIC},w} \geq L_{\text{LMIC}} \\;\\&\\; i=\text{LMIC}\\\\
Z_{R,s,w}             & K_{s,\text{LMIC},w} \geq L_{\text{LMIC}} \\;\\&\\; K_{s,\text{LIC},w} < L_{\text{LIC}} \\;\\&\\; i=\text{LIC}\\\\
0 & K_{s,\text{LIC},w} \geq L_{\text{LIC}} 
\end{array}\right.
\end{equation}$$

The logic of this reads as follows:

- $A_3=0.5$ billion doses per year from reserved capacity go exclusively
  to HIC, which is expressed as a fraction of the total reservation,
  $M_{R,s}$
- Any remaining reserved capacity doses are allocated according to
  population
- Once HIC reach their total demand, doses from reserved capacity are
  split proportional to population between UMIC, LMIC and LIC, and so on

For $x\in\lbrace E,B\rbrace$,

$$\begin{equation}
k_{s,x,i,w} = \left\\{ \begin{array}{lr}
Z_{x,s,w}            & K_{s,\text{HIC},w} < L_{\text{HIC}} \\;\\&\\; i=\text{HIC} \\\\
0                     & K_{s,\text{HIC},w} < L_{\text{HIC}} \\;\\&\\; i\neq\text{HIC} \\\\
Z_{x,s,w}            & K_{s,\text{HIC},w} \geq L_{\text{HIC}} \\;\\&\\; K_{s,\text{UMIC},w} < L_{\text{UMIC}} \\;\\&\\; i=\text{UMIC} \\\\
0                     & K_{s,\text{HIC},w} \geq L_{\text{HIC}} \\;\\&\\; K_{s,\text{UMIC},w} < L_{\text{UMIC}} \\;\\&\\; i\neq\text{UMIC} \\\\
Z_{x,s,w}            & K_{s,\text{UMIC},w} \geq L_{\text{UMIC}} \\;\\&\\; K_{s,\text{LMIC},w} < L_{\text{LMIC}} \\;\\&\\; i=\text{LMIC} \\\\
0                     & K_{s,\text{UMIC},w} \geq L_{\text{UMIC}} \\;\\&\\; K_{s,\text{LMIC},w} < L_{\text{LMIC}} \\;\\&\\; i\neq\text{LMIC} \\\\
Z_{x,s,w}            & K_{s,\text{LMIC},w} \geq L_{\text{LMIC}} \\;\\&\\; i=\text{LIC} \\\\
0                     & K_{s,\text{LMIC},w} \geq L_{\text{LMIC}} \\;\\&\\; i\neq\text{LIC} 
\end{array}\right.
\end{equation}$$

The logic of this reads as follows:

- Until HIC demand is reached, all doses from unreserved capacity go to
  HIC. None go to UMIC, LMIC and LIC.
- Once HIC demand has been met and until UMIC demand is reached, all
  doses from unreserved capacity go to UMIC. None go to HIC, LMIC and
  LIC.
- Once HIC and UMIC demand have been met and until LMIC demand is
  reached, all doses from unreserved capacity go to LMIC. None go to
  HIC, UMIC and LIC.
- Once HIC, UMIC and LMIC demand have been met, all remaining doses from
  unreserved capacity go to LIC. None go to LMIC, UMIC and HIC.

Total supply of first-schedule doses in each year period is

$$\begin{equation}
A_{x,s,y}^{(1)} = \sum_i\sum_{w\in y}k_{s,x,i,w}.
\end{equation}$$

We assume, for every second dose of SSV, a booster will be given one
year later for $N^{\text{(boost)}} = 2$ years.

Thus

$$\begin{equation}
A_{s,y}^{(2)} = \left\\{ \begin{array}{lr}
\frac{1}{2}\sum_{x}A_{x,s,y-1}^{(1)}            & y=2 \\\\
\frac{1}{2}\sum_{x}\left(A_{x,s,y-1}^{(1)} + A_{x,s,y-2}^{(1)}\right)            & y>2 
\end{array}\right.
\end{equation}$$

and

$$\begin{equation}
A_{x,s,y}^{(2)} = \min(A_{s,y}^{(2)}, M_R - A_{R,s,y}^{(1)})
\end{equation}$$

for $x=R$ and

$$\begin{equation}
A_{x,s,y}^{(2)} = \max(A_{s,y}^{(2)} - A_{R,s,y}^{(2)}, 0)
\end{equation}$$

for $x\in\lbrace E, B \rbrace$.

Then

$$\begin{equation}
A_{x,s,y} = A_{x,s,y}^{(1)} + A_{x,s,y}^{(2)}.
\qquad(4.1)
\end{equation}$$

<div class="figure">

<img src="README_files/figure-gfm/procurement-1.png" alt="Doses procured by country income level"  />
<p class="caption">
<span id="fig:procurement"></span>Figure 4.2: Doses procured by country
income level
</p>

</div>

## 4.4 Delivery

**These values do not look correct**

Delivery is written as $h_{s,i,w}^{(j)}$ doses delivered in scenario
$s$, income group $i$ and week $w$ for schedule $j$, which is $j=1$ for
first-dose SSV, $j=2$ for second-dose SSV, and $j=2+k$ for
$k=\lbrace 1,...,N^{(boost)}\rbrace$.

The second dose is prioritised over the first, and follows the first by
four weeks, so

$$
h_{s,i,w}^{(2)} = h_{s,i,w-4}^{(1)}.
$$

and

$$
h_{s,i,w}^{(2)} = h_{s,i,w-4}^{(1)}.
$$

Available doses are $K_{s,i,w}$, the cumulative doses allocated to
income group $i$ by week $w$, minus doses given so far. First doses stop
being given once $L_i/2$ is reached.

$$\begin{equation}
h_{s,i,w}^{(1)} =
\max\left\lbrace 0
\right\rbrace
\end{equation}$$

<div class="figure">

<img src="README_files/figure-gfm/scendelivery-1.png" alt="Cumulative vaccine coverage (second SSV dose) by country income level"  />
<p class="caption">
<span id="fig:scendelivery"></span>Figure 4.3: Cumulative vaccine
coverage (second SSV dose) by country income level
</p>

</div>

# 5 BPSV delivery

## 5.1 Timing

The duration of the Phase III trial is $W_3^{(B)} = 18$ weeks. The time
to manufacturing transition is $I_R = 12$ weeks, and the time to
manufacturing scale-up $C_R = 10$ weeks; these are the same as the
reserved-capacity times for SSV.

Facility transition occurs in week 1. Thus manufacturing begins in week
$1+I_R = 13$ and dose distribution begins in week $1+W_3^{(B)} = 19$.

## 5.2 Production

The number of doses, in billions, that are made in week $w$ is:

$$\begin{equation}
Z_{w} = \left\\{\begin{array}{lr}0 & w < I_R \\\\ 
\frac{1}{52}\frac{w-I_x+1}{C_x}M_{x,s} & w\in[I_R, I_R+C_R) \\\\ 
\frac{1}{52}M_{x,s}  & w-1\geq I_R+C_xR
\end{array}\right.
\end{equation}$$

<div class="figure">

<img src="README_files/figure-gfm/bpsvsupply-1.png" alt="BPSV doses made available from manufacturing per scenario. Weeks are in reference to the sequencing of the pathogen."  />
<p class="caption">
<span id="fig:bpsvsupply"></span>Figure 5.1: BPSV doses made available
from manufacturing per scenario. Weeks are in reference to the
sequencing of the pathogen.
</p>

</div>

## 5.3 Allocation

Doses are all allocated in proportion to the eligible population.

<div class="figure">

<img src="README_files/figure-gfm/bpsvprocurement-1.png" alt="BPSV doses procured by country income level"  />
<p class="caption">
<span id="fig:bpsvprocurement"></span>Figure 5.2: BPSV doses procured by
country income level
</p>

</div>

## 5.4 Delivery

**These values do not look correct**

<div class="figure">

<img src="README_files/figure-gfm/bpsvdeliveryplot-1.png" alt="BPSV vaccine coverage by country income level"  />
<p class="caption">
<span id="fig:bpsvdeliveryplot"></span>Figure 5.3: BPSV vaccine coverage
by country income level
</p>

</div>

# 6 Parameter samples

![](README_files/figure-gfm/parsamples-1.png)<!-- -->![](README_files/figure-gfm/parsamples-2.png)<!-- -->![](README_files/figure-gfm/parsamples-3.png)<!-- -->![](README_files/figure-gfm/parsamples-4.png)<!-- -->![](README_files/figure-gfm/parsamples-5.png)<!-- -->![](README_files/figure-gfm/parsamples-6.png)<!-- -->![](README_files/figure-gfm/parsamples-7.png)<!-- -->![](README_files/figure-gfm/parsamples-8.png)<!-- -->![](README_files/figure-gfm/parsamples-9.png)<!-- -->![](README_files/figure-gfm/parsamples-10.png)<!-- -->![](README_files/figure-gfm/parsamples-11.png)<!-- -->![](README_files/figure-gfm/parsamples-12.png)<!-- -->![](README_files/figure-gfm/parsamples-13.png)<!-- -->![](README_files/figure-gfm/parsamples-14.png)<!-- -->![](README_files/figure-gfm/parsamples-15.png)<!-- -->![](README_files/figure-gfm/parsamples-16.png)<!-- -->![](README_files/figure-gfm/parsamples-17.png)<!-- -->![](README_files/figure-gfm/parsamples-18.png)<!-- -->![](README_files/figure-gfm/parsamples-19.png)<!-- -->![](README_files/figure-gfm/parsamples-20.png)<!-- -->![](README_files/figure-gfm/parsamples-21.png)<!-- -->![](README_files/figure-gfm/parsamples-22.png)<!-- -->![](README_files/figure-gfm/parsamples-23.png)<!-- -->![](README_files/figure-gfm/parsamples-24.png)<!-- -->![](README_files/figure-gfm/parsamples-25.png)<!-- -->![](README_files/figure-gfm/parsamples-26.png)<!-- -->![](README_files/figure-gfm/parsamples-27.png)<!-- -->![](README_files/figure-gfm/parsamples-28.png)<!-- -->

# 7 Contributors

Model: Peter Windus, Andy Torkelson

Data: Peter Windus, Andy Torkelson, Damian Walker

Documentation: Peter Windus, Andy Torkelson, Rob Johnson

R code: Rob Johnson

# 8 References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Banks2023" class="csl-entry">

Banks, C, RD Estanislao, SJ De los Reyes, JE De Guzman, LB
Sumpaico-Tanchanco, B Makani-Lim, R Archer, and L Boonstoppel. 2023.
“The Cost of Delivering COVID-19 Vaccines in the Philippines.” Geneva:
ThinkWell.
<https://thinkwell.global/wp-content/uploads/2024/03/Cost-of-delivering-COVID19-vaccines-in-the-Philippines_final-report_19-Dec-2023.pdf>.

</div>

<div id="ref-CEPI2021" class="csl-entry">

CEPI. 2021. “CEPI 2022–2026 Strategy.”
<https://static.cepi.net/downloads/2023-12/CEPI-2022-2026-Strategy-v3-Jan21_0.pdf>.

</div>

<div id="ref-CEPI2022" class="csl-entry">

———. 2022. “Delivering Pandemic Vaccines in 100 Days.”
<https://static.cepi.net/downloads/2024-02/CEPI-100-Days-Report-Digital-Version_29-11-22.pdf>.

</div>

<div id="ref-Edoka2024" class="csl-entry">

Edoka, Ijeoma, Lineo Marie Matsela, Khumo Modiba, Yolandie Luther,
Sharlene Govender, Thapelo Maotoe, Heena Brahmbhatt, Pedro T. Pisa,
Gesine Meyer-Rath, and Jacqui Miot. 2024. “Costs of the COVID-19
Vaccination Programme: Estimates from the West Rand District of South
Africa, 2021/2022.” *BMC Health Services Research* 24 (1): 857.
<https://doi.org/10.1186/s12913-024-11251-1>.

</div>

<div id="ref-Glennerster2023" class="csl-entry">

Glennerster, Rachel, Christopher M. Snyder, and Brandon Joel Tan. 2023.
*Calculating the Costs and Benefits of Advance Preparations for Future
Pandemics*. *IMF Economic Review*. Vol. 71. Palgrave Macmillan UK.
<https://doi.org/10.1057/s41308-023-00212-z>.

</div>

<div id="ref-Gouglas2018" class="csl-entry">

Gouglas, Dimitrios, Tung Thanh Le, Klara Henderson, Aristidis Kaloudis,
Trygve Danielsen, Nicholas Caspersen Hammersland, James M Robinson,
Penny M Heaton, and John-Arne Røttingen. 2018. “Estimating the Cost of
Vaccine Development Against Epidemic Infectious Diseases: A Cost
Minimisation Study.” *The Lancet Global Health* 6 (12): e1386–96.
<https://doi.org/10.1016/S2214-109X(18)30346-2>.

</div>

<div id="ref-Griffiths2021" class="csl-entry">

Griffiths, Ulla, Alex Adjagba, Marcia Attaran, Raymond Hutubessy,
Nathalie Van De Maele, Karene Yeung, Wei Aun, et al. 2021. “Costs of
Delivering COVID-19 Vaccine in 92 AMC Countries.” February. COVAX
Working Group. <https://www.who.int/publications/i/item/10665337553>.

</div>

<div id="ref-Kazaz2021" class="csl-entry">

Kazaz, Burak. 2021. “Incentivizing COVID-19 Vaccine Developers to Expand
Manufacturing Capacity.” Center for Global Development.
<https://www.cgdev.org/sites/default/files/incentivizing-covid-19-vaccine-developers-expand-manufacturing-capacity.pdf>.

</div>

<div id="ref-LinksbridgeSPC2025" class="csl-entry">

Linksbridge SPC. 2025. “Global Vaccine Market Model.”
<https://4550bf57-cdn.agilitycms.cloud/help-guides/Introduction%20to%20GVMM%20v6.1.pdf>.

</div>

<div id="ref-Moi2024" class="csl-entry">

Moi, Flavia, Laura Boonstoppel, Rachel Archer, and Pierre Akilimali.
2024. “The Cost of Delivering COVID-19 Vaccines in the Democratic
Republic of the Congo.” Geneva: ThinkWell.
<https://thinkwell.global/wp-content/uploads/2024/04/DRC-C19-costing-study-report_final.pdf>.

</div>

<div id="ref-Namalela2025" class="csl-entry">

Namalela, Tozé, Flavia Moi, Amélia Dipuve, Pedro Marizane Pota, José
Guambe, Maria Tereza Couto, and Laura Boonstoppel. 2025. “The Cost of
Delivering COVID-19 Vaccines in Mozambique: A Bottom-up Costing Study.”
*BMC Health Services Research* 25 (1): 521.
<https://doi.org/10.1186/s12913-025-12671-3>.

</div>

<div id="ref-Nguyen2024" class="csl-entry">

Nguyen, Van Minh, Flavia Moi, Laura Boonstoppel, Hong Thi Duong, Chien
Chinh Vien, and Minh Van Hoang. 2024. “The Cost of Delivering COVID-19
Vaccines in Vietnam.” *BMC Health Services Research* 24 (1): 779.
<https://doi.org/10.1186/s12913-024-11202-w>.

</div>

<div id="ref-Noh2024" class="csl-entry">

Noh, Dave Haeyun, Roopa Darwar, Belinda V. Uba, Shiva Gab-deedam, Stella
Yani, Akolade Jimoh, Ndadilnasiya Waziri, et al. 2024. “Cost of COVID-19
Vaccine Delivery in Nine States in Nigeria via the U.S. Government
Initiative for Global Vaccine Access.” *BMC Health Services Research* 24
(1): 1232. <https://doi.org/10.1186/s12913-024-11645-1>.

</div>

<div id="ref-Nonvignon2022" class="csl-entry">

Nonvignon, Justice, Richmond Owusu, Brian Asare, Alex Adjagba, Yap Wei
Aun, Karene Hoi Ting Yeung, Joycelyn Naa Korkoi Azeez, et al. 2022.
“Estimating the Cost of COVID-19 Vaccine Deployment and Introduction in
Ghana Using the CVIC Tool.” *Vaccine* 40 (12): 1879–87.
<https://doi.org/10.1016/j.vaccine.2022.01.036>.

</div>

<div id="ref-Orangi2022" class="csl-entry">

Orangi, Stacey, Angela Kairu, Anthony Ngatia, John Ojal, and Edwine
Barasa. 2022. “Examining the Unit Costs of COVID-19 Vaccine Delivery in
Kenya.” *BMC Health Services Research* 22 (1): 439.
<https://doi.org/10.1186/s12913-022-07864-z>.

</div>

<div id="ref-Oyatoye2023" class="csl-entry">

Oyatoye, I. 2023. “Costs and Financing Gap of Delivering COVID-19
Vaccine to 133 Low- and Middle-Income Countries.” Cape Town, South
Africa.
<https://immunizationeconomics.org/wp-content/uploads/2024/01/Ibironke-Oyatoye-Costs-and-financing-gap-of-C19v-delivery-in-LMICS_Final.pdf>.

</div>

<div id="ref-Pfizer2023" class="csl-entry">

Pfizer. 2023. “Pfizer and the European Commission Enter into
Manufacturing Reservation Agreement for <span class="nocase">mRNA-based
Vaccines</span> to Help Protect Against Future Pandemics.”
<https://www.pfizer.com/news/announcements/pfizer-and-european-commission-enter-manufacturing-reservation-agreement-mrna>.

</div>

<div id="ref-Ruisch2025" class="csl-entry">

Ruisch, Anika, Simon Ntopi, Ishani Mathur, Maeve Conlin, Anna McCaffrey,
Damian G. Walker, and Christian Suharlim. 2025. “The Cost of Delivering
COVID-19 Vaccines in Four Districts in Malawi.” *Cost Effectiveness and
Resource Allocation* 23 (1): 36.
<https://doi.org/10.1186/s12962-025-00610-2>.

</div>

<div id="ref-Tumusiime2024" class="csl-entry">

Tumusiime, Cathbert, Rachel Archer, Charlotte Muheki, Paul Kiggundu,
Richard Ssemujju, Angellah Nakyanzi, Prossy Kiddu Namyalo, et al. 2024.
“The Cost of Delivering COVID-19 Vaccines in Kampala, Uganda.” Uganda:
ThinkWell.
<https://immunizationeconomics.org/wp-content/uploads/2024/05/thinkwell-report-uganda-final.pdf>.

</div>

<div id="ref-U.S.BLS2025" class="csl-entry">

U.S. BLS. 2025. “CPI Inflation Calculator.”
<https://www.bls.gov/data/inflation_calculator.htm>.

</div>

<div id="ref-VaccinesEurope2023" class="csl-entry">

Vaccines Europe. 2023. “Vaccines Europe Analysis of Vaccine Production
Lead Times.”
<https://www.cgdev.org/sites/default/files/incentivizing-covid-19-vaccine-developers-expand-manufacturing-capacity.pdf>.

</div>

<div id="ref-Vaughan2025" class="csl-entry">

Vaughan, Kelsey, Onalenna T. Mokena, Goabaone Rankgoane-Pono, Moses
Keetile, and Ulla Kou Griffiths. 2025. “Costs of Delivering COVID-19
Vaccine in Botswana During the Height of the Pandemic: A Retrospective
Study.” *BMC Health Services Research* 25 (1): 405.
<https://doi.org/10.1186/s12913-025-12455-9>.

</div>

<div id="ref-Vaughan2023" class="csl-entry">

Vaughan, K, E Smith, C Schütte, F Moi, and L Boonstoppel. 2023. “The
Cost of Delivering COVID-19 Vaccines in
<span class="nocase">C<span class="nocase">ô</span>te</span> d’Ivoire.”
ThinkWell & Genesis Analytics.
<https://thinkwell.global/wp-content/uploads/2023/09/Cote-dIvoire-final-report_FINAL.pdf>.

</div>

<div id="ref-Wong2019" class="csl-entry">

Wong, Chi Heem, Kien Wei Siah, and Andrew W Lo. 2019. “Estimation of
Clinical Trial Success Rates and Related Parameters.” *Biostatistics* 20
(2): 273–86. <https://doi.org/10.1093/biostatistics/kxx069>.

</div>

<div id="ref-Yesmin2024" class="csl-entry">

Yesmin, Afroja, Flavia Moi, Tarek Hossain, Rachel A. Archer, Monjurul
Islam, and Laura Boonstoppel. 2024. “The Cost of COVID-19 Vaccine
Delivery in Bangladesh.” *Human Vaccines & Immunotherapeutics* 20 (1):
2411820. <https://doi.org/10.1080/21645515.2024.2411820>.

</div>

<div id="ref-Yeung2023" class="csl-entry">

Yeung, Karene Hoi Ting, Eunkyoung Kim, Wei Aun Yap, Chansay
Pathammavong, Lauren Franzel, Yu Lee Park, Peter Cowley, Ulla Kou
Griffiths, and Raymond Christiaan W. Hutubessy. 2023. “Estimating the
Delivery Costs of COVID-19 Vaccination Using the COVID-19 Vaccine
Introduction and Deployment Costing (CVIC) Tool: The Lao People’s
Democratic Republic Experience.” *BMC Medicine* 21 (1): 248.
<https://doi.org/10.1186/s12916-023-02944-1>.

</div>

</div>
