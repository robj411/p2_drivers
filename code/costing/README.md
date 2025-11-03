The Costing Model
================

- [1 Overview](#1-overview)
- [2 Parameters](#2-parameters)
- [3 Model details](#3-model-details)
  - [3.1 Preparedness cost equation (annual calculation,
    2025-2039):](#31-preparedness-cost-equation-annual-calculation-2025-2039)
  - [3.2 Response cost equation (annual calculation,
    2040-2045):](#32-response-cost-equation-annual-calculation-2040-2045)
  - [3.3 Risk-adjusted R&D cost per candidate calculation
    -](#33-risk-adjusted-rd-cost-per-candidate-calculation--)
  - [3.4 Procurement cost
    calculation:](#34-procurement-cost-calculation)
  - [3.5 Delivery Cost Equation:](#35-delivery-cost-equation)
  - [3.6 Vaccination Scenarios](#36-vaccination-scenarios)
- [4 Results](#4-results)
- [5 Attributions / Authors](#5-attributions--authors)

<!-- # Figures (temporary) {.unlisted .unnumbered} -->

This document describes the costing model that is used in the CEPI
application.

# 1 Overview

# 2 Parameters

| Math notation | Description | Distribution | Parameter 1 | Parameter 2 | Parameter 3 | Parameter 4 | Parameter 5 |
|:--:|:--:|:--:|----|----|----|----|----|
| $T_{0; 365}^{(S)}$ | SSV preclinical duration (365); weeks | Constant | 14 |  |  |  |  |
| $T_{0; 200}^{(S)}$ | SSV preclinical duration (200DM); weeks | Constant | 5 |  |  |  |  |
| $T_{0; 100}^{(S)}$ | SSV preclinical duration (100DM); weeks | Constant | 5 |  |  |  |  |
| $T_{1; 365}^{(S)}$ | SSV phase I duration (365); weeks | Constant | 19 |  |  |  |  |
| $T_{1; 200}^{(S)}$ | SSV phase I duration (200DM); weeks | Constant | 7 |  |  |  |  |
| $T_{1; 100}^{(S)}$ | SSV phase I duration (100DM); weeks | Constant | 0 |  |  |  |  |
| $T_{2; 365}^{(S)}$ | SSV phase II duration (365); weeks | Constant | 19 |  |  |  |  |
| $T_{2; 200}^{(S)}$ | SSV phase II duration (200DM); weeks | Constant | 0 |  |  |  |  |
| $T_{2; 100}^{(S)}$ | SSV phase II duration (100DM); weeks | Constant | 0 |  |  |  |  |
| $T_{3; 365}^{(S)}$ | SSV phase III duration (365); weeks | Constant | 16 |  |  |  |  |
| $T_{3; 200}^{(S)}$ | SSV phase III duration (200DM); weeks | Constant | 15 |  |  |  |  |
| $T_{3; 100}^{(S)}$ | SSV phase III duration (100DM); weeks | Constant | 8 |  |  |  |  |
| $V_{L; 0}$ | Cost of vaccine delivery at start up (0–10%) in LIC; USD per dose | Triangular | 1 | 1.5 | 2 |  |  |
| $V_{L; 11}$ | Cost of vaccine delivery during ramp up (11–30%) in LIC; USD per dose | Triangular | 0.75 | 1 | 1.5 |  |  |
| $V_{L; 31}$ | Cost of vaccine delivery getting to scale (31–80%) in LIC; USD per dose | Triangular | 1 | 2 | 4 |  |  |
| $V_{LM; 0}$ | Cost of vaccine delivery at start up (0–10%) in LMIC; USD per dose | Triangular | 3 | 4.5 | 6 |  |  |
| $V_{LM; 11}$ | Cost of vaccine delivery during ramp up (11–30%) in LMIC; USD per dose | Triangular | 2.25 | 3 | 4.5 |  |  |
| $V_{LM; 31}$ | Cost of vaccine delivery getting to scale (31–80%) in LMIC; USD per dose | Triangular | 1.5 | 2 | 2.5 |  |  |
| $V_{UM; 0}$ | Cost of vaccine delivery at start up (0–10%) in UMIC; USD per dose | Triangular | 6 | 9 | 12 |  |  |
| $V_{UM; 11}$ | Cost of vaccine delivery during ramp up (11–30%) in UMIC; USD per dose | Triangular | 4.5 | 6 | 9 |  |  |
| $V_{UM; 31}$ | Cost of vaccine delivery getting to scale (31–80%) in UMIC; USD per dose | Triangular | 3 | 4 | 5 |  |  |
| $V_{H; 0}$ | Cost of vaccine delivery at start up (0–10%) in HIC; USD per dose | Triangular | 30 | 40 | 75 |  |  |
| $V_{H; 11}$ | Cost of vaccine delivery during ramp up (11–30%) in HIC; USD per dose | Triangular | 30 | 40 | 75 |  |  |
| $V_{H; 31}$ | Cost of vaccine delivery getting to scale (31–80%) in HIC; USD per dose | Triangular | 30 | 40 | 75 |  |  |
| $M_G$ | Global annual manufacturing volume; billion doses | Constant | 15 |  |  |  |  |
| $M_C$ | Current annual manufacturing volume; billion doses | Constant | 6 |  |  |  |  |
| $F$ | Facility transition start; weeks before vaccine approval | Constant | 7 |  |  |  |  |
| $I_R$ | Weeks to initial manufacturing, reserved infrastructure | Constant | 12 |  |  |  |  |
| $I_E$ | Weeks to initial manufacturing, existing and unreserved infrastructure | Constant | 30 |  |  |  |  |
| $I_B$ | Weeks to initial manufacturing, built and unreserved infrastructure | Constant | 48 |  |  |  |  |
| $C_R$ | Weeks to scale up to full capacity, reserved infrastructure | Constant | 10 |  |  |  |  |
| $C_E$ | Weeks to scale up to full capacity, existing and unreserved infrastructure | Constant | 16 |  |  |  |  |
| $C_B$ | Weeks to scale up to full capacity, built and unreserved infrastructure | Constant | 16 |  |  |  |  |
| $P_0$ | Probability of success; preclinical | Multinomial | 0.41 | 0.57 |  |  |  |
| $P_1$ | Probability of success; Phase I | Multinomial | 0.33 | 0.9 |  |  |  |
| $P_2$ | Probability of success; Phase II | Multinomial | 0.33 | 0.79 |  |  |  |
| $P_3$ | Probability of success; Phase III | Uniform | 0.4 | 0.8 |  |  |  |
| $T_0^{(e)}$ | Cost, preclinical, experienced manufacturer; USD | Exponential | 24213683 | 1700000 | 1.4e+08 |  |  |
| $T_0^{(i)}$ | Cost, preclinical, inexperienced manufacturer; USD | Inverse Gaussian | 7882792 | 13455907 | 1700000 | 3.7e+07 |  |
| $T_1^{(e)}$ | Cost, Phase I, experienced manufacturer; USD | Inverse Gaussian | 15339198 | 8076755 | 1900000 | 7e+07 |  |
| $T_1^{(i)}$ | Cost, Phase I, inexperienced manufacturer; USD | Inverse Gamma | 2.28 | 9799081 | 1e+06 | 3e+07 |  |
| $T_2^{(e)}$ | Cost, Phase II, experienced manufacturer; USD | Log normal | 28297339 | 24061641 | 3800000 | 1.4e+08 |  |
| $T_2^{(i)}$ | Cost, Phase II, inexperienced manufacturer; USD | Inverse Gaussian | 17124622 | 35918793 | 4400000 | 5.4e+07 |  |
| $T_3^{(e)}$ | Cost, Phase III, experienced manufacturer; USD | Inverse Gamma | 1.31 | 51397313 | 1.5e+07 | 9.1e+08 |  |
| $T_3^{(i)}$ | Cost, Phase III, inexperienced manufacturer; USD | Beta prime | 4.89 | 1.69 | 11400026 | 2500000 | 4e+08 |
| $L$ | Licensure; USD | Constant | 287750 |  |  |  |  |
| $T_0^{(B)}$ | BPSV preclinical duration; years | Multinomial | 1 | 2 |  |  |  |
| $T_1^{(B)}$ | BPSV Phase I duration; years | Multinomial | 1 | 2 |  |  |  |
| $T_2^{(B)}$ | BPSV Phase II duration; years | Constant | 2 |  |  |  |  |
| $T_3^{(B)}$ | BPSV Phase III duration; years | Multinomial | 2 | 3 | 4 |  |  |
| $L^{(B)}$ | Licensure duration; years | Constant | 2 |  |  |  |  |
| $G$ | BPSV cost of goods supplied; USD per dose | Constant | 4.68 |  |  |  |  |
| $A$ | Advanced capacity reservation fee; USD per dose per year | Constant | 0.53 |  |  |  |  |
| $S_R$ | SSV procurement price, reserved capacity; USD per dose | Constant | 6.29 |  |  |  |  |
| $S_U$ | SSV procurement price, reactive capacity; USD per dose | Constant | 18.94 |  |  |  |  |
| $E$ | Enabling activities; million USD per year | Constant | 700 |  |  |  |  |
| $I$ | Inflation (2018–2025) | Constant | 0.28 |  |  |  |  |
| $r$ | Discount rate | Uniform | 0.02 | 0.06 |  |  |  |
| $M_p$ | Profit margin | Constant | 0.2 |  |  |  |  |
| $M_f$ | Fill/finish cost | Constant | 0.14 |  |  |  |  |

Notation and parametric assumptions for inputs to the costing model.
Parameters are used as follows: uniform distributions go from Parameter
1 to Parameter 2. Triangular distributions go from Parameter 1 to
Parameter 3 with a peak at Parameter 2. Multinomial distributions have
equally probable values listed individually. Exponential distributions
have as a mean Parameter 1 and are truncated at Parameters 2 and 3.
Inverse Gaussian distributions have as a mean Parameter 1, as a shape
Parameter 2, and are truncated at Parameters 3 and 4. Log normal
distributions have as a mean Parameter 1, as a standard deviation
Parameter 2, and are truncated at Parameters 3 and 4. Inverse Gamma
distributions have shape Parameter 1, scale Parameter 2, and are
truncated at Parameters 3 and 4. Beta Prime distributions have shape
Parameters 1 and 2, scale Parameter 3, and are truncated at Parameters 4
and 5.

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-4.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-5.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-6.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-7.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-8.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-9.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-10.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-11.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-12.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-13.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-14.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-15.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-16.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-17.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-18.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-19.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-20.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-21.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-22.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-23.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-24.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-25.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-26.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-27.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-28.png)<!-- -->

# 3 Model details

## 3.1 Preparedness cost equation (annual calculation, 2025-2039):

(BPSV R&D + BPSV Stockpile + SARS-X Reserved capacity + Enabling
activities) / (1 + discount rate) ^ (year – 2025)

## 3.2 Response cost equation (annual calculation, 2040-2045):

(BPSV R&D + SARS-X R&D + BPSV Procurement + SARS-X Procurement + BPSV
Delivery + SARS-X Delivery) / (1 + discount rate) ^ (year – 2025)

## 3.3 Risk-adjusted R&D cost per candidate calculation -

Sum of the cost of each phase multiplied by the likelihood of phase
occurrence (probability of success for previous phases)

Probability of Occurrence (PoO) = 1 \* PoS (PhaseN-1) …

\$ (Preclin) \* PoO (Preclin) + \$ (Ph1) \* PoO (Ph1) + \$ (Ph2) \* PoO
(Ph2) + \$ (Ph3) \* PoO (Ph3) + \$ (License) \* PoO (License)

Probabilities of success for preclinical, Phase I, Phase II, and Phase
III are $P_0$, $P_1$, $P_2$ and $P_3$. Then probabilities of occurrence
are:

$$\hat{P}_i = \left\{\begin{array}{lr}1 & i=0 \\ \prod_{j=0}^{i-1}P_j & i\in\{1,2,3\} \\ \prod_{j=0}^{3}P_j & i=L \end{array}\right.$$

and the cost of each phase is $T_i$ (a (weighted) sum of experienced and
inexperienced manufacturers?). Then the total cost is

$$D_{\text{R\&D}} = \sum_{i=0}^3 \hat{P}_iT_i + \hat{P}_LL$$

## 3.4 Procurement cost calculation:

Scenario 1: Annual demand under 6.6B

Annual demand \* \$6.29 \* 1.14 \*1.2

Scenario 2: Annual demand over 6.6B

Annual demand \* \$18.94

If we write annual demand as $A_{\cdot,y}$, then we would have

$$D_{\text{SSV},y} = \min\{A_{SSV,y},6.6e9\}\cdot S_R\cdot(1+M_p)\cdot(1+M_f)  + \max\{A_{SSV,y}-6.6e9,0\}\cdot S_U$$

$$D_{\text{BPSV},y} = A_{BPSV,y}\cdot G$$

## 3.5 Delivery Cost Equation:

WB status demand/0.8 \* 0.1 \* (0-10% cost) + WB status demand/0.8 \*
0.2 \* (11-30% cost) + WB status demand/0.8 \* 0.5 \* (30-80% cost)

I assume this means, for total populations $N_i^{(T)}$ in income group
$i\in\{\text{LIC, LMIC, UMIC, HIC}\}$, and delivery cost $D$:

$$D_{\text{SSV}} = \sum_{i}N_i^{(T)}(0.1V_{i; 0} + 0.2V_{i; 11} + 0.5V_{i; 31}) $$

Are the same percentages applied to BPSV, which goes only to people aged
65 or older, with populations $N_i^{(65)}$? Or will coverage be reached
after the “ramp up” phase?

$$D_{\text{BPSV}} = \sum_{i}N_i^{(65)}(0.1V_{i; 0} + 0.2V_{i; 11} + 0.5V_{i; 31}) $$

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

Literature review of global and country-specific delivery rates.

## 3.6 Vaccination Scenarios

Vaccine supply at each investment level is modelled by producing doses
sourced from preexisting advance capacity reservations or from the
private market. The production volumes for each of these supply
mechanisms are interconnected and based on the following assumptions:

- Current annual global vaccine manufacturing capacity is 9B doses and
  is used in the response to the SARS-X pandemic
- A portion of this capacity is reserved in each level for pandemic
  response following the advance capacity reservation volumes in each
  investment level (see Table XX)
- Advance capacity reservations allow manufacturing facilities to
  transition and scale production of SSV doses faster because facilities
  are paid to maintain a high level of transition preparedness and
  produce vaccines on the same platform as the approved SSV
- Advance capacity reservation agreements guarantee a
  population-proportional distribution to all World Bank income groups
  as part of their terms and conditions (an exception is made for the
  500m doses available in all levels, which comes from an existing
  advance capacity reservation held by the European Union)
- Private manufacturing will respond to meet any additional gap in SSV
  demand above the volume covered by advance capacity reservations
- Private manufacturing of the SARS-X vaccine will occur in all existing
  vaccine production facilities not under a reserved capacity agreement
- Additional production capacity of 6B doses will be built in response
  to a global pandemic, like during the Covid-19 response
- In BAU and half of scenarios, private manufacturers prioritise
  delivering SARS-X doses to HICs, like during the Covid-19 response.
  Private manufacturers provide population-proportional distribution to
  all income groups, following the same logic as advanced capacity
  reservations, in the remaining scenarios (see Table 2)
- HICs initially receive SARS-X doses quicker than other World Bank
  income groups due to existing advance capacity reservations and
  preferential treatment from private manufacturers, in select
  scenarios. When HICs receive sufficient doses for vaccinating eligible
  populations, the prioritized doses are reallocated to UMICs. This
  process repeats for LLMICs once all UMICs doses are delivered

| Category | Reserved capacity | Private response (existing capacity) | Private response (built capacity) |
|:---|----|:---|:---|
| Annual manufacturing volune | By scenario (0.5–9B) | 9B minus reserved volume | 6B |
| Facility transition start | 7 weeks before vaccine approval | 7 weeks before vaccine approval | 7 weeks before vaccine approval |
| Weeks to initial manufacturing | 12 | 30 | 48 |
| Scale-up weeks to full capacity | 10 | 16 | 16 |

Manufacturing response timeline assumptions

<!-- | Weeks from transition start | 0-11 | 12-21 | 22-29 | 30-45  | 46-47 | 48-63 | 64+ | -->
<!-- |---|---|---|---|---|---|---|---| -->
<!-- | Reserved Capacity (%)  || Scaling from 0-100 | 100 | 100 | 100 | 100 | 100 |  -->
<!-- | Private Capacity (Existing; %)  || | | Scaling from 0-100 | 100 | 100 | 100 |  -->
<!-- | Private Capacity (Response; %)  | | | |  | | Scaling from 0-100 | 100 | -->

| Weeks from transition start | Reserved Capacity (%) | Private Capacity (Existing; %) | Private Capacity (Response; %) |
|:---|:---|:---|:---|
| 0–11 |  |  |  |
| 12–21 | Scaling from 0-100 |  |  |
| 22–29 | 100 |  |  |
| 30–45 | 100 | Scaling from 0-100 |  |
| 46–47 | 100 | 100 |  |
| 48–63 | 100 | 100 | Scaling from 0-100 |
| 64+ | 100 | 100 | 100 |

Vaccine Production Timeline

# 4 Results

# 5 Attributions / Authors

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Banks2023" class="csl-entry">

Banks, C, RD Estanislao, SJ De los Reyes, JE De Guzman, LB
Sumpaico-Tanchanco, B Makani-Lim, R Archer, and L Boonstoppel. 2023.
“The Cost of Delivering COVID-19 Vaccines in the Philippines.” Geneva:
ThinkWell.
<https://thinkwell.global/wp-content/uploads/2024/03/Cost-of-delivering-COVID19-vaccines-in-the-Philippines_final-report_19-Dec-2023.pdf>.

</div>

<div id="ref-Edoka2024" class="csl-entry">

Edoka, Ijeoma, Lineo Marie Matsela, Khumo Modiba, Yolandie Luther,
Sharlene Govender, Thapelo Maotoe, Heena Brahmbhatt, Pedro T. Pisa,
Gesine Meyer-Rath, and Jacqui Miot. 2024. “Costs of the COVID-19
Vaccination Programme: Estimates from the West Rand District of South
Africa, 2021/2022.” *BMC Health Services Research* 24 (1): 857.
<https://doi.org/10.1186/s12913-024-11251-1>.

</div>

<div id="ref-Griffiths2021" class="csl-entry">

Griffiths, Ulla, Alex Adjagba, Marcia Attaran, Raymond Hutubessy,
Nathalie Van De Maele, Karene Yeung, Wei Aun, et al. 2021. “Costs of
Delivering COVID-19 Vaccine in 92 AMC Countries.” February. COVAX
Working Group. <https://www.who.int/publications/i/item/10665337553>.

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
