100-day mission: Model description
================

-   [1 Simulation rules](#1-simulation-rules)
-   [2 Epi model](#2-epi-model)
    -   [2.1 Disease state transitions](#21-disease-state-transitions)
    -   [2.2 Vaccination state
        transitions](#22-vaccination-state-transitions)
    -   [2.3 Contact rates](#23-contact-rates)
        -   [2.3.1 Matrix $A$: community
            contacts](#231-matrix-a-community-contacts)
        -   [2.3.2 Matrix $B$: Worker-to-worker
            contacts](#232-matrix-b-worker-to-worker-contacts)
        -   [2.3.3 Matrix $C$: Consumers-to-worker
            contacts](#233-matrix-c-consumers-to-worker-contacts)
    -   [2.4 Social distancing](#24-social-distancing)
-   [3 Socio-economic costs](#3-socio-economic-costs)
    -   [3.1 Lost lives](#31-lost-lives)
    -   [3.2 Lost economic output](#32-lost-economic-output)
    -   [3.3 Lost education](#33-lost-education)

# 1 Simulation rules

-   Countries are instantiated with two random variables: the response
    time, and their importation time
-   The response time is the time at which the reporting country reports
    having seen X hospital cases, where X is a random number between 1
    and 20
-   The importation time is a random number between 0 and 20 days, where
    0 days would be equivalent to the spillover, or origin, country
-   The simulation starts at the minimum between the response time and
    the importation time
-   At the response time, the BPSV, if present, is given to people aged
    65 and older; testing begins; social distancing begins; economic
    closures, if in use, are implemented
-   At the importation time, five people are moved from compartment S to
    compartment E
-   If closures are being implemented, the rules in Tables
    <a href="#tab:rulesreactive">1.1</a> and
    <a href="#tab:ruleselimination">1.2</a> are followed
-   The SARS-X–specific vaccine is rolled out starting on day 107 or 372
    after the response time, depending on the investment assumption
-   All people aged 15 and over are eligible for vaccination, and we
    assume 80% take it up
-   Distribution rate increases linearly to a maximum of 1% of the
    population per day, at which is stays until 80% coverage is reached
-   When vaccine rollout is complete, closures, testing and social
    distancing end
-   When the doubling time is more than 30 days and there are fewer than
    1,000 people in hospital, the simulation ends.

| From/to            | No closures                                                                                                             | Light closures                                                    | Heavy closures                                            |
|:-------------------|:------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------|:----------------------------------------------------------|
| **No closures**    |                                                                                                                         |                                                                   | t = response time OR Hospital occupancy &gt; 95% capacity |
| **Light closures** | (Growth rate &lt; 0.025 OR Hospital occupancy &lt; 25% capacity) AND vaccine rollout complete OR $R(D(\textbf{1}) < 1$) |                                                                   | Hospital occupancy &gt; 95% capacity                      |
| **Heavy closures** |                                                                                                                         | Hospital occupancy &lt; 25% capacity AND t &gt; 7 + response time |                                                           |

Table 1.1: State transition rules for reactive closure strategies

| From/to            | No closures                                        | Light closures                            | Heavy closures                                            |
|:-------------------|:---------------------------------------------------|:------------------------------------------|:----------------------------------------------------------|
| **No closures**    |                                                    |                                           | t = response time OR Hospital occupancy &gt; 95% capacity |
| **Light closures** | Vaccine rollout complete OR $R(D(\textbf{1}) < 1$) |                                           | Rt &gt; 1.2                                               |
| **Heavy closures** | Vaccine rollout complete OR $R(D(\textbf{1}) < 1$) | Rt &lt; 0.95 AND t &gt; 7 + response time |                                                           |

Table 1.2: State transition rules for the elimination strategy

# 2 Epi model

## 2.1 Disease state transitions

<div class="figure">

<img src="README_files/figure-gfm/statetransitions-1.png" alt="Disease state transitions"  />

<p class="caption">

Figure 2.1: Disease state transitions

</p>

</div>

Possible transitions between disease states are shown in Figure
<a href="#fig:statetransitions">2.1</a>. Transition rates are functions
of time $t$, vaccination status $v$, and group identity $g$ (where the
groups are the 45 sectors and the four age groups).

The rate of infection of susceptible individuals, $k_1(v,t)$, is defined
as

``` math
k_0(t) = \rho(t)\beta\left(D(x)\cdot I^{(eff)}\right), 
```

$$\begin{equation}
k_1(v,t) = \eta_{A,v}k_0(t)
\qquad(2.1)
\end{equation}$$

with

``` math
 I^{(eff)}=\sum_{v=0}^2(\epsilon (1-p_3+\delta p_3)I_{A,v}+(1-p_4+\delta p_4)I_{S,v}). 
```

Here, $\eta_{A,v}$ is the relative probability to be infected given
vaccine status $v$; $\rho(t)$ is the time-dependent modifier of the rate
of infection, $\beta$, which captures the impact of social distancing;
$D(x)$ is the contact matrix between groups and depends on the economic
configuration $x$; $\epsilon$ is the reduction in infectiousness from
asymptomatic relative to symptomatic individuals; $p_3$ and $p_4$ are
the proportions of asymptomatic and symptomatic infectious days,
respectively, spent self isolating; $\delta$ is the infectiousness
remaining given the rate to self isolate; and $I_{\cdot,\cdot}$ is the
vector of number of infectious asymptomatic ($I_{A,\cdot}$) and
symptomatic ($I_{S,\cdot}$) people who are unvaccinated ($I_{\cdot,0}$),
vaccinated with the BPSV ($I_{\cdot,1}$), or vaccinated with the
specific vaccine ($I_{\cdot,2}$).

``` math
 k_2 = (1-p_S)/\sigma 
```

is the rate to asymptomatic infectiousness, where $p_S$ is the
probability to become symptomatic, and $\sigma$ is the expected duration
of the latent period before the onset of infectiousness;

``` math
 k_3 = 1/\gamma_A  
```

is the rate of recovery from asymptomatic infection;

``` math
 k_4 = p_S/ \sigma; 
```

is the rate of symptom onset;

``` math
k_5 =  (1-p_H) / \gamma_I 
```

is the rate of recovery from symptomatic infection, where $p_H$ is the
probability to be hospitalised, and
$\gamma_I = p_H\gamma_H + (1-p_H)\gamma_R$ is the expected time to be in
compartment $I_S$: $\gamma_H$ is the expected duration before
hospitalisation and $\gamma_R$ is the expected duration before recovery.

``` math
p_H=\eta_{H,v}\hat{p}_H
```

is the baseline probability to be hospitalised ($`\hat{p}_H`$) adjusted
by the vaccine effect protecting against hospitalisation
($`\eta_{H,v}`$). Then

``` math
k_6 = p_H/\gamma_I
```

is the rate of hospitalisation following symptomatic infection.

``` math
k_7 = (1-p_D) / \lambda_H
```

is the rate of recovery of hospitalised patients, where
$`p_D=\hat{p}_Df_H(H)`$ is the baseline probability to die given
hospitalisation, adjusted by a factor encoding the increase in fatality
rate as hospital occupancy increases,
$`f_H(H)=\max\{1,1+1.87(H-H_{\text{max}})/H_{\text{max}}\}`$.
$\lambda_H = p_D\lambda_D + (1-p_D)\lambda_R$ is the expected time to be
in compartment $H$: $\lambda_D$ is the expected duration before death
and $\lambda_R$ is the expected duration before recovery. $p_D$ is the
probability to die given hospitalisation. Finally,

``` math
k_8 = p_D/\lambda_H
```

is the rate of death following hospitalisation.

## 2.2 Vaccination state transitions

In our model, $v=0$ refers to unvaccinated people, $v=1$ to people who
have received a full schedule of BPSV, and $v=2$ to people who have
received a full schedule of the specific vaccine. How we model
transitions between vaccination states is shown in Figure
<a href="#fig:vaccinetransitions">2.2</a>.

$k_9$ and $k_{15}$ represent the rates of BPSV vaccination of
unvaccinated susceptible and recovered people, and $k_{16}$ and $k_{17}$
represent the rates of vaccinating BPSV-vaccinated susceptible and
recovered people. $k_{13}$ and $k_{14}$ represent the rates of
vaccinating people directly with the specific vaccine. $k_{10}$ and
$k_{18}$ are the rates of seroconversion to vaccine-induced immunity,
and $k_{12}(t)=k_1(0,t)$ and $k_{19}(t)=k_1(1,t)$ are the rates of
infection of just-vaccinated people, which returns them to the
epidemiological pathway of the lower vaccination level.

<div class="figure">

<img src="README_files/figure-gfm/vaccinetransitions-1.png" alt="Vaccine state transitions"  />

<p class="caption">

Figure 2.2: Vaccine state transitions

</p>

</div>

## 2.3 Contact rates

The configuration $x$ and the proportion of workers working from home
$q$ determine the scaling of exposure to infection between different
groups for different reasons:

-   Worker absence due to sector closure
-   Worker absence due to working from home
-   Student absence due to school closure
-   Customer absence due to sector closure: impact on workers
-   Customer absence due to sector closure: impact on customers

We construct contact matrix $D(x)$ as the sum of four matrices: $A(x)$
(community contacts), $B(x)$ (worker-to-worker contacts), $C(x)$
(consumer-to-worker contacts), and $\hat{C}(x)$ (worker-to-consumer
contacts). We construct peacetime matrices ($x=\textbf{1}$) beginning
with a “target matrix,” which the four matrices should add up to, which
is taken from (Walker et al. 2020). By sampling relevant values, we
decompose the whole matrix into its component parts. To incorporate
closures, each matrix is transformed independently, before they are all
added together again.

Matrix $D(\textbf{1})$ is estimated using as a basis a contact matrix
from (Walker et al. 2020). These are 16-by-16 matrices, ($D^{(16)}$),
for five-year age bands $a$ up to age group 75+. We map the matrix to a
four-by-four matrix $D^{(4)}$ corresponding to the four age groups $g$
used in the DAEDALUS model, using population sizes, $\hat{P}_a$:

``` math
D_{gg'}^{(4)} = \frac{\sum_{a\in g}\hat{P}_{a}\sum_{a'\in g'}D^{(16)}_{a,a'}}{\sum_{a\in g}\hat{P}_{a}}.
```

Using $P_g$ to represent the population sizes of the DAEDALUS age
groups,

``` math
P_g=\sum_{a\in g}\hat{P}_a,
```

We get to the matrix $D(\textbf{1})$ by broadcasting the four-by-four
matrix to the 49-by-49 one. Contacts from all groups $i$ to working
groups $j$ depend on the age group of the group ($`g(i)`$), and the
fraction of the age-population represented in group $j$, where $w_{j}$
is the number of people in group $j$:

``` math
D_{ij}(\textbf{1}) = D^{(4)}_{g(i),g(j)}\frac{w_{j}}{P_{g(j)}}
```

for $i$ and $j$ including all groups (working and non-working). Each
group $i$ contains people that belong to only one age group $g$. We
refer to the age group of the people in group $i$ as $g(i)$. Then
$P_{g(j)}$ is the number of people in the age group of group $j$, so
$P_{g(j)}=w_{j}$ for age groups 0 to 4, 5 to 19 and 65+, and
$P_{g(j)}=\sum_{j\in\{1,...,N,N+3\}}w_{j}$ for ages 20 to 64.

In setting up a country, we sample values for $D^{(16)}$ (from which we
get $`D(\textbf{1})`$). At the same time, we sample the proportion of
contacts that come from workplaces, and workplace-related contacts. From
these, we get $B(\textbf{1})$ and $C(\textbf{1})$, constructing the
matrices and normalising.

Matrix B is diagonal owing to lack of data regarding between-sector
contacts (Haw et al. 2022). Note that $B_{ii}(\textbf{1})=0$ for $i>N$.
Consumer-to-worker contacts (matrix $C$) describe contacts experienced
by workers from consumers per sector. Note that $C_{ij}(\textbf{1})=0$
for $i>N$. Matrix $\hat{C}(\textbf{1})$ is the complement of matrix
$C(\textbf{1})$, computed by multiplying through by population,
transposing, and dividing again by population.

With $D(\textbf{1})$, $C(\textbf{1})$, $B(\textbf{1})$ and
$\hat{C}(\textbf{1})$, we learn $A(\textbf{1})$.

$A$ is decomposed into its constituent parts, representing intra- and
inter-household interactions ($L$), school interactions ($S$),
hospitality interactions ($H$) and travel interactions ($T$):

``` math
A(\textbf{1})=A^{(L)} + A^{(S)}(\textbf{1}) + A^{(H)}(\textbf{1}) + A^{(T)}(\textbf{1})
```

Values for $A^{(S)}(\textbf{1})$ come from sampled values representing
the fractions of contacts that come from school. School contacts are
estimated separately in two age groups (pre-school age: 0—4; school age:
5—19): $A^{(S)}(\textbf{1})$ has entries of zero for groups $g$ not in
school, and values for $g$=0 to 4 years old and $g$=5 to 19 year olds.

Likewise, $A^{(T)}(\textbf{1})$ is also sampled as a fraction of total
contacts. $A_{ij}^{(T)}(\textbf{1})\geq 0$ for $i=1,...,N$.
$A_{ij}^{(T)}(\textbf{1})=0$ for $i>N$.

Finally, $A^{(H)}(\textbf{1})$ is sampled as a fraction of
$A(\textbf{1})- A^{(S)}(\textbf{1}) - A^{(T)}(\textbf{1})$, which leaves
$A^{(L)}$.

### 2.3.1 Matrix $A$: community contacts

-   Any contact made at home, in a vehicle or other private place,
    retail outlet, public transport, leisure facilities, with loved ones
    in a closed place (“Chez des proches en lieux clos”), open place
    (park, street)
-   Disaggregated by age group (0 – 4; 5 – 19; 20 – 64; 65+)
-   Consumer-consumer contacts are added with respect to the proportion
    the hospitality and education sectors are opened

We construct $A(x)$ from its constituent parts, representing intra- and
inter-household interactions ($L$), school interactions ($S$),
hospitality interactions ($H$) and travel interactions ($T$):

``` math
A(x)=A^{(L)} + A^{(S)}(x) + A^{(H)}(x) + A^{(T)}(x).
```

School contacts under $x$ are simply the scaled values. $x_{S}$ is the
extent to which schools are open, so that the number of contacts per
person scales superlinearly with school closure.

$$\begin{equation}
A_{ii}^{(S)}(x)=x_{S}^2A_{ii}^{(S)}(\textbf{1}).
\qquad(2.2)
\end{equation}$$

Matrix $A^{(T)}$ counts contacts between working people, representing
travel. We assume that transport contacts only add to the infection risk
if the sector is open and the workers travel to and from their
workplace. Again, the value for configuration $x$ is the value for
$\textbf{1}$ scaled accordingly:

$$\begin{equation}
A_{ij}^{(T)}(x) = x_{j}(1-q_i)(1-q_j)A_{ij}^{(T)}(\textbf{1}).
\qquad(2.3)
\end{equation}$$

$q_i$ is the proportion of workers from sector $i$ working from home,
and $(1-q_i)(1-q_j)$ scales contacts between workers superlinearly to
approximate the reduced transmission between commuting workers: there
should be fewer contacts per person on average, and there should be
fewer people having these contacts.

Also in this equation, $x_{j}$ scales the numbers of contacts linearly
with respect to sector closure. At the same time, the number of people
in the compartments will be reduced by their sector closure, $x_{i}$.
This, in combination with the scaled contacts, leads to superlinear
scaling.

Matrix $A^{(H)}(x)$ gives the contacts made in the hospitality sector:

$$\begin{equation}
A^{(H)}(x) = x_{H}^2A^{(H)}(\textbf{1})
\qquad(2.4)
\end{equation}$$

The value $x_{H}$ is the workforce-weighted average extent to which the
hospitality sectors are open, so that the number of contacts per person
scales superlinearly according to closure:

``` math
x_{H} = \frac{\sum_ix_{i}w_i}{\sum_iw_i}
```

where we sum over only the hospitality sectors.

### 2.3.2 Matrix $B$: Worker-to-worker contacts

-   Contacts made at work (office, studio, etc.) and which are reported
    to be made (almost) every day, or a few times per week  
-   Disaggregated by sector
-   Individuals who stated that they are in employment
-   Individuals who are of working age (20 – 64)

$$\begin{equation}
B_{ii}(x) = x_{i}(1-q_i)^2B_{ii}(\textbf{1}),
\qquad(2.5)
\end{equation}$$

for the $i=1,...,N$ working groups, with the number of contacts adjusted
according to at-home working ($q_i$) and sector openness ($x_{i}$). As
before, there is superlinear scaling of contacts with respect to working
from home. There is linear scaling with respect to sector closure: that
is, there are fewer contacts per person, but we do not approximate there
being fewer people having them. This is because the latter is accounted
for in the movement of people out of the group upon its closure.

### 2.3.3 Matrix $C$: Consumers-to-worker contacts

-   Contacts made at work (office, studio, etc.) and which are reported
    to be a few times per month, a few times per year or less often, for
    the first time
-   Disaggregated by sector
-   Individuals who stated that they are in employment  
-   Individuals who are of working age (20 – 64)  
-   If more than 20 contacts are made by the individual, the survey
    respondent could state the total number of contacts made instead of
    listing all individual contacts. If this was the case, this number
    was used instead of the sum of individual contacts made

$$\begin{equation}
C_{ij}(x) = x_{i}(1-q_i)C_{ij}(\textbf{1}),
\qquad(2.6)
\end{equation}$$

for $j=1,...,N+3$.

Here, there is linear scaling of $C_{ij}(\textbf{1})$ with respect to
working from home, and linear scaling with respect to sector closure,
which becomes superlinear scaling for sectors as individuals are moved
out of the compartment, as with matrix $B(x)$.

## 2.4 Social distancing

We parametrise the effects of ‘social distancing’ in the model using
Google’s mobility data (Figure <a href="#fig:smoothmobility">2.3</a>).
These changes in mobility were consequences of both government mandates
and individual’s choices. As we cannot separate the two, we consider a
range of possibilities, based on the range of mobility changes observed
for a given level of stringency (Figure
<a href="#fig:mobilitydrop">2.4</a>). In our model, the mandated
economic configuration leads to a change in contacts. We associate the
reduction in contacts, which translates as a relative reduction in
transmission, with the reduction in mobility.

<div class="figure">

<img src="README_files/figure-gfm/smoothmobility.png" alt="Mobility trajectories in 2020 for all countries, with points showing the point at which the largest drop was observed. Trajectories are averaged over &quot;Retail and recreation&quot;, &quot;Transit stations&quot; and &quot;Workplaces&quot; and smoothed with a spline of 80 knots." width="2809" />

<p class="caption">

Figure 2.3: Mobility trajectories in 2020 for all countries, with points
showing the point at which the largest drop was observed. Trajectories
are averaged over “Retail and recreation,” “Transit stations” and
“Workplaces” and smoothed with a spline of 80 knots.

</p>

</div>

<div class="figure">

<img src="README_files/figure-gfm/mobilitydrop.png" alt="The largest drop in mobility from Figure  ef{fig:smoothmobility} is plotted against the stringency on that date." width="3000" />

<p class="caption">

Figure 2.4: The largest drop in mobility from Figure
ef{fig:smoothmobility} is plotted against the stringency on that date.

</p>

</div>

-   We want to write mobility as a function of mandate and some epi
    outcome, e.g. deaths: $m = (1-b)f(y,g) + b$ where $m$ is mobility,
    $y$ is deaths per million, $g$ is government mandate, and
    $`0 < b < 1`$ is the baseline.
-   We want mobility to drop monotonically with both the mandate and the
    epi outcome: $\frac{df}{dy}<0$, $\frac{df}{dg}<0$.
-   We want a maximum mobility of 1 when both the mandate and the epi
    outcome are 0: $f(0,0)=1$.
-   We want mobility to approach $b$ when the mandate and the epi
    outcome become large: $\lim_{x\to 10^6, g\to 1}f(y,g)= 0$.
-   We want to allow for the possibility of redundancy between the two
    variables: $f(0,0)/f(0,g) > f(x,0)/f(y,g)$ and
    $f(0,0)/f(y,0) > f(0,g)/f(y,g)$ for $y,g>0$.

A simple model to achieve these criteria is:
$$f(y,g) = \frac{1}{1+k_1y+k_2g}$$ with $k_1, k_2>0$.

However, we might also want a model that can be parametrised with a
distribution whose uncertainty covers the whole range of possible
eventualities. The equivalent model with compounded effects would be
$$f_1(y,g) = \frac{1}{1+k_1y}\frac{1}{1+k_2g}.$$ The equivalent model
with completely overlapping effects would be
$$f_2(y,g) = \frac{1}{1+\max(k_1y,k_2g)}.$$ Then we could include ‘model
uncertainty’ via some parameter $\beta\sim\mathcal{U}(0,1)$, defining
$$f(y,g) = (f_1(y,g))^\beta(f_2(y,g))^{(1-\beta)}.$$

<div class="figure">

<img src="README_files/figure-gfm/mobilityfitted.png" alt="Fit of model to data." width="2096" />

<p class="caption">

Figure 2.5: Fit of model to data.

</p>

</div>

<div class="figure">

<img src="README_files/figure-gfm/mobilityposterior.png" alt="Posterior distribution for parameters $k_1$ and $b$." width="2096" />

<p class="caption">

Figure 2.6: Posterior distribution for parameters $k_1$ and $b$.

</p>

</div>

<div class="figure">

<img src="README_files/figure-gfm/mobilitycurves.png" alt="Sampled curves for four levels of mitigation. Data shown as points." width="2096" />

<p class="caption">

Figure 2.7: Sampled curves for four levels of mitigation. Data shown as
points.

</p>

</div>

# 3 Socio-economic costs

We assign monetary values to YLLs and to years of education in order to
add health and education costs of mitigation strategies to the costs of
economic closures. We define the total socio-economic costs TSC of an
epidemic as the sum of the individual costs:

$$\begin{equation}
\text{TSC} = f_{\text{education}}(X,\text{VSY}) + f_{\text{lives}}(D,l,\text{VLY},r) + f_{\text{GDP}}(Z_0-Z,r),
\label{eq:swf}
\end{equation}$$

where the arguments are the number of school years lost ($X$) and the
value of a school year (VSY); the number of deaths per age group due to
COVID-19 ($D$), the remaining life expectancy per age group ($l$), the
value of a life year (VLY), and discount rate $r$; and the lost output
over the period due to reduced economic activity ($Z_0-Z$).

## 3.1 Lost lives

To value lives lost, we make use of the expected remaining life years
per age group (Global Burden of Disease Collaborative Network 2021).
These are used to estimate the expected number of years of life lost per
death, and to estimate the value of a life year. We map the remaining
life expectancy $l_a$ for the GBD age groups $a$ to $l_g$ for the model
age groups $g$ as a population-weighted average, taking into account the
size of each age group, $N_a$. For the expected number of life years
lost per death, we take into account also the probability to die given
infection, $P(D|I,a)$:
$$l_g^{\text{(death)}} = \frac{\sum_{a\in g}N_al_aP(D|I,a)}{\sum_{a\in g}N_aP(D|I,a)}; $$
$$l_g^{\text{(life)}} = \frac{\sum_{a\in g}N_al_a}{\sum_{a\in g}N_a}; $$

Expected life years remaining with discounting taken into account can be
written

$$\begin{equation}
\hat{l}_g=\sum_{y=1}^{l_g}\frac{1}{(1+r)^{y}}
\end{equation}$$

for discount rate $r>0$. The discounted number of years lost given $D_g$
deaths due to COVID-19 for each age group is $$\begin{equation}
Y=\sum_gD_g\hat{l}_g^{\text{(death)}}.
\end{equation}$$

The VLY used by policy makers should reflect the value that members of
the society place on reductions of their own mortality. We rely on the
intrinsic rather than instrumental interpretation of the valuation of
life (Cutler and Summers 2020), and we use existing estimates of the
value of a statistical life (VSL) to estimate VLY. We interpret the VSL
as a population-weighted average (Ananthapavan et al. 2021; Robinson,
Sullivan, and Shogren 2021), where each age group has a VSL defined by
the number of expected life years remaining, and where each discounted
year has the same value:

$$\begin{equation}
\text{VSL}=\frac{\sum_gN_g\hat{l}_g^{\text{(life)}}}{\sum_gN_g}\text{VLY}.
\end{equation}$$

Finally, we define the value of lives lost as

$$\begin{equation}
f_{\text{lives}}(D,l^{\text{(death)}},\text{VLY},r)=\text{VLY}*Y.
\end{equation}$$

## 3.2 Lost economic output

We measure the cost of economic closures in terms of lost gross value
added (GVA): the GDP generated by an economic configuration is the
maximum GVA (denoted $z_i$ for each sector $i$) multiplied by the
respective sector openings, summed over the period. The maximum possible
GDP (which is with no closures) is $Z_0=\frac{T}{365}\sum_{i}z_i$, and
we use pre-pandemic output to define the maximum possible values.

$x_{i}(d)$ is the proportion of the workforce contributing to economic
production out of the total workforce $N_i$ on day $d$. The workforce
can be additionally depleted due to sickness, hospitalisation and death,
leaving a smaller fraction ($`\hat{x}_{i}(d)`$) to contribute to
production.

$$\begin{equation}
\hat{x}_{i}(d)=x_{i}(d)-\sum_{v}\left(I_{S,v,i}(d)+H_{v,i}(d)+D_{v,i}(d)\right)/N_i
\end{equation}$$

for vaccine status $v$, infectious and symptomatic $I_S$, hospitalised
$H$, deceased $D$, with total population $N_i$ for sector $i$. Then the
total GVA for the $T$-day period is:

$$\begin{equation}
Z = \frac{1}{365}\left(\sum_{i\neq\text{ed}}^{\mathcal{N}}\sum_{d=1}^Tz_i\hat{x}_{i}(d) + Tz_{\text{ed}}\right)
\qquad(3.1)
\end{equation}$$

and the GDP loss compared to the maximum is $Z_0-Z$. All economic
sectors contribute GVA according to the level they are open for
production, except for the education sector which contributes its
maximum possible monthly GVA, $z_{\text{ed}}$ per month.

## 3.3 Lost education

Education loss is quantified as the effective number of years of
schooling lost. In each month $\tau=1,...,T$ of the projection period,
education is divided into a fraction of in-person teaching
$(w_{\text{ed},\tau})$ and remote teaching $(1-w_{\text{ed},\tau})$. We
compute the total amount of education completed as a weighted sum of the
school opening and its complement, where remote teaching contributes a
fraction $s$ of the educational value of in-person teaching. The total
educational loss is the amount of teaching that was remote
$(1-w_{\text{ed},\tau})$ multiplied by the value lost in remote teaching
$(1-s)$, multiplied by the number of people of school age and the number
of years that each period $\tau$ represents:
`math\label{eq:ed_cost} X=\frac{1}{12}\sum_{\tau=1}^T(1-s)(1-w_{\text{ed},\tau})\sum_{g\in\text{school age}}N_g,`
where $N_g$ is the number of people in age group $g$. We consider only
the second age group out of four (5 to 19 years) to be school-age
students.

The VSY is the value one year of education adds to the economy over a
person’s lifetime. This value can be estimated using existing
projections that the loss of an entire school year will cost 202% of
GDP, assuming that growth will be impacted by having a less productive
workforce for the following 52 years (a cohort of 12 years with careers
expected to last 40 years) . Finally,
`mathf_{\text{education}}(X,\text{VSY}) = X*\text{VSY}.`

For the effectiveness of remote teaching, we use $s=\frac{1}{3}$ , where
infrastructure coverage, household access and effectiveness of teaching
were considered. %, based on an estimated 96% of areas having the
necessary infrastructure, access 86% of students being able to access
the material, and effectiveness of teaching of 40% . This aligns with an
analogous estimate for the Philippines of 37% .

Using the projections that the loss of a year of education for the
current school cohort is worth 202% of current GDP , we obtain a value
of \$34,000 for a person–year of schooling in Indonesia. We contrast
this with a projection of educational losses in the Philippines , for
whom one year of remote teaching (at 37% the effectiveness of in-person
teaching) was estimated to have cost the economy 10.7 trillion pesos,
which equates one full year lost to 89% of the country’s GDP in 2019.
Applied to Indonesia, the value of a year of education per student would
be \$15,000. These percentages align well with the estimate that for
middle-income countries a whole year of lost education costs 73% of GDP,
if we assume that the closures affect 90% of students. There, the
authors estimate only the cost to lifetime wages without considering
consequent changes to demand, consumption or growth.

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Ananthapavan2021" class="csl-entry">

Ananthapavan, Jaithri, Marj Moodie, Andrew J. Milat, and Rob Carter.
2021. “<span class="nocase">Systematic review to update ‘value of a
statistical life’ estimates for Australia</span>.” *International
Journal of Environmental Research and Public Health* 18 (11).
<https://doi.org/10.3390/ijerph18116168>.

</div>

<div id="ref-Cutler2020" class="csl-entry">

Cutler, David M., and Lawrence H. Summers. 2020. “<span
class="nocase">The COVID-19 pandemic and the \$16 trillion
virus</span>.” *JAMA* 324 (15). <https://doi.org/10.1257/pol.20170046>.

</div>

<div id="ref-GlobalBurdenofDiseaseCollaborativeNetwork2021"
class="csl-entry">

Global Burden of Disease Collaborative Network. 2021. “<span
class="nocase">Global Burden of Disease Study 2019 (GBD 2019) Reference
Life Table</span>.” Seattle, United States of America: Institute for
Health Metrics; Evaluation (IHME).

</div>

<div id="ref-Haw2020" class="csl-entry">

Haw, David, Giovanni Forchini, Patrick Doohan, Paula Christen, Matteo
Pianella, Rob Johnson, Sumali Bajaj, et al. 2022. “<span
class="nocase">Optimizing social and economic activity while containing
SARS-CoV-2 transmission using DAEDALUS</span>.” *Nature Computational
Science* 2: 223–33. <https://doi.org/10.25561/83928>.

</div>

<div id="ref-Robinson2021" class="csl-entry">

Robinson, Lisa A., Ryan Sullivan, and Jason F. Shogren. 2021. “<span
class="nocase">Do the benefits of COVID-19 policies exceed the costs?
Exploring uncertainties in the age–VSL relationship</span>.” *Risk
Analysis* 41 (5): 761–70. <https://doi.org/10.1111/risa.13561>.

</div>

<div id="ref-Walker2020" class="csl-entry">

Walker, Patrick G. T., Charles Whittaker, Oliver J. Watson, Marc
Baguelin, Peter Winskill, Arran Hamlet, Bimandra A. Djafaara, et al.
2020. “<span class="nocase">The impact of COVID-19 and strategies for
mitigation and suppression in low- and middle-income countries</span>.”
*Science* 369 (6502): 413–22. <https://doi.org/10.1126/science.abc0035>.

</div>

</div>
