100-day mission: Model description
================

-   [1 Disease state transitions](#1-disease-state-transitions)
-   [2 Vaccination state transitions](#2-vaccination-state-transitions)
-   [3 Contact rates](#3-contact-rates)
    -   [3.1 Matrix $A$: community
        contacts](#31-matrix-a-community-contacts)
    -   [3.2 Matrix $B$: Worker-to-worker
        contacts](#32-matrix-b-worker-to-worker-contacts)
    -   [3.3 Matrix $C$: Consumers-to-worker
        contacts](#33-matrix-c-consumers-to-worker-contacts)
    -   [3.4 Reduced contact rates as a function of economic
        closures](#34-reduced-contact-rates-as-a-function-of-economic-closures)
        -   [3.4.1 Method 1](#341-method-1)
        -   [3.4.2 Method 2](#342-method-2)
        -   [3.4.3 Method 3](#343-method-3)
-   [4 Definitions of socio-economic
    costs](#4-definitions-of-socio-economic-costs)
    -   [4.0.1 Valuing lost lives](#401-valuing-lost-lives)
    -   [4.0.2 Valuing lost education](#402-valuing-lost-education)

# 1 Disease state transitions

<div class="figure">

<img src="README_files/figure-gfm/statetransitions-1.png" alt="Disease state transitions"  />

<p class="caption">

Figure 1.1: Disease state transitions

</p>

</div>

Possible transitions between disease states are shown in Figure
<a href="#fig:statetransitions">1.1</a>. Transition rates are functions
of time $t$, vaccination status $v$, and group identity $g$ (where the
groups are the 45 sectors and the four age groups).

The rate of infection of susceptible individuals, $k_1(v,t)$, is defined
as

``` math
k_0(t) = \rho(t)\beta\left(D(x)\cdot I^{(eff)}\right), 
```

$$\begin{equation}
k_1(v,t) = \eta_{A,v}k_0(t)
\qquad(1.1)
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
p_H=\eta_{H,v}f_H(H)\hat{p}_H
```

is the baseline probability to be hospitalised ($`\hat{p}_H`$) adjusted
by the vaccine effect protecting against hospitalisation
($`\eta_{H,v}`$) and a function encoding the increase in fatality rate
as hospital occupancy increases, $`f_H(H)`$. Then

``` math
k_6 = p_H/\gamma_I
```

is the rate of hospitalisation following symptomatic infection.

``` math
k_7 = (1-p_D) / \lambda_H
```

is the rate of recovery of hospitalised patients, where $p_D$ is the
probability to die given hospitalisation, and
$\lambda_H = p_D\lambda_D + (1-p_D)\lambda_R$ is the expected time to be
in compartment $H$: $\lambda_D$ is the expected duration before death
and $\lambda_R$ is the expected duration before recovery. $p_D$ is the
probability to die given hospitalisation. Finally,

``` math
k_8 = p_D/\lambda_H
```

is the rate of death following hospitalisation.

# 2 Vaccination state transitions

In our model, $v=0$ refers to unvaccinated people, $v=1$ to people who
have received a full schedule of BPSV, and $v=2$ to people who have
received a full schedule of the specific vaccine. How we model
transitions between vaccination states is shown in Figure
<a href="#fig:vaccinetransitions">2.1</a>.

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

Figure 2.1: Vaccine state transitions

</p>

</div>

# 3 Contact rates

In the DAEDALUS model (Haw et al. 2022), the economic configuration is
manifest in the epidemiological model via changes of compartment
membership, and scaling of contact rates. The configuration $x_{i}$ and
the proportion of workers working from home $p_i$ determine the scaling
of contacts between workers, between consumers and workers, between
school students, and between consumers in commercial settings.

We construct contact matrix $D(x)$ as the sum of three matrices: $A(x)$
(community contacts), $B(x)$ (worker-to-worker contacts), $C(x)$
(consumer-to-worker contacts), and $\hat{C}(x)$ (worker-to-consumer
contacts). We construct peacetime matrices ($x=\textbf{1}$) beginning
with a “target matrix,” which the four matrices should add up to, which
is taken from (Walker et al. 2020). By sampling relevant values, we
decompose the whole matrix into its component parts. To incorporate
closures, each matrix is transformed independently, before they are all
added together again.

Matrix $D(\textbf{1})$ is estimated using as a basis the contact matrix
from (Walker et al. 2020). This is a 16-by-16 matrix, ($D^{(16)}$), for
five-year age bands $a$ up to age group 75+. We map it to a four-by-four
matrix $D^{(4)}$ corresponding to the four age groups $g$ used in the
DAEDALUS model, using population sizes, $\hat{P}_a$:

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
$P_{g(j)}=\sum_{j\in\{1,...,N,N+4\}}w_{j}$ for ages 20 to 64.

## 3.1 Matrix $A$: community contacts

-   Any contact made at home, in a vehicle or other private place,
    retail outlet, public transport, leisure facilities, with loved ones
    in a closed place (“Chez des proches en lieux clos”), open place
    (park, street)
-   Disaggregated by age group (0 – 4; 5 – 19; 20 – 64; 65+)
-   Consumer-consumer contacts are added with respect to the proportion
    the hospitality and education sectors are opened

Community contacts (matrix $A$) are any contacts made by people who are
not at work. This includes contacts in the household, during travel to
and from the workplace and non-work-related travel, outside spaces,
leisure activities (e.g. meeting friends), retail outlets
(e.g. supermarkets), and contacts made in the hospitality or service
sectors. The columns of the community matrix $A$ are weighted by the
size of the workforce (measured in headcounts) in each sector. The value
of row sums depends on the extent to which given sectors are open. If
all economic sectors were fully closed, the total contact is reduced to
roughly 40% of total contact when all sectors are open.

We decompose Matrix $A$ into its constituent parts, representing intra-
and inter-household interactions ($L$), school interactions ($S$),
hospitality interactions ($H$) and travel interactions ($T$):

``` math
A=A^{(L)} + A^{(S)} + A^{(H)} + A^{(T)}
```

School contacts are estimated separately in two age groups (pre-school
age: 0 – 4; school age: 5 – 19). Diagonal matrix $A^{(S)}$ counts the
contacts in schools. It has entries of zero for groups $g$ not in
school, and sampled values for $g$=0 to 4 years old and $g$=5 to 19 year
olds. Then

$$\begin{equation}
A_{ii}^{(S)}=x_{S}^2s_{g(i)}.
\qquad(3.1)
\end{equation}$$

The value $x_{S}$ is the extent to which schools are open, so that the
number of contacts per person scales superlinearly according to closure.
Values for $A_{ii}$ are taken as random variables ranging from 0 to
$D_{ii}$.

Matrix $A^{(T)}$ counts contacts between working people, representing
travel. We assume that transport contacts only add to the infection risk
if the sector is open and the workers travel to and from their
workplace.

$$\begin{equation}
A_{ij}^{(T)} = x_{j\tau}(1-p_{i\tau})(1-p_{j\tau})\frac{2.5w_{j}}{\sum_gP_g}
\qquad(3.2)
\end{equation}$$

for $i=1,...,N$. $A_{ij}^{(T)}=0$ for $i>N$.

$p_{i\tau}$ is the proportion of workers from sector $i$ working from
home during period $\tau$, and $(1-p_{i\tau})(1-p_{j\tau})$ scales
contacts between workers superlinearly to approximate the reduced
transmission between commuting workers: there should be fewer contacts
per person on average, and there should be fewer people having these
contacts. We reduce the transmission rates within the groups as a proxy
for moving the individuals out of the group.

Also in this equation, $x_{j\tau}$ scales the numbers of contacts
linearly with respect to sector closure. At the same time, the number of
people in the compartments will be reduced by their sector closure,
$x_{i\tau}$. This, in combination with the scaled contacts, leads to
superlinear scaling.

Matrix $A^{(H)}$ gives the contacts made in the hospitality sector. Each
age group makes an average of 0, 0.5, 1 and 1.5 total contacts for age
groups 0-4, 5-19, 20-64, and 65+, respectively . These contacts are made
in proportion to population, so we can write

$$\begin{equation}
A_{ij}^{(H)} = x_{H,\tau}\frac{A_{g(i)}^{(H0)}w_{j}}{\sum_{j'}w_{j'}}
\qquad(3.3)
\end{equation}$$

with $A^{(H0)} = \{0,0.5,1,1.5\}$, and $A_{ij}^{(H)}=0$ for
$g(i)\neq g(j)$.

The value $x_{H,\tau}$ is the workforce-weighted average extent to which
the hospitality sectors are open in the period $\tau$, so that the
number of contacts per person scales linearly according to closure:

``` math
x_{H,\tau} = \frac{\sum_ix_{i\tau}w_i}{\sum_iw_i}
```

where we sum over only the hospitality sectors.

Finally, $A^{(L)} = A - (A^{(S)} + A^{(H)} + A^{(T)}.)$

## 3.2 Matrix $B$: Worker-to-worker contacts

-   Contacts made at work (office, studio, etc.) and which are reported
    to be made (almost) every day, or a few times per week  
-   Disaggregated by sector
-   Individuals who stated that they are in employment
-   Individuals who are of working age (20 – 64)

Worker-to-worker contacts (matrix B) describe the at-work contacts in
sectors, i.e. the number of contacts per day reported by an individual
actively working in the same sector, which we denote $b_i$ (values given
per sector in Supplementary Table 1). Here “actively working” refers to
one period, i.e. two months in our application. Matrix B is diagonal
owing to lack of data regarding between-sector contacts. Worker-worker
contacts are defined by those contacts recorded to have happened at work
and frequently (reported as a contact made almost every day). At work
contacts at low frequency are classified as worker-consumer contacts.
%At-home working ($p_{i\tau}$) is considered, and community contact
rates apply for contacts between working household members.

Then

``` math
B_{ii} = x_{i\tau}(1-p_{i\tau})^2b_i,
(\#eq:worker)
```

for the $i=1,...,N$ working groups, with the number of contacts adjusted
according to at-home working ($p_{i\tau}$) and sector openness
($x_{i\tau}$). Note that $B_{ii}=0$ for $i>N$. As before, there is
superlinear scaling of contacts with respect to working from home. There
is linear scaling with respect to sector closure: that is, there are
fewer contacts per person, but we do not approximate there being fewer
people having them. This is because the latter is accounted for in the
movement of people out of the group upon its closure.

## 3.3 Matrix $C$: Consumers-to-worker contacts

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
-   Matrix $\hat{C}$ is the complement of matrix $C$, computed by
    multiplying through by population, transposing, and dividing again
    by population.

Consumer-to-worker contacts (matrix $C$) describe contacts experienced
by workers from consumers per sector, denoted $c_i$. As for $B$, the
columns are weighted by sector population, though the row sums are
sector specific. Contacts experienced by workers from consumers are
defined by those contacts recorded to have happened at work less
frequently than every day (i.e. recorded as a few times a week, a few
times a month, a few times a year or less often, or for the first time).
%For a technical description as to how the contact matrices are
constructed, please see the Supplementary Material.

Then

``` math
C_{ij} = x_{i\tau}(1-p_{i\tau})\frac{c_{i}w_{j}}{\sum_{j'}^{N+4}w_{j'}},
(\#eq:ctow)
```

for $j=1,...,N+4$. $C_{ij}=0$ for $i>N$.

Here, there is linear scaling of $C_{ij}$ with respect to working from
home, and linear scaling with respect to sector closure, which becomes
superlinear scaling for sectors as individuals are moved out of the
compartment, as with matrix $B$.

## 3.4 Reduced contact rates as a function of economic closures

There are three methods we use to model reduction in transmission due to
sector closures:

1.  Moving individuals into a different compartment, and scaling the
    number of contacts linearly
2.  Scaling contacts superlinearly only
3.  Scaling contacts linearly only

and there are five scenarios where transmission is reduced:

-   Worker absence due to sector closure (uses method 1)
-   Worker absence due to working from home (uses method 2)
-   Student absence due to school closure (uses method 3)
-   Customer absence due to sector closure: impact on workers (uses
    method 1)
-   Customer absence due to sector closure: impact on customers (uses
    method 3)

### 3.4.1 Method 1

Recall that the number of infections is proportional to susceptible
times contact rate times proportion of population infectious
($S_iD_{ij}I_j/P_j$). When the number of people is reduced, the number
of infections is reduced via both the number of people at risk ($S_i$)
and the number of contacts per person ($D_{ij}$). This is precisely how
we model worker absence due to sector closure: we reduce the number of
contacts per person by scaling the contact rates according to the
fraction who remain in the workplace, and the number of people who are
susceptible is reduced by moving workers out of their sector compartment
following to the closure schedule. The result is superlinear scaling
with respect to closure. The same method applies for worker contacts
made with the general population.

### 3.4.2 Method 2

For the other type of worker absence, working from home, we do not move
any people from their compartmental group, so we compensate by reducing
the number of contacts further by the amount the compartment would have
reduced to, i.e. by the proportion of workers who remain,
$(1-p_{i\tau})$, if $p_{i\tau}$ is the fraction of workers from sector
$i$ who work from home in time period $\tau$. The same method applies
for worker contacts while travelling.

### 3.4.3 Method 3

When workers work from home, each member of the population should
interact with fewer workers each, which we approximate by scaling the
number of contacts made between workers and consumers down by
($1-p_{i\tau}$). This scaling is linear, reflecting that only the number
of contacts per person changes, not the number of people making
contacts. Likewise, the linear scalings reflect fewer contacts per
person in schools and hospitality settings, but we retain the same
numbers of people who make contacts at all.

# 4 Definitions of socio-economic costs

We assign monetary values to YLLs and to years of education in order to
add health and education costs of mitigation strategies to the costs of
economic closures. This enables comparisons between configurations,
which forms the basis of the decision framework.

We define the total socio-economic costs TSC of an economic
configuration $\{x_{i,\tau}\}$ as the sum of the individual costs:
`math\text{TSC} = f_{\text{education}}(X,\text{VSY}) + f_{\text{lives}}(D,l,\text{VLY},r) + f_{\text{GDP}}(Z_0-Z,r), \label{eq:swf}`

where the arguments are the number of school years lost ($X$) and the
value of a school year (VSY); the number of deaths per age group due to
COVID-19 ($D$), the remaining life expectancy per age group ($l$), the
value of a life year (VLY), and discount rate $r$; and the lost output
over the period due to reduced economic activity ($Z_0-Z$).

The cost includes both immediate costs, and mid- to long-term costs, all
attributable to choices made about closures over a short time period.
The losses are added up over all years for which the loss is expected to
persist, using the same discounting rate $r$ of 3%, which is implicit in
the VSY and applied to the other losses. As the counterfactual
no-pandemic scenario has a cost of 0, Equation represents the sum total
socio-economic costs of the pandemic over the period in question.  
The individual components $f_{\text{lives}}$, $f_{\text{GDP}}$, and
$f_{\text{education}}$ are described in full in Sections , , and ,
respectively.

### 4.0.1 Valuing lost lives

To value lives lost, we make use of the expected remaining life years
per age group . These are used to estimate the expected number of years
of life lost per death, and to estimate the value of a life year. We map
the remaining life expectancy $l_a$ for the GBD age groups $a$ to $l_g$
for the model age groups $g$ as a population-weighted average, taking
into account the size of each age group, $N_a$. For the expected number
of life years lost per death, we take into account also the probability
to die given infection, $P(D|I,a)$ (Table ):
`mathl_g^{\text{(death)}} = \frac{\sum_{a\in g}N_al_aP(D|I,a)}{\sum_{a\in g}N_aP(D|I,a)};`
`mathl_g^{\text{(life)}} = \frac{\sum_{a\in g}N_al_a}{\sum_{a\in g}N_a};`

Expected life years remaining with discounting taken into account can be
written $$\begin{align}
\hat{l}_g&=\sum_{y=1}^{l_g}\frac{1}{(1+r)^{y}}%\\
\end{align}$$ for discount rate $r>0$. The discounted number of years
lost given $D_g$ deaths due to COVID-19 for each age group is

``` {.mathy=\\sum_gd_g\\hat{l}_g^{\\text{(death)}}.```}
%We use this notation to define the cost as a simple product of the VLY and the discounted number of years for input into Equation \ref{eq:swf}.

The VLY used by policy makers should reflect the value that members of the society place on reductions of their own mortality. %We do not take into account any health costs apart from lives lost due to COVID-19. 
We rely on the intrinsic rather than instrumental interpretation of the valuation of life \cite{Cutler2020}, and we use existing estimates of the value of a statistical life (VSL) to estimate VLY. We interpret the VSL as a population-weighted average \citemref{Ananthapavan2021,Robinson2021}, where each age group has a VSL defined by the number of expected life years remaining, and where each discounted year has the same value: 
  ```math\text{VSL}=\frac{\sum_gN_g\hat{l}_g^{\text{(life)}}}{\sum_gN_g}\text{VLY}.```

Finally, we define the value of lives lost as 
```mathf_{\text{lives}}(D,l^{\text{(death)}},\text{VLY},r)=\text{VLY}*Y.```
In the applications, we report the total YLLs,  ```math\text{YLL}=\sum_gD_gl_g^{\text{(death)}}, ``` in their natural units of years, and, separately, the monetized life years lost, $f_{\text{lives}}(D,l^{\text{(death)}},\text{VLY},r)$, in \$.

We identify three estimates for the VSL in Indonesia. For the example in Section \ref{illustration}, we use the value of \$950,000 for a life of a working-age individual \cite{Wulandari}. %The author describes this value as ``high'', contrasting it with several estimates for comparable regions.
With 38 years remaining life expectancy for this group, we get an estimate for the VLY of \$42,000. For the scenarios, we use a range of values from \$20,000 to \$80,000 that encompass upper and lower limits for the VLY. The upper limit comes from an estimate of the VSL being worth 160 times GNI per capita (based on purchasing power parity \citemref{Robinson2021}), which we evaluate to be \$1,910,000, giving the VLY as \$80,000. For the lower limit, we use an estimate of \$592,000 for the VSL \citemref{Viscusi2017}. 


### Valuing lost economic output

%% define GVA

We measure the cost of economic closures in terms of lost gross value added (GVA): the GDP generated by an economic configuration is the maximum GVA per month (denoted $z_i$ for each sector $i$) multiplied by the respective sector openings, summed over the period. The maximum possible GDP (which is with no closures for $T$ months) is $Z_0=T\sum_{i}z_i$, and we use pre-pandemic output to define the maximum possible values.

Recall that $w_{i,\tau}$ denotes the proportion of the workforce contributing to economic production out of the total workforce $N_i$. The workforce can be additionally depleted due to sickness, hospitalisation and death, leaving a smaller fraction ($\hat{w}_{i,\tau}$) to contribute to production.\footnote{$\hat{w}_{i,\tau}=w_{i,\tau}-\sum_{v}\mathbb{E}_{d\in \tau}\left(I_{S,v,i}(d)+H_{v,i}(d)+D_{v,i}(d)\right)/N_i$ for vaccine status $v$, infectious and symptomatic $I_S$, hospitalised $H$, deceased $D$, with total population $N_i$ for sector $i$. $\mathbb{E}_{d\in \tau}$ represents that we compute the average number for all days $d$ in month $\tau$.} The output for the remaining workforce is $\hat{x}_{i,\tau}=\hat{w}_{i,\tau}$. Then the total GVA for the period is: 
  ```math\label{eq:economiccost}
Z = \sum_{i\neq \text{ed}}^{\mathcal{N}}\sum_{\tau=1}^Tz_i\hat{x}_{i,\tau} + Tz_{\text{ed}}
```

and the GDP loss compared to the maximum is $Z_0-Z$. All economic
sectors contribute GVA according to the level they are open for
production, except for the education sector which contributes its
maximum possible monthly GVA, $z_{\text{ed}}$ per month. %Note that the
economic loss due to absence {(symptomatic infection, hospitalisation,
and death)} is negligible compared to other costs (Figure ).

### 4.0.2 Valuing lost education

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

<div id="ref-Haw2020" class="csl-entry">

Haw, David, Giovanni Forchini, Patrick Doohan, Paula Christen, Matteo
Pianella, Rob Johnson, Sumali Bajaj, et al. 2022. “<span
class="nocase">Optimizing social and economic activity while containing
SARS-CoV-2 transmission using DAEDALUS</span>.” *Nature Computational
Science* 2: 223–33. <https://doi.org/10.25561/83928>.

</div>

<div id="ref-Walker2020" class="csl-entry">

Walker, Patrick G. T., Charles Whittaker, Oliver J. Watson, Marc
Baguelin, Peter Winskill, Arran Hamlet, Bimandra A. Djafaara, et al.
2020. “<span class="nocase">The impact of COVID-19 and strategies for
mitigation and suppression in low- and middle-income countries</span>.”
*Science* 369 (6502): 413–22. <https://doi.org/10.1126/science.abc0035>.

</div>

</div>
