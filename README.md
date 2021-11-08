
<!-- README.md is generated from README.Rmd. Please edit that file -->

### Estimating a reproductive rate when offspring ages are uncertain: A novel multievent mark-recapture model applied to an endangered beluga whale population

#### Himes Boor, G.K., T.L. McGuire, A.J. Warlick, R.L. Taylor, S.J. Converse, J.R. McClung, A.D. Stephens

For questions regarding the code contact [Gina Himes
Boor](mailto:gkhimesboor@montana.edu)

#### Citation:

Himes Boor, G.K., T.L. McGuire, A.J. Warlick, R.L. Taylor, S.J.
Converse, J.R. McClung, A.D. Stephens. Estimating a reproductive rate
when offspring ages are uncertain: A novel multievent mark-recapture
model applied to an endangered beluga whale population.

------------------------------------------------------------------------

#### Abstract

1.  Understanding survival and reproductive rates of a population is
    critical to determining its long-term dynamics and viability. For
    some hard-to-study species, obtaining the mark-recapture data
    necessary to estimate these vital rates can be difficult or
    impossible. This is particularly true for reproductive rates in
    populations where offspring are not always detected, remain with
    parents for multiple years, or cannot be aged with certainty.
    Existing reproductive-rate models cannot currently accommodate these
    situations.
2.  To address this, we developed a novel multievent mark-recapture
    model that uses all available adult and adult-offspring sightings,
    including sightings with older offspring of uncertain age. We
    applied this model to 13 years of photo-ID data from a critically
    endangered population of beluga whales (Delphinapterus leucas) to
    estimate reproductive rates, as well as survival rates for breeding
    females, non-breeders (male and female), and calves. We also
    evaluated model performance using simulated data.
3.  Our results indicated that breeding female belugas that did not give
    birth the previous year have a reproductive rate of 0.279 (95%
    credible interval: 0.226 – 0.34) and exhibit a higher survival rate
    (0.962; 0.945 – 0.975) than non-breeders (0.931; 0.917 – 0.944).
    Non-breeders become breeders at a rate of about 0.072 (0.059 –
    0.085) per year. Survival of calves ≤1 year was 0.926 (0.85 –
    0.964), while apparent survival (i.e., surviving and staying with
    mother) of older calves (≥2 years old) was 0.492 (0.398 – 0.595).
    Our simulation analysis indicated our model is asymptotically
    unbiased. Even at moderate sample sizes and detection rates, the
    ecological parameters of interest were estimated with good precision
    and minimal bias.
4.  This work yields an important new development in multievent
    mark-recapture modeling, allowing estimation of reproductive rates
    and juvenile survival for populations with extended adult-offspring
    associations and uncertain offspring ages. Additionally, our model
    provides new robust demographic rate estimates for an endangered
    population that indicate low reproduction and non-breeder survival
    are likely both contributing to the population’s continued decline.
    Breeder survival may also be a factor, but additional work is
    necessary to understand the relative contributions of each
    demographic rate and to identify exogenous stressors.

------------------------------------------------------------------------

#### Code:

[CIBW\_ME\_Model-ms\_code.R](https://github.com/gkhimesboor/HimesBoor_etal_2021_MEE/blob/9976e998f1fb87a36f7e40a1775889215e14dc55/scripts/CIBW_ME_Model-ms_code.R):
R code for the multievent mark-recapture model developed to estimate
survival and reproductive rates from the Cook Inlet beluga whale
photo-ID data, and that can be adapted for other data with similar
characteristics (i.e., mark-recapture study involving marked adults that
are sometimes observed with their unmarked offspring of unknown age)

#### Data:

[ms\_SH\_data.csv](https://github.com/gkhimesboor/HimesBoor_etal_2021_MEE/blob/9976e998f1fb87a36f7e40a1775889215e14dc55/inputs/ms_SH_data.csv):
formatted Cook Inlet beluga whale photo-ID mark-recapture data for
running the multievent model described above. The data were collected by
Dr. Tamara McGuire and colleagues at [The Cook Inlet Beluga Whale
Photo-ID Project](https://www.cookinletbelugas.com/) and should not be
used outside of this analysis without express permission from
[Dr. McGuire](mailto:tamaracookinletbeluga@gmail.com).

#### Additional Required Files:

[start\_mat-ms\_SH\_data.csv](https://github.com/gkhimesboor/HimesBoor_etal_2021_MEE/blob/9976e998f1fb87a36f7e40a1775889215e14dc55/inputs/start_mat-ms_SH_data.csv):
formatted starting latent matrix required by JAGS to run the multievent
model (see
[Description.txt](https://github.com/gkhimesboor/HimesBoor_etal_2021_MEE/blob/9976e998f1fb87a36f7e40a1775889215e14dc55/inputs/DESCRIPTION.txt)
for more information about this file)

------------------------------------------------------------------------

#### Funding

Funding for development of this model came from the [North Pacific
Research Board](https://www.nprb.org/) (project \# 1718). Cook Inlet
beluga photo-ID data were collected using a wide range of funding
sources.
