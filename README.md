
<!-- README.md is generated from README.Rmd. Please edit that file -->

### Estimating reproductive and juvenile survival rates when offspring ages are uncertain: a novel multievent mark-resight model with beluga whale case study

#### Himes Boor, G.K., T.L. McGuire, A.J. Warlick, R.L. Taylor, S.J. Converse, J.R. McClung, A.D. Stephens

For questions regarding the code, contact [Gina Himes
Boor](mailto:gkhimesboor@montana.edu)

#### Citation:

Himes Boor, G.K., T.L. McGuire, A.J. Warlick, R.L. Taylor, S.J.
Converse, J.R. McClung, A.D. Stephens. 2023. Estimating reproductive and
juvenile survival rates when offspring ages are uncertain: a novel
multievent mark-resight model with beluga whale case study. Methods in
Ecology and Evolution 14(2):631-642.

------------------------------------------------------------------------

#### Abstract

1.  Understanding the survival and reproductive rates of a population is
    critical to determining its long-term dynamics and viability.
    Mark-resight models are often used to estimate these demographic
    rates, but estimation of survival and reproductive rates is
    challenging, especially for wide-ranging, patchily distributed, or
    cryptic species. In particular, existing mark-resight models cannot
    accommodate data from populations in which offspring remain with
    parents for multiple years, are not always detected, and cannot be
    aged with certainty.
2.  Here we describe a Bayesian multievent mark-resight modelling
    framework that uses all available adult and adult-offspring
    sightings (including sightings with older offspring of uncertain
    age) to estimate reproductive rates and survival rates of adults and
    juveniles. We extend existing multievent mark-resight models that
    typically only incorporate adult breeding state uncertainty by
    additionally accounting for age uncertainty in unmarked offspring
    and uncertainty in the duration of the mother-offspring association.
    We describe our model in general terms and with a simple
    illustrative example, then apply it in a more complex empirical
    setting using thirteen years of photo-ID data from a critically
    endangered population of beluga whales (Delphinapterus leucas). We
    evaluated model performance using simulated data under a range of
    sample sizes, and adult and offspring detection rates.
3.  Applying our model to the beluga data yielded precise estimates for
    all demographic rates of interest despite substantial uncertainty in
    calf ages, including non-breeder survival and reproductive rates
    lower than that estimated for other beluga populations. Simulations
    suggested our model yields asymptotically unbiased parameter
    estimates with good precision and low bias even with moderate sample
    sizes and detection rates.
4.  This work represents an important new development in multievent
    mark-resight modeling, allowing estimation of reproductive and
    juvenile survival rates for populations with extended adult –
    offspring associations and uncertain offspring ages (e.g., some
    marine mammals, elephants, bears, great apes, bats, and birds). Our
    model facilitated estimation of robust demographic rates for an
    endangered beluga population that were previously inestimable (e.g.,
    non-breeder and juvenile survival, reproductive rate) and that will
    yield new insights into this population’s continued decline.

------------------------------------------------------------------------

#### CODE:

**Model:**  
[CIBW_ME_Model-ms_code.R](scripts/CIBW_ME_Model-ms_code.R): R code for
the multievent mark-recapture model developed to estimate survival and
reproductive rates from sighting history data of marked adults that are
sometimes observed with their unmarked offspring of unknown age. The
model was originally developed for photo-ID data from Cook Inlet beluga
whales, but it can be adapted for other populations with similar
characteristics (i.e., extended parental care of offspring that can not
always be aged with certainty).

**Simulations:**  
[CIBW_ME_sim.R](scripts/CIBW_ME_sim.R): R code for the accompanying
simulation analysis for the multievent mark-recapture model examining
model performance across varying detection rates and sample population
sizes.

#### DATA:

[ms_SH_data.csv](inputs/ms_SH_data.csv): formatted Cook Inlet beluga
whale photo-ID mark-recapture data for running the multievent model
described above. The data were collected by Dr. Tamara McGuire and
colleagues at [The Cook Inlet Beluga Whale Photo-ID
Project](https://www.cookinletbelugas.com/) and should not be used
outside of this analysis without express permission from
[Dr. McGuire](mailto:tamaracookinletbeluga@gmail.com).

#### Additional Required Files:

[start_mat-ms_SH_data.csv](inputs/start_mat-ms_SH_data.csv): formatted
starting latent matrix required by JAGS to run the multievent model (see
[Inputs README file](inputs/README.md) for more information about this
file)

------------------------------------------------------------------------

#### Funding

Funding for development of this model came from the [North Pacific
Research Board](https://www.nprb.org/) (project \# 1718). Cook Inlet
beluga photo-ID data were collected using a variety of funding sources.
