
<!-- README-Inputs.md is generated from README-Inputs.Rmd. Please edit that file -->

### Description of Input Files

This folder contains the formatted data to run the
[CIBW\_ME\_Model-ms\_code.R](scripts/CIBW_ME_Model-ms_code.R) multievent
mark-recapture model as well as the starting latent matrix required by
JAGS. See manuscript for detailed description of data.

-   [ms\_SH\_data.csv](inputs/ms_SH_data.csv): Formatted (.csv) data to
    run the multievent mark-recapture model
-   [start\_mat-ms\_SH\_data.csv](inputs/start_mat-ms_SH_data.csv):
    Formatted (.csv) “starting” latent matrix corresponding to the
    ms\_SH\_data.csv file

The starting latent matrix is one possible iteration of the true states
described by the observed sighting history data. JAGS requires this
starting matrix in order to run the model. Other starting latent
matrices could be used as long as the true states in the matrix could
give rise to the observed data.
