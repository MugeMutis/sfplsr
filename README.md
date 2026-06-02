
SFPLSR

The goal of SFPLSR is to provide a comprehensive framework for spatial
functional partial least squares regression.

Abstract

This study proposes a spatial functional partial least squares framework
for estimating a spatial autoregressive scalar-on-function regression
model, where the underlying structure is inherently
infinite-dimensional. The method embeds spatial dependence into the
dimension reduction process by constructing latent components that
maximize the covariance between the scalar response and the functional
predictor under spatial autocorrelation, yielding a finite-dimensional
representation of the model. This provides a unified approach that
combines the dimension-reduction advantages of classical functional
partial least squares with explicit modeling of spatial dependence at
the component extraction stage. The finite-sample performance of the
proposed approach is evaluated through Monte Carlo simulations and an
empirical application using weather data from the United Kingdom. The
results show that the proposed approach consistently outperforms
existing methods in terms of estimation and prediction accuracy,
particularly under strong spatial dependence.

Installation

You can install the development version of SFPLSR from GitHub with:

R \# install.packages(“devtools”)
devtools::install_github(“MugeMutis/SFPLSR”)

Main Functions & Workflow

The core functionality of the SFPLSR package is centered around an
end-to-end workflow for spatial functional data analysis, from data
generation to model estimation and prediction:

spatial_data_generation: Generates synthetic spatial functional data for
simulation studies and framework validation. It constructs
infinite-dimensional functional predictors paired with a scalar response
variable, explicitly incorporating spatial autoregressive dependence and
spatial autocorrelation structures based on a spatial weights matrix.

sfplsr: Implements the core spatial functional partial least squares
regression framework. This function performs spatial dimension reduction
by extracting latent components that maximize the covariance between the
scalar response and the functional covariates while accounting for
spatial autocorrelation.

predict_sfplsr: Computes out-of-sample predictions for the estimated
spatial autoregressive functional partial least squares regression.

Datasets

UK_weather_data: A comprehensive environmental and meteorological
dataset compiled from the NASA Prediction Of Worldwide Energy Resources
(POWER) project. It contains spatial and functional observations for 231
meteorological stations homogeneously distributed across the United
Kingdom. For each station, the dataset provides:

Spatial Attributes: Precise geographical coordinates (Latitude and
Longitude) to construct spatial weights matrices and model spatial
dependency.

Functional Covariates: Daily temperature trajectories (curves) captured
over the years 2024 and 2025, representing the infinite-dimensional
functional predictors.

Scalar Response: Surface shortwave downward flux (mean solar radiation)
values for the years 2024 and 2025, serving as the scalar response
variable in the spatial autoregressive scalar-on-function regression
framework.

Contact Muge Mutis - <muge.mutis@yildiz.edu.tr>
