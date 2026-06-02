
# SFPLSR

The **SFPLSR** package provides a comprehensive framework for spatial
functional partial least squares regression. It is designed for modeling
scalar-on-function regression problems in the presence of spatial
dependence.

## Abstract

This study proposes a spatial functional partial least squares framework
for estimating a spatial autoregressive scalar-on-function regression
model, where the underlying structure is inherently
infinite-dimensional. The method embeds spatial dependence into the
dimension reduction process by constructing latent components that
maximize the covariance between the scalar response and the functional
predictor under spatial autocorrelation, yielding a finite-dimensional
representation of the model. The proposed approach provides a unified
framework that combines the dimension-reduction advantages of classical
functional partial least squares with explicit modeling of spatial
dependence during component extraction. Its finite-sample performance is
evaluated through Monte Carlo simulations and an empirical application
using weather data from the United Kingdom. The results demonstrate that
the proposed method consistently outperforms existing alternatives in
terms of estimation and prediction accuracy, particularly under strong
spatial dependence.

## Installation

You can install the development version of **SFPLSR** directly from
GitHub:

    install.packages("devtools"); pak::pak("MugeMutis/SFPLSR")

## Main Functions

The package provides an end-to-end workflow for spatial functional data
analysis, from data generation to model estimation and prediction.

#### spatial_data_generation()

Generates synthetic spatial functional data for simulation studies and
methodological validation. The function constructs infinite-dimensional
functional predictors together with a scalar response variable while
explicitly incorporating spatial autoregressive dependence through a
spatial weights matrix.

#### sfplsr()

Implements the proposed Spatial Functional Partial Least Squares
Regression (SFPLSR) framework. The method extracts latent components
that maximize the covariance between the scalar response and the
functional predictors while accounting for spatial autocorrelation.

#### predict_sfplsr()

Computes out-of-sample predictions from a fitted SFPLSR model.

## Dataset

#### UK_weather_data()

The package includes a real-world environmental and meteorological
dataset compiled from the NASA Prediction Of Worldwide Energy Resources
(POWER) project.

The dataset contains observations from 231 meteorological stations
distributed across the United Kingdom and includes:

- Spatial Information: Latitude and longitude coordinates for each
  station. Suitable for constructing spatial weights matrices and
  modeling spatial dependence.

- Functional Covariates: Daily temperature trajectories observed during
  2024 and 2025. Treated as infinite-dimensional functional predictors.

- Scalar Response: Surface shortwave downward flux (mean solar
  radiation) measurements for 2024 and 2025. Used as the scalar response
  variable in the spatial autoregressive scalar-on-function regression
  model.

## References

M. Mutis, U. Beyaztas, H. L. Shang (2026). *Spatial functional partial
least squares regression*.

## Contact

Muge Mutis - <muge.mutis@yildiz.edu.tr>
