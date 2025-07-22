# Fresnel Diffraction Simulation

This project models Fresnel (near field) diffraction patterns of light passing through various apertures using two numerical integration techniques: **Gaussian quadrature** and **Monte Carlo integration**.

Developed as part of a physics and numerical methods course, the code simulates diffraction from **rectangular** and **circular** apertures, and visualizes 1D and 2D intensity distributions on the screen. This project blends **scientific computing**, **numerical analysis**, and **data visualization** in Python.

## Project Functions

-  Simulates near-field diffraction using Fresnel's approximation
-  Models both rectangular and circular apertures
-  Compares Gaussian quadrature vs. Monte Carlo integration for handling complex integrals
-  Demonstrates trade-offs between numerical accuracy and computational efficiency
-  Generates 1D and 2D intensity plots with `matplotlib`

## Technical Skills

- Numerical integration with `scipy.integrate.dblquad`
- Random sampling with Monte Carlo methods
- Python scripting and scientific visualization
- Performance analysis and error tolerance tuning

## Report

A formal report detailing the theory, implementation, results, and conclusions is included here:  
 [`Fresnel_diffraction_from_an_aperture_report.pdf`](./Fresnel_diffraction_from_an_aperture_report.pdf)

## Plot Outputs

<p align="center">
  <img src="images/fresnel_1d.png" width="450"/>
  <img src="images/mc_diffraction.png" width="450"/>
</p>


## How to Run

Make sure you have Python 3 and the following packages:
- numpy
- matplotlib
- scipy
