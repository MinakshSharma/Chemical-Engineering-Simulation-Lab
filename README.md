# Diffusion-Reaction Analysis using MATLAB

## Overview

This repository contains the MATLAB implementation developed as part of the **CH230: Chemical Engineering Simulation Lab** course at the **Indian Institute of Technology Ropar**.

The project focuses on the numerical solution of a one-dimensional reaction–diffusion equation with spatially varying diffusivity. The numerical solution is compared with the corresponding analytical solution, and cubic spline interpolation is used for smooth visualization of the concentration profile.

---

## Objectives

- Solve the governing reaction–diffusion equation numerically.
- Compare numerical and analytical solutions.
- Visualize concentration profiles using spline interpolation.
- Understand the application of numerical methods in chemical engineering.

---

## Numerical Methods Used

- Finite Difference Method (FDM)
- Gaussian Elimination with Partial Pivoting
- Cubic Spline Interpolation
- Analytical Solution using Bessel Function Approximation

---

## Software

- MATLAB

---

## Repository Contents

```
.
├── project.matlab.m                 # MATLAB source code
├── Diffusion-Reaction Analysis (Revised).pdf   # Project report
└── README.md
```

---

## Results

The MATLAB implementation:

- Computes the numerical concentration profile using the Finite Difference Method.
- Implements Gaussian Elimination with Partial Pivoting without relying on MATLAB's built-in linear system solvers.
- Uses Cubic Spline Interpolation for smooth concentration profiles.
- Compares the numerical solution with the analytical solution obtained using Bessel function approximations.

---

## Team Members

- Dhruv Kaushal (2024CHB1066)
- Cheemakurthi Purandhar Sai (2024CHB1062)
- Ganedi Durga Tejaswi (2024CHB1067)
- **Minaksh Sharma (2024CHB1072)**
- Sunitha Nayak (2024CHB1086)

---

## Disclaimer

This repository contains the coursework completed as part of the CH230 Chemical Engineering Simulation Lab at IIT Ropar. It is intended for educational purposes only.
