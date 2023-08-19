# Computational Approach for Structure Generation of Anisotropic Particles (CASGAP)

## Introduction
`casgpap` package originally developed by Nitant Gupta, and Arthi Jayaraman is a method that generates structures composed of anisotropic particles which adhere to the user provided distributions of particle parameters. In the current implementation, all particles are spheroidal with semi-axis lengths $a$ and $c$, such that it is symmetric about the $c-$ axis. The spheroids are characterized by a volumetric radius $R=\sqrt{a^2 c}$ and an aspect ratio $\gamma=c/a$.  The orientation of the $c-$ axis is characterized by a quaternion $q={\theta_w, \mathbf{W}}$ details of which can be found in the paper [1].

The volumetric radius $R$ and aspect ratio $\gamma$ are assumed to follow log-normal distributions, with mean values $R_\mu$ and $\gamma_\mu$, and their standard deviations $R_\sigma$ and $\gamma_\sigma$, respectively. The orientation of the spheroids are distributed using the von Mises-Fisher distribution, with $\kappa$ as the orientational anisotropy parameter, details are provided in the paper [1].

The method is subdivided into two phases. Phase I builds a population of particles that follow the prescribed distributions of parameters. Phase II then assembles the particles sequentially to avoid overlaps.

The original code was written in MATLAB. With the help of Sri Vishnuvardhan Reddy Akepati at University of Delaware, a python version of the same code is also made available.

## Citation

__If you use this code, please cite the reference:__

Original Article on CASGAP:  

[1] Gupta, Nitant; Jayaraman, Arthi. Computational Approach for Structure Generation of Anisotropic Particles (CASGAP) with Targeted Distributions of Particle Design and Orientational Order. Submitted to Nanoscale 2023 (accepted). https://doi.org/10.1039/D3NR02425C
