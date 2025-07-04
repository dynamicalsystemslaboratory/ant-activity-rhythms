**README: Journal Article Materials**

**Title:** *The role of colony size on activity rhythms of ants*\
**Authors:** Michael Napoli, Rifat Sipahi, Maurizio Porfiri\
**Corresponding Author:** Maurizio Porfiri\
**Date:** July 3, 2025

---

### Directory Structure

* `data/` Contains data necessary to generate figures/perform analysis.
    * `doering2024/` Experimental data from [28] on the activity of *T. rudis* ants.
    * `results/` Data created during the analysis and saved for figure recreation.
        * `exponent-calibration/` Exponents computed during model calibration.
        * `parametric-gamm/` Oscillation period and burst height computed over the $N$-$\gamma$ space.

* `figures/` Scripts used to perform analysis and generate the figures included in the manuscript.
    * `figure1.ipynb` Proportion of active ants as a function of time in a colony of *T. rudis*.
    * `figure2.ipynb` Analysis of the role of colony size on rhythmic activity, based on experimental results by Doering et al. [28].
    * `figure4.ipynb` Example of model predictions for different colony sizes, keeping all other parameters fixed.
    * `figure5.ipynb` Critical delay, $\tau_c$, and corresponding frequency, $\omega_c$, for colonies with varying size, $N$.
    * `figure6.ipynb` Scaling of the oscillation period and burst height in the model and for varying $N$ and $\gamma$.
    * `figure7.ipynb` Scaling of oscillation period and burst height from the compartmental model, calibrated on experimental data.

* `FOUR/` Custom library used to compute the discrete Fourier transform.

---

### Notes

* All scripts are written in Python. Necessary libraries and versions include:
    * Numpy `1.26.0`
    * Matplotlib `3.9.2`
    * Statsmodels `0.14.2`
    * Joblib `1.4.2`
    * Scipy `1.13.1`
