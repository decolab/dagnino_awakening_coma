# Re-awakening the Brain for Disorders of Consciousness (DoC)

**Matlab implementation of the Leading Eigenvector Dynamics Analysis (LEiDA) framework and Whole-brain Hopf model**

This repository contains the source code used in the study:

> Dagnino, P. C., Escrichs, A., López-González, A., Gosseries, O., Annen, J., Sanz Perl, Y., Kringelbach, M. L., Laureys, S., & Deco, G. (2024). *Re-awakening the brain: Forcing transitions in disorders of consciousness by external in silico perturbation*. PLoS Computational Biology, 20(5), e1011350. [https://doi.org/10.1371/journal.pcbi.1011350](https://doi.org/10.1371/journal.pcbi.1011350)

This study explores the probabilistic metastable substates (PMS) of different diagnostic caterogories of DoC (model-free) and how in silico perturbation of whole-brain Hopf models can drive transitions in patients with DoC (model-based).

---

## Purpose & Scope

- Model transitions between low- and high-consciousness brain states via in silico perturbations.
- Simulate whole-brain activity using effective connectivity derived from empirical fMRI and structural connectomes.
- Assess perturbation outcomes across diagnostic categories (e.g., UWS, MCS).
- Investigate the regional sensitivity to perturbation that is driving the transition.

---

## Repository Structure

- `hopf_simulations/`: MATLAB scripts for simulating perturbation protocols.
- `leida/`: MATLAB scripts for calculating the PMS.

> Note: Original patient data are not included due to privacy. Use placeholder or simulated data to test scripts.

---

## ⚙️ Requirements

- MATLAB R2018a or later
- resting-state fMRI and dMRI connectome matrices (preprocessed and in time series format)


---

## Usage Instructions

1. Clone the repository and add paths to MATLAB.

2. Load the empirical connectome and fMRI baseline state.

3. Run the LEiDA scripts.

4. Run the whole-brain model scripts

5. Analyze state transitions and the regions most sensitive to perturbation.

---

## Data Sources

- This work uses clinical neuroimaging data from patients with DoC.
- All patient data were anonymized and handled in accordance with institutional ethics guidelines.

---

## Citation

If you use this code, please cite:

Dagnino, P. C., Escrichs, A., López-González, A., Gosseries, O., Annen, J., Sanz Perl, Y., Kringelbach, M. L., Laureys, S., & Deco, G. (2024). *Re-awakening the brain: Forcing transitions in disorders of consciousness by external in silico perturbation*. PLoS Computational Biology, 20(5): e1011350. [https://doi.org/10.1371/journal.pcbi.1011350](https://doi.org/10.1371/journal.pcbi.1011350)

---

## License

This project is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License.

---

## Contact

For support or scientific questions, please contact the first author:

Paula Clara Dagnino\
[paulinaclara.dagnino@upf.edu](mailto\:paulinaclara.dagnino@upf.edu)\
CNS Department, University Pompeu Fabra

