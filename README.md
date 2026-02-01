[READ ME.txt](https://github.com/user-attachments/files/24991893/READ.ME.txt)
# TFM – PHITS and Python codes for ICF transport simulations

This repository contains the input files and analysis scripts developed for the
Master's Thesis focused on the modelling and Monte Carlo simulation of particle
and radiation transport in indirect-drive inertial confinement fusion (ICF).

The simulations were performed using the PHITS Monte Carlo code, and the results
were post-processed using custom Python scripts.

---

## Repository structure

### PHITS/
PHITS input files defining the physical and geometrical model at stagnation,
including sources and tallies for different fusion products:

- Neutron transport simulations (14.1 MeV source).
- Alpha particle transport simulations (3.5 MeV source).


---

### python/
Python scripts developed for post-processing and analysis of PHITS outputs:

- `PHITS_Source_generator.py`  
  Generates the particle source distribution used in PHITS simulations.

- `Neutron_spectrum.ipynb`  
  Reads and analyses neutron energy spectra from PHITS tallies.

- `Photon_spectrum.ipynb`  
  Reads and analyses photon energy spectra from PHITS tallies.

- `t_product_neutron.py`  
  Analyses neutron production from tallies

- `t_product_photons.py`  
  Analyses photon production from tallies

---

## Software and versions

- PHITS: version 3.34
- Python: version 3.13

