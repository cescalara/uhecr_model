# Impact of using the UHECR arrival energies to constrain source associations

[![DOI](https://zenodo.org/badge/166797043.svg)](https://zenodo.org/badge/latestdoi/166797043)

Code used in [arXiv:1811.06464](https://arxiv.org/abs/1811.06464). The directories here are organised to follow the paper and 
each section is documented using jupyter notebooks. The main part of the project is a joint hierarchical model for the UHECR 
energies and arrival directions which in implemented in [Stan](https://mc-stan.org). The model includes UHECR production, 
propagation and detection effects.

If you use this code in your work please cite the [paper](https://arxiv.org/abs/1811.06464) and software using the above DOI.

## Dependencies

* [Stan](https://mc-stan.org) and [PyStan](https://pystan2.readthedocs.io/en/latest/)
* [CRPropa 3](https://github.com/CRPropa/CRPropa3) - for fits to simulations run using CRPropa 3
* [basemap](https://matplotlib.org/basemap/users/installing.html) - for sky plots
* [fancy](https://github.com/cescalara/fancy) - a python wrapper to make things (a bit) cleaner
* [stan_utility](https://github.com/grburgess/stan_utility) - utility package for pystan
* [Jupyter](https://jupyter.org) - for notebooks

## Update

This code is forked over from `cescalara/uhecr_model`. The code is updated so that it runs with the latest versions of:
- basemap
- h5py
- seaborn
- matplotlib
- pystan2 (this code does not work with pystan3)
- CRPropa3

All requirements for the specific package versions is given in `uhecr_project_env.yml` and can be easily replicated by creating a new conda environment
with this file. 

The utility packages `cescalara/fancy` and `grburgess/stan_utility` have also been updated and forked over to `uhecr_project/fancy`, `uhecr_project/stan_utility`. 
Please use these packages with this code instead.

## Questions?

If you have any questions feel free to open an issue or get in touch (capel.francesca@gmail.com).
