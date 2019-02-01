# Impact of using the UHECR arrival energies to constrain source associations

Code used in [arXiv:1811.06464](https://arxiv.org/abs/1811.06464). The directories here are organised to follow the paper and 
each section is documented using jupyter notebooks. The main part of the project is a joint hierarchical model for the UHECR 
energies and arrival directions which in implemented in [Stan](https://mc-stan.org). The model includes UHECR production, 
propagation and detection effects.

If you use this code in your work please cite the above paper.

## Dependencies

* [Stan](https://mc-stan.org) and [PyStan](https://pystan.readthedocs.io/en/latest/)
* [CRPropa 3](https://github.com/CRPropa/CRPropa3) - for fits to simulations run using CRPropa 3
* [basemap](https://matplotlib.org/basemap/users/installing.html) - for sky plots
* [fancy](https://github.com/cescalara/fancy) - a python wrapper to make things (a bit) cleaner
* [Jupyter](https://jupyter.org) - for notebooks

## Questions?

If you have any questions feel free to open an issue or get in touch (capel.francesca@gmail.com).
