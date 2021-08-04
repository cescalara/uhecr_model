## Tests and Checks

- modify_exposure: Modifies the exposure by some factor. Requires rerun of precomputation tables and simulation.
- random_seed: runs simulations for different random seeds, and also looks at the normalized cumulative distribution.
- random_ra: randomize the right ascension of the UHECR dataset, as the exposure is independent of right ascension (m = m(dec)). Should return different associations (if it were the same, then the simulation is disregarding the UHECR positions).
- random_sources: simulate the UHECRs using a known source catalogue, then fit using a randomly generated source catalogue of same number. The source associations should be different (if it were similar, this would mean that the simulation returns the same results disregarding the sources).