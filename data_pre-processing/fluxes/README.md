# Fluxes Data Pre-processing
To replicate the process used to calculate metabolic fluxes, the following scripts must be executed in order:
- processData: estimates intracellular and carotenoid metabolic fluxes.
- modelFVA: uses results from processData and a genome scale constraint-based model of _Saccharomyces cerevisiae_ to estimate credible bounds of fluxes in the metabolic pathway, in particular SQLS, which is equivalent to ERG9b in the kinetic models.
- estimateERG9Flux: estimates ERG9b based on results of processData and modelFVA.

The scripts also require access to the biomass_and_metabolites_data folder

## output
- BiomassConversionStatistics: conversion factors statistics (biomass/L)/OD600 for each of the strains.
- BiomassConversionsStatistics: conversion factors statistics (gDCW/L)OD600 for each of the chemostat cultivations.
- BiomassMeasurementsOutliers: outliers of biomass measurements for chemostats.
- BiomassMeasurementsStatistics: OD600 statisics of chemostat cultivations.
- DilutionRatesStatistics: dilution rate statistics of chemostat cultivations.
- MetablitesProductionStatistics: metabolite titer statistics if chemostat cultivations.
- FluxesBiomassStatistics: flux statistics of metabolites in mmol/gDCW/h.
- FluxesVolumeStatisticss: flux statistics of metabolites in mmol/L/h (intracellular volume).
- FVAMeasuredResults: Results of FVA for the fluxes that were measured.
- FVAEstimatedResults: Results of relevant fluxes that were not measured and belong to the mevalonate pathway, synthesis of carotenoids or the beginning of the sterol pathway.
- FVAReaction: stoichiometry of reactions in the pathway.
- ERG9BiomassStatistics: estimated flux of ERG9b in mmol/gDCW/h.
- ERG9VolumeStatistics: estimated flux of ERG9b in mmol/L/h (intracellular volume).

## cobra_model
This folder contains the original constraint-based model (iMM904_contextualized), previously published by [Lopez et al. (2020)](https://doi.org/10.3389/fbioe.2020.578793), and the code (transforModel) used to generate a version that deletes reactions not included in beta-car strains (iMM904_b_car).