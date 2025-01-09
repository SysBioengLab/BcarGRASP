# BcarGRASP
## Introduction
BcarGRASP is set of kinetic models of beta-carotene production in recombinant _Saccharomyces cerevisiae_ strains. The data, code, and inputs used to generate the models are also included. 
Models were created using the General Reaction Assembly and Sampling Platform (GRASP) computational implementation in Matlab R2021a. They also required the use of COBRA v3.0, NLOPT, parallel computing toolbox, and the optimization toolbox of Matlab.
Publication: (pending)
Requirements:
- GRASP (https://github.com/biosustain/GRASP/tree/main)
- COBRA Toolbox v3.0 (https://github.com/opencobra/cobratoolbox)
- NLOPT (https://github.com/stevengj/nlopt)
- Parallel Computing Toolbox (Matlab)
- Optimization Toolbox (Matlab)
- AriaMx Software v2.0 ([Agilent Technologies, Inc.](https://www.agilent.com/en/product/real-time-pcr-%28qpcr%29/real-time-pcr-%28qpcr%29-instruments/ariamx-software-download), only for processing qPCR raw data)

## Folders
The data, code, inputs and models are divided according to the methodology described in the publication.
Each major folder has their own README files with details about their contents.

### biomass_and_metabolites_data
Excel files regarding chemostat cultivations biomass and metabolite measurements. Used for the calculation of metabolic fluxes.
### rt-qpcr_data
Excel and Aria files regarding chemostat cultivations transcript measurements. Used for the estimation of relative mRNA transcripts.
### data_pre-processing
Code used to process biomass, metabolite and transcript data to generate input for the models.
### kinetic_models
Kinetic models of beta-carotene production. The folder is a equivalent to the GRASP Toolbox, with additional code, modifications, inputs and outputs to generate the models. This [shortcut](kinetic_models/ABC-GRASP/GRASP-main/io/output) leads to the kinetic models for the different type of structures and growth rates employed

## Nomenclature
The experimental conditions follow two interchangeable nomenclatures across the files:
- XY: where X stands for the beta-car strain (2, 3 or 4) and Y the growth rate (A for low, and B for high)
- XDZ: where X stands for the beta-car strain (2, 3 or 4) and Z for the dilution rated of the chemostat cultivation (010 for 0.1 1/h, and 025 for 0.25 1/h).
