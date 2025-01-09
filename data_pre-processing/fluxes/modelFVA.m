%% This scripts estimates the non measured fluxes for GRASP (SQLS)
% SQLS is used as an approximation for flux ERG9b in the GRASP models
% Requires the Cobra Toolbox (Heirendt et al., 2019)

% It requires as inputs:
% - IMM904_b_car.mat: Constraint-based model included in the publication of
% Lopez et al. (2020) previously treated to delete beta-ionone related
% reactions
% - DilutionRatesStatistics: estimated dilution rates of the chemostat
% cultivations
% - FluxesBiomassStatistics: estimated fluxes per biomass of the
% metabolites included in the model
% Other input variables are defined to select the fluxes to be modified

% It generates in the 'output' folder the outputs:
% - FVAMeasuredResults: Results of FVA for the fluxes that were measured
% - FVAEstimatedResults: Results of relevant fluxes that were not measured
% and belong to the mevalonate pathway, synthesis of carotenoids or the
% beginning of the sterol pathway
% - FVAReaction: stoichiometry of reactions in the pathway

%% Clean variables
clc,clearvars,close all

%% Load data and initialize variables
% Initialize CobraToolbox
% initCobraToolbox

% Load GEM of the beta-ionone producing yeast
load('cobra_model\iMM904_b_car.mat')
changeCobraSolver("matlab");

% Show new reactions added (10 in total)
model.rxnNames(end-9:end)

% Define Variables
strainConditions = {'2D01';'2D025';'3D01';'3D025';'4D01';'4D025'};                                  % names of the chemostat conditions
measuredMetabolites = {'Glucose';'Ethanol';'Acetate';'Glycerol';'Lycopene';'BetaCarotene'};         % metabolite names
fluxesMetabolites = {'EX_glc(e)';'EX_etoh(e)';'EX_ac(e)';'EX_glyc(e)';'EX_lyc(c)';'EX_b_car(c)'};   % names of the fluxes associated with the metabolites
alpha = 0.05;                                                                                       % significance to estimate confidence intervals

% Reaction List
listRxns = {'ACACT1';'ACACT1m';'HMGCOAS';'HMGCOASm';'HMGCOAR';'MEVK1';'MEVK2';'MEVK3';'MEVK4';
    'PMEVK';'DPMVD';'IPDDI';'DMATT';'GRTT';'FRTT';'SQLS';'PHYS';
    'PHYD_z_car(c)_forming';'PHYD_lyc(c)_forming';'LICC_g_car(c)_forming';'LICC_b_car(c)_forming';
    'EX_phy(c)';'EX_z_car(c)';'EX_lyc(c)';'EX_g_car(c)';'EX_b_car(c)';'EX_ergst(e)'};               % list of relevant reactions for the model

% Load tables with relevant values or their names
DilutionRatesTable = readtable('output\DilutionRatesStatistics.xlsx');
FluxesBiomassFilename = 'output\FluxesBiomassStatistics.xlsx';

%% Modify model parameters

% Fix rxn bounds for maximum yield calculation
model = changeRxnBounds(model,'ATPM',0.7,'l');            	% based on last yeastGEM
model = changeRxnBounds(model,'EX_co2(e)',0,'l');
model = changeRxnBounds(model,'EX_o2(e)',-1000,'l');
model = changeRxnBounds(model,'EX_b_car(c)',0,'l');
model = changeRxnBounds(model,'EX_lyc(c)',0,'l');
model = changeRxnBounds(model,'biomass_SC5_notrace',0,'l');
model = changeRxnBounds(model,'EX_glc(e)',-10,'b');
model = changeRxnBounds(model,'EX_co2(e)',1000,'u');
model = changeRxnBounds(model,'ATPM',1000,'u');
model = changeRxnBounds(model,'biomass_SC5_notrace',1000,'u');
model = changeRxnBounds(model,'EX_phy(c)',0,'b');
model = changeRxnBounds(model,'EX_z_car(c)',0,'b');
model = changeRxnBounds(model,'EX_lyc(c)',1000,'u');
model = changeRxnBounds(model,'EX_g_car(c)',0,'b');
model = changeRxnBounds(model,'EX_b_car(c)',1000,'u');

model = changeRxnBounds(model,'EX_ergst(e)',0,'b');
model = changeRxnBounds(model,'MEVK2',0,'b');
model = changeRxnBounds(model,'MEVK3',0,'b');
model = changeRxnBounds(model,'MEVK4',0,'b');
model = changeRxnBounds(model,'HMGCOASm',0,'b');
model = changeRxnBounds(model,'ACACT1m',0,'b');

model = changeRxnBounds(model,'EX_lanost(e)',0,'b');
model = changeRxnBounds(model,'EX_lanostest_SC(e)',0,'b');
model = changeRxnBounds(model,'EX_zymst(e)',0,'b');
model = changeRxnBounds(model,'EX_zymstest_SC(e)',0,'b');
model = changeRxnBounds(model,'EX_fecost(e)',0,'b');
model = changeRxnBounds(model,'EX_fecostest_SC(e)',0,'b');
model = changeRxnBounds(model,'EX_epist(e)',0,'b');
model = changeRxnBounds(model,'EX_epistest_SC(e)',0,'b');
model = changeRxnBounds(model,'EX_ergst(e)',0,'b');
model = changeRxnBounds(model,'EX_ergstest_SC(e)',0,'b');

%% Extract Dilution Rate information and calculate CIs

dilutionRatesMean = DilutionRatesTable{:,2};
dilutionRatesSE = DilutionRatesTable{:,3};
dilutionRatesDF = DilutionRatesTable{:,5};

dilutionRatesT = tinv(1-alpha/2,dilutionRatesDF);
dilutionRatesLCI = dilutionRatesMean-dilutionRatesSE.*dilutionRatesT;
dilutionRatesUCI = dilutionRatesMean+dilutionRatesSE.*dilutionRatesT;

%% Extract CIs of measured fluxes

measuredFluxesMean = zeros(length(strainConditions),length(measuredMetabolites));
measuredFluxesLCI = zeros(length(strainConditions),length(measuredMetabolites));
measuredFluxesUCI = zeros(length(strainConditions),length(measuredMetabolites));

for i=1:length(measuredMetabolites)
    FluxTable = readtable('output\FluxesBiomassStatistics.xlsx','Sheet',measuredMetabolites{i,1});
    measuredFluxesMean(:,i) = FluxTable{1:6,2};
    measuredFluxesLCI(:,i) = FluxTable{1:6,4};
    measuredFluxesUCI(:,i) = FluxTable{1:6,5};
end

%% Modify Model and Perform FVA for CI Estimation of ERG9 (SQLS)

% Define variables for FVA with the limits of the CIs
fvaMeasuredLCI = zeros(length(measuredMetabolites),length(strainConditions));
fvaMeasuredUCI = zeros(length(measuredMetabolites),length(strainConditions));

fvaEstimatedLCI = zeros(length(listRxns),length(strainConditions));
fvaEstimatedUCI = zeros(length(listRxns),length(strainConditions));

% Perfrom FVA with the limits of the CIs
for i=1:length(strainConditions)
    model = changeRxnBounds(model,fluxesMetabolites,measuredFluxesLCI(i,:)','l');
    model = changeRxnBounds(model,fluxesMetabolites,measuredFluxesUCI(i,:)','u');
    model = changeObjective(model,'EX_ergst(e)');
    model = changeRxnBounds(model,'biomass_SC5_notrace',dilutionRatesLCI(i,:),'l');
    model = changeRxnBounds(model,'biomass_SC5_notrace',dilutionRatesUCI(i,:),'u');

    [fvaLCI,fvaUCI] = fluxVariability(model,1);

    for j=1:length(measuredMetabolites)
        fvaMeasuredLCI(j,i) = fvaLCI(strcmp(model.rxns,fluxesMetabolites{j,1}));
        fvaMeasuredUCI(j,i) = fvaUCI(strcmp(model.rxns,fluxesMetabolites{j,1}));
    end
    
    for j=1:length(listRxns)
        fvaEstimatedLCI(j,i) = fvaLCI(strcmp(model.rxns,listRxns{j,1}));
        fvaEstimatedUCI(j,i) = fvaUCI(strcmp(model.rxns,listRxns{j,1}));
    end
end

% Define variables for FVA with the mean values of the fluxes
fvaMeasuredMeanL = zeros(length(measuredMetabolites),length(strainConditions));
fvaMeasuredMeanU = zeros(length(measuredMetabolites),length(strainConditions));

fvaEstimatedMeanL = zeros(length(listRxns),length(strainConditions));
fvaEstimatedMeanU = zeros(length(listRxns),length(strainConditions));

% Perfrom FVA with the mean values to check feasibility
for i=1:length(strainConditions)
    model = changeRxnBounds(model,fluxesMetabolites,measuredFluxesMean(i,:)','b');
    model = changeObjective(model,'EX_ergst(e)');
    model = changeRxnBounds(model,'biomass_SC5_notrace',dilutionRatesMean(i,:),'b');

    [fvaMeanL,fvaMeanU] = fluxVariability(model,1);

    for j=1:length(measuredMetabolites)
        fvaMeasuredMeanL(j,i) = fvaMeanL(strcmp(model.rxns,fluxesMetabolites{j,1}));
        fvaMeasuredMeanU(j,i) = fvaMeanU(strcmp(model.rxns,fluxesMetabolites{j,1}));
    end
    
    for j=1:length(listRxns)
        fvaEstimatedMeanL(j,i) = fvaMeanL(strcmp(model.rxns,listRxns{j,1}));
        fvaEstimatedMeanU(j,i) = fvaMeanU(strcmp(model.rxns,listRxns{j,1}));
    end
end

%% Prepare and Export Tables with the fluxes

% Table for Measured Fluxes
for i = 1:length(measuredMetabolites)
    fvaLCI = fvaMeasuredLCI(i,:)';
    fvaUCI = fvaMeasuredUCI(i,:)';
    fvaMCI = (fvaLCI+fvaUCI)./2;

    fvaMeanL = fvaMeasuredMeanL(i,:)';
    fvaMeanU = fvaMeasuredMeanU(i,:)';
    fvaMeanM = (fvaMeanL + fvaMeanU)./2;

    fvaTable = table(strainConditions,fvaLCI,fvaMCI,fvaUCI,fvaMeanL,fvaMeanM,fvaMeanU);
    fvaTable = renamevars(fvaTable,["strainConditions","fvaLCI","fvaMCI","fvaUCI","fvaMeanL","fvaMeanM","fvaMeanU"],["StrainCondition","Lower Median 95% CI","Median 95% CI","Upper 95% CI","Mean Lower CI 95%","Mean Median CI 95%","Mean Upper CI 95%"]);
    writetable(fvaTable,'output\FVAMeasuredResults.xlsx','Sheet',measuredMetabolites{i,1})
end

% Table for Estimated Fluxes
for i = 1:length(listRxns)
    fvaLCI = fvaEstimatedLCI(i,:)';
    fvaUCI = fvaEstimatedUCI(i,:)';
    fvaMCI = (fvaLCI+fvaUCI)./2;

    fvaMeanL = fvaEstimatedMeanL(i,:)';
    fvaMeanU = fvaEstimatedMeanU(i,:)';
    fvaMeanM = (fvaMeanL + fvaMeanU)./2;

    fvaTable = table(strainConditions,fvaLCI,fvaMCI,fvaUCI,fvaMeanL,fvaMeanM,fvaMeanU);
    fvaTable = renamevars(fvaTable,["strainConditions","fvaLCI","fvaMCI","fvaUCI","fvaMeanL","fvaMeanM","fvaMeanU"],["StrainCondition","Lower Median 95% CI","Median 95% CI","Upper 95% CI","Mean Lower CI 95%","Mean Median CI 95%","Mean Upper CI 95%"]);
    writetable(fvaTable,'output\FVAEstimatedResults.xlsx','Sheet',listRxns{i,1})
end

listRxns{end+1,:} = 'biomass_SC5_notrace';
stoichiometryRxns = [printRxnFormula(model, 'rxnAbbrList',listRxns)];

stoichiometryTable = table(listRxns,stoichiometryRxns);
stoichiometryTable = renamevars(stoichiometryTable,["listRxns","stoichiometryRxns"],["Model Rxn","Reaction Stoichiometry"]);
writetable(stoichiometryTable,'output\FVAReaction.xlsx','Sheet','Stoichiometry');
