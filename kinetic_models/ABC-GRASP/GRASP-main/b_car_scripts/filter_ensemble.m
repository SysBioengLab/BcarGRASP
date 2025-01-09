%% Filter an ensemble based on discrepancy scores
% Uses as input the direction fo a 'backup' file used for sampling an
% ensemble with non-reference conditions in the Rejection mode of GRASP
% It also requires the ID of the model, the number of particles to be
% considered from the model, and the number of particles to remain after
% rejection

% The output includes a filtered model located at the GRASP-main/io/output
% direction which contains the defined number of particles to be selected
% as a posterior ensemble, and an unfiltered model located at
% GRASP-main/io/unfiltered_output, with the complete set of model particles
% considered for the rejection process. The ensemble variable is named as
% ensembleFiltered, while the unfiltered 

% The script also performs Control and Response analysis on the filtered
% (posterior) ensemble, and the results are stored in GRASP-main/io/output/cr
% with the same name as the model ID, but with an additional '_crResults'
% at the end of the filename

%% Clean variables and state matlab standards

clear, close all
rng('default');                 % for reproducibility
format longE
digits(256)                     % to prevent errors in calculations

%% Path

addpath(fullfile('..', 'matlab_code', 'analysisFxns'), ...
        fullfile('..', 'matlab_code', 'ensembleFxns'), ...
        fullfile('..', 'matlab_code', 'patternFxns'));

%% Load and Define Inputs
% Modify this section to filter different ensembles

load('../../ensemble_sampling/Detailed_D010/backup_b_car_rejection_model_detailed_D010.mat')    % model to be loaded
modelID = 'b_car_rejection_model_detailed_D010';                                                % model ID
modelsN = 10000;                                                                                % number of particles to be considered for the prior (first N particles)
modelsS = 100;                                                                                  % number of particles to be selected for the posterior

%% Select only the first N models

counter = 0;
for i = 1:length(isModelValidBackup)
    if counter >= modelsN
        isModelValidBackup(i,1) = false;
    end
    if isModelValidBackup(i,1)
        counter = counter + 1;
    end
end

%% Dispose unvalid models

preFilter = isModelValidBackup;

fmins = fminsBackup(preFilter,:);
tolScores = tolScoreBackup(preFilter,:);
times = timesBackup(preFilter,:);

simFluxes = simFluxesBackup(preFilter,:);
xopts = xoptBackup(preFilter,:);
models = modelsBackup(:,preFilter);
strucIdx = structIdxBackup(preFilter,:);

%% Calculate metric

meanTolScores = mean(tolScores,2);

%% Plot Scores

figure()
histogram(meanTolScores,100)

%% Filter

[~,indexesTolScores] = sort(meanTolScores);

filter = indexesTolScores(1:modelsS,:)';

fminsFiltered = fmins(filter,:);
tolScoresFiltered = tolScores(filter,:);
simFluxesFiltered = simFluxes(filter,:);
xoptsFiltered = xopts(filter,:);
modelsFiltered = models(:,filter);
strucIdxFiltered = strucIdx(filter,:);

maxTolScoreFiltered = max(mean(tolScoresFiltered,2));

%% Create Unfiltered Ensemble

ensembleUnfiltered = ensemble;
ensembleUnfiltered.populations(1).strucIdx = strucIdx;                                                                           % model structures
ensembleUnfiltered.populations(1).tolScore = tolScores;                                                                           % tolerance score
ensembleUnfiltered.populations(1).xopt = xopts;                                                                               % optimal value found
ensembleUnfiltered.populations(1).simFluxes = simFluxes;                                                                          % simulated fluxes
ensembleUnfiltered.populations(1).models = models;

%% Create Filtered Ensemble

ensembleFiltered = ensemble;
ensembleFiltered.populations(1).strucIdx = strucIdxFiltered;                                                                           % model structures
ensembleFiltered.populations(1).tolScore = tolScoresFiltered;                                                                           % tolerance score
ensembleFiltered.populations(1).xopt = xoptsFiltered;                                                                               % optimal value found
ensembleFiltered.populations(1).simFluxes = simFluxesFiltered;                                                                          % simulated fluxes
ensembleFiltered.populations(1).models = modelsFiltered;

%% Perform MCA Analysis to check results

t0 = tic;
saveMCAMatrices = 1;
mcaResults = controlAndResponseAnalysis(ensembleFiltered,saveMCAMatrices);
categories = {};       
plotControlAndResponseAnalysis(mcaResults, ensembleFiltered, categories);
t1 = toc(t0);

%% Save Results
outputEnsembleUnfiltered = fullfile('..', 'io', 'unfiltered_output', [modelID, '.mat']);
outputEnsembleFiltered = fullfile('..', 'io', 'output', [modelID, '.mat']);
outputMCA = fullfile('..', 'io', 'cr', [modelID, '_crResults.mat']);
save(outputEnsembleUnfiltered,"ensembleUnfiltered")
save(outputEnsembleFiltered,"ensembleFiltered")
save(outputMCA,"mcaResults")

%% Display worst filtered score:

disp(['The cut-off score was: ',num2str(maxTolScoreFiltered)])








