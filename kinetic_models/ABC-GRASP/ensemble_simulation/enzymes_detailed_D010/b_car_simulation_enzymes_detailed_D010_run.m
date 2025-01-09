%% Simulation of detailed models generated at low growth rate when enzyme abundances are perturbated
% It uses the filtered ensemble (generated using the script
% 'filter_ensemble') at the GRASP-main/io/output direction to simulate
% model particles
% It also requires an Excel file at GRASP-main/io/input with the quantities
% that are being perturbated
% The output is a 'backup' file, with the details of each simulated model
% particle
% The script must be run from the folder it is located

%% Clean variables and state matlab standards

clear, close all
rng('default');                 % for reproducibility
format longE
digits(256)                     % to prevent errors in calculations

p = gcp('nocreate');            % deletes previous parallel session if any
if ~isempty(p)
    delete(gcp)
end

%% Path

addpath(fullfile('..','..','.local/matlab') )                                 % in Linux
addpath(fullfile('..','..','GRASP-main', 'matlab_code', 'analysisFxns'), ...
        fullfile('..','..','GRASP-main', 'matlab_code', 'ensembleFxns'), ...
        fullfile('..','..','GRASP-main', 'matlab_code', 'patternFxns'))

%% Input Parameters

ensembleDir = fullfile('..','..','GRASP-main','io','output','b_car_rejection_model_detailed_D010.mat');     % direction of the model ensemble
strucIdx = 1;                                                                                               % always only one structure
dataDir = fullfile('..','..','GRASP-main','io','input','b_car_enzymes_simulation.xlsx');                    % direction of the input file with perturbances
conditionsN = 32;                                                                                           % number of conditions to be simulated
parallelSampling = true;                                                                                    % if parallel process wants to be used
coresN = 16;                                                                                                % number of cores for parallel process
backupDir = 'backup_b_car_simulation_enzymes_model_detailed_D010.mat';                                      % direction and name of backup file
solver = 'NLOPT'; %FMINCON                                                                                  % solver type, either NLOPT or FMINCON

%% Perform Simulation

tInitial = tic;
simulateFluxes(ensembleDir,strucIdx,dataDir,conditionsN,parallelSampling,coresN,solver,backupDir)           % perform simulations
toc(tInitial)
