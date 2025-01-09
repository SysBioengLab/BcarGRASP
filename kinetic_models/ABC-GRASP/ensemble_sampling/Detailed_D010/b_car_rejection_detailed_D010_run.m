%% Sample and select models using ABC rejection for low growth and detailed structures
% It requires the Excel file located in the GRASP-main/io/input folder as
% input
% It generates an unfiltered version of the ensemble in the same folder, it
% also creates a 'backup' file with all the data of particle simulation for the
% non-reference experimental conditions. This file is later used to filter
% the ensembles from the prior and generate the posterior distribution
% This script must be run from the folder it is placed

%% Clean variables and state matlab standards

clear, close all
rng('default');                 % for reproducibility
format longE
digits(256)                     % To prevent error in calculations

p = gcp('nocreate');            % Deletes the previous parallel sessions if any
if ~isempty(p)
    delete(gcp)
end

%% Path

addpath(fullfile('..','..','.local/matlab'))                                    % Used in Linux
addpath(fullfile('..','..','GRASP-main', 'matlab_code', 'analysisFxns'), ...
        fullfile('..','..','GRASP-main', 'matlab_code', 'ensembleFxns'), ...
        fullfile('..','..','GRASP-main', 'matlab_code', 'patternFxns'));


%% Input

maxNumberOfSamples = 1e5;                                               % Sets the maximum number of samples
eigThreshold = 1e-5;                                                    % Sets the threshold for eigenvalues
        
modelID = 'b_car_rejection_model_detailed_D010';                        % Name of the ensemble model
%diary([modelID,'.txt'])                                                % Starts diary function
diary off

inputFile  = fullfile('..','..','GRASP-main', 'io', 'input', modelID);  % Location of the input file
outputFile = fullfile([modelID, '.mat']);                               % Ouput file

%% Ensemble model

tInitial = tic;
ensemble = buildEnsemble(inputFile, outputFile, maxNumberOfSamples, eigThreshold);  % Generate ensemble
toc(tInitial)