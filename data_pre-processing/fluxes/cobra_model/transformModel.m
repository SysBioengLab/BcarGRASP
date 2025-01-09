%% This script updates the model for only Bcar producing strains
% Does not include genes to form beta-ionone using beta-carotene
% Based in the code of Pedro Saa presented in Lopez et al. (2020)
% It employs the iMM904 GEM contextualized for the study as input
% Requires installation of the Cobra Toolbox (Heirendt et al. 2019)
% The output is the the re-contextualized constraint-based model 
% IMM904_b_car, and it is stored in the same folder as the script

% It also generates excels with the reactions of the models to review them
% before using other scripts

%% Clean variables and console, close windows
clc,clearvars,close all

%% Load data and initialize variables
% Initialize CobraToolbox
initCobraToolbox

% Load GEM of the beta-ionone producing yeast
load('iMM904_contextualized.mat')

% Write the model as excel to simplify visualization
if usejava('desktop') % This line of code is to avoid execution of example in non gui-environments
writeCbModel(model,'fileName', 'iMM904_b_ionone.xls','format','xls')
end

% Show new reactions added (25 in total)
model.rxnNames(end-24:end)
model.rxns(end-24:end)
printRxnFormula(model, 'rxnAbbrList',model.rxns(end-24:end));
pre_indexes = [23	22	17	16	15	14	8	7	6	5	4	3	2	1	0];
indexes = size(model.rxns,1) - pre_indexes;

new_model = removeRxns(model, model.rxns(indexes));

% Show reactions actually included (only 10)
new_model.rxnNames(end-9:end)
new_model.rxns(end-9:end)
printRxnFormula(new_model, 'rxnAbbrList',new_model.rxns(end-9:end));

if usejava('desktop') % This line of code is to avoid execution of example in non gui-environments
writeCbModel(new_model,'fileName', 'iMM904_b_car.xls','format','xls')
writeCbModel(new_model,'iMM904_b_car.mat')
end




