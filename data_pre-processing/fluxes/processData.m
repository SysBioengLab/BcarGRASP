%% Process fluxes data: EX_lyc and EX_b_car
% This script processed several data sources to estimate the fluxes EX_lyc
% and EX_b_car of the kinetic models.

% It requires the names of the files to be used as input
% - BiomassConversion: Measurements to estimate biomass/volume conversion 
% from OD600
% - BiomassMeasurements: Measurements of OD600 for chemostat cultivations
% - RPMFluxConversion: Measurements for conversion from RPM to feed flow
% - RPMMeasurements: RPM of the feed flow for the chemostat cultivations
% - MetabolitesMeasurements: Concentration of metabolites in the chemostat
% cultivations (per volume)
% - MolWeights: Molecular weight of metabolites
% - DensityMeanStd: Assumed mean density and its std of yeast cells
% - MoistWeightConversion: Assumed mean moisture content and its std of
% yeast cells

% It produces several outputs, all of them stored in tables at 'output'
% filefolder.
% - BiomassConversionStatistics: mean conversion factors (biomass/L)/OD600
% for each of the strains
% - BiomassConversionsStatistics: same file as before but for chemostats
% - BiomassMeasurementsOutliers: outliers of biomass measurements for
% chemostats
% - BiomassMeasurementsStatistics: mean OD600 of chemostats
% - DilutionRatesStatistics: mean dilution rate of chemostats
% - MetablitesProductionStatistics: titers of metabolites in the chemostat
% cultivations
% - FluxesBiomassStatistics: mean fluxes of metabolites per biomass
% - FluxesVolumeStatisticss: mean fluxes of metabolites per intracellular
% volume

%% Clean variables and close windows
clear, close all
rng('default');                 % for reproducibility
format longE

%% Inputs: Import Data or state Filenames

strainConditions = {'2D01';'2D025';'3D01';'3D025';'4D01';'4D025'};                                                      % names of chemostat cultivations
metabolitesCulture = {'Glucose';'Ethanol';'Acetate';'Glycerol';'Lycopene';'BetaCarotene'};                              % names of metabolites measured in the chemostats
metabolitesMedia = {'GlucoseMedia'};                                                                                    % names of metabolites measured in the feed media
BiomassConversionTable = readtable('..\..\biomass_and_metabolites_data\BiomassConversion.xlsx','Sheet','Data');
BiomassMeasurementsFile = '..\..\biomass_and_metabolites_data\BiomassMeasurements.xlsx';
RPMFluxConversionTable = readtable('..\..\biomass_and_metabolites_data\RPMFluxConversion.xlsx','Sheet','Data');
RPMMeasurementsTable = readtable('..\..\biomass_and_metabolites_data\RPMMeasurements.xlsx','Sheet','Data');
MetabolitesMeasurementsFile = '..\..\biomass_and_metabolites_data\MetabolitesMeasurements.xlsx';
MolWeightsTable = readtable('..\..\biomass_and_metabolites_data\MolWeights.xlsx','Sheet','Data');
DensityTable = readtable('..\..\biomass_and_metabolites_data\DensityMeanStd.xlsx','Sheet','Data');
MoistWeightConversionTable = readtable('..\..\biomass_and_metabolites_data\MoistWeightConversion.xlsx','Sheet','Data');

% Define number of data to be simulated

simulationN = 1e7;                                                          % Samples for Monte Carlo simulation
alpha = 0.05;                                                               % alpha for confidenc intervals


%% Process Biomass Conversion from OD600 to gDCW/L

biomassConversionUnprocessed = BiomassConversionTable{:,2:end}';

% Calculate average of each biological replicate
biomassConversionBiologicalReplicatesN = max(biomassConversionUnprocessed(4,:));

biomassConversion = zeros(3,biomassConversionBiologicalReplicatesN);
biomassConversionUnprocessedOutliers = false(3,biomassConversionBiologicalReplicatesN);
for ix = 1:biomassConversionBiologicalReplicatesN
    tempBiomassConversion = biomassConversionUnprocessed(1:3,biomassConversionUnprocessed(4,:)==ix);
%     tempBiomassConversionOutliers = isoutlier(tempBiomassConversion,'median',2);
%     biomassConversionUnprocessedOutliers(1:3,biomassConversionUnprocessed(4,:)==ix) = tempBiomassConversionOutliers;
%     tempBiomassConversion(tempBiomassConversionOutliers) = nan;
    biomassConversion(:,ix) = mean(tempBiomassConversion,2,'omitnan');
end

% Check for Outliers
biomassConversionOutliers = isoutlier(biomassConversion,'median',2);
biomassConversionFiltered = biomassConversion;
biomassConversionFiltered(biomassConversionOutliers) = nan;

% Calculate Statistics
biomassConversionMean = mean(biomassConversionFiltered,2,'omitnan');
biomassConversionSD = std(biomassConversionFiltered,0,2,'omitnan');
biomassConversionN = sum(~biomassConversionOutliers,2);
biomassConversionSE = biomassConversionSD./sqrt(biomassConversionN);
biomassConversionDF = biomassConversionN-1;

% Export Table
biomassConversionHeaders = BiomassConversionTable.Properties.VariableNames(:,2:end-1)';
biomassConversionStatisticsTable = table(biomassConversionHeaders,biomassConversionMean,biomassConversionSD,biomassConversionSE,biomassConversionN,biomassConversionDF);
biomassConversionStatisticsTable = renamevars(biomassConversionStatisticsTable,["biomassConversionHeaders","biomassConversionMean","biomassConversionSD","biomassConversionSE","biomassConversionN","biomassConversionDF"],["Strain","Mean [(gDCW/L)/OD600]","SD [(gDCW/L)/OD600]","SE [(gDCW/L)/OD600]","N Measurements","Freedoom Degrees"]);
writetable(biomassConversionStatisticsTable,"output\BiomassConversionStatistics.xlsx",'Sheet','Statistics');

% Redistribute the information to match the conditions
biomassConversionMean = reshape(repmat(biomassConversionMean',2,1),[],1);
biomassConversionSD = reshape(repmat(biomassConversionSD',2,1),[],1);
biomassConversionN = reshape(repmat(biomassConversionN',2,1),[],1);
biomassConversionSE = reshape(repmat(biomassConversionSE',2,1),[],1);
biomassConversionDF = reshape(repmat(biomassConversionDF',2,1),[],1);

biomassConversionStatisticsTable = table(strainConditions,biomassConversionMean,biomassConversionSD,biomassConversionSE,biomassConversionN,biomassConversionDF);
biomassConversionStatisticsTable = renamevars(biomassConversionStatisticsTable,["strainConditions","biomassConversionMean","biomassConversionSD","biomassConversionSE","biomassConversionN","biomassConversionDF"],["StrainCondition","Mean [(gDCW/L)/OD600]","SD [(gDCW/L)/OD600]","SE [(gDCW/L)/OD600]","N Measurements","Freedoom Degrees"]);
writetable(biomassConversionStatisticsTable,"output\BiomassConversionsStatistics.xlsx",'Sheet','Statistics');



%% Process Biomass Measurements in OD600

% Define number of conditions
strainConditionsN = length(strainConditions);

% Extract information in tables
BiomassMeasurementsTables = cell(strainConditionsN,1);
for i = 1:strainConditionsN
    BiomassMeasurementsTables{i,1} = readtable(BiomassMeasurementsFile,'Sheet',strainConditions{i,1});
end

% Define vectors for outliers, mean, sd, se, and n
biomassMeasurementsOutliers = cell(strainConditionsN,1);
biomassMeasurementsMean = zeros(strainConditionsN,1);
biomassMeasurementsSD = zeros(strainConditionsN,1);
biomassMeasurementsN = zeros(strainConditionsN,1);
biomassMeasurementsSE = zeros(strainConditionsN,1);
biomassMeasurementsDF = zeros(strainConditionsN,1);

% Calculate Statistics and Save Tables with Outliers for recognition of
% erased measurements
for i = 1:strainConditionsN
    measurements = BiomassMeasurementsTables{i,1}{:,2};
    biomassMeasurementsOutliers{i,1} = isoutlier(measurements,'median');
    measurementsFiltered = measurements;
    measurementsFiltered(biomassMeasurementsOutliers{i,1}) = nan;
    biomassMeasurementsMean(i,1) = mean(measurementsFiltered,'omitnan');
    biomassMeasurementsSD(i,1) = std(measurementsFiltered,0,'omitnan');
    biomassMeasurementsN(i,1) = sum(~biomassMeasurementsOutliers{i,1});
    biomassMeasurementsSE(i,1) = biomassMeasurementsSD(i,1)./sqrt(biomassMeasurementsN(i,1));
    biomassMeasurementsDF(i,1) = biomassMeasurementsN(i,1)-1;
    BiomassMeasurementsTables{i,1}.Outliers = biomassMeasurementsOutliers{i,1};
    writetable(BiomassMeasurementsTables{i,1},'output\BiomassMeasurementsOutliers.xlsx','Sheet',strainConditions{i,1})
end

% Export Table
biomassMeasurementsStatisticsTable = table(strainConditions,biomassMeasurementsMean,biomassMeasurementsSD,biomassMeasurementsSE,biomassMeasurementsN,biomassMeasurementsDF);
biomassMeasurementsStatisticsTable = renamevars(biomassMeasurementsStatisticsTable,["strainConditions","biomassMeasurementsMean","biomassMeasurementsSD","biomassMeasurementsSE","biomassMeasurementsN","biomassMeasurementsDF"],["StrainCondition","Mean [OD600]","SD [OD600]","SE [OD600]","N Measurements","Freedom Degrees"]);
writetable(biomassMeasurementsStatisticsTable,"output\BiomassMeasurementsStatistics.xlsx",'Sheet','Statistics');

%% Process Dilution Rates (Version considering Intercept as zero)
% 
% % Calculate flux in ml/h
% rpmConversion = RPMFluxConversionTable{:,2};
% volumeConversion = RPMFluxConversionTable{:,3};
% timeConversion = RPMFluxConversionTable{:,4};
% fluxConversion = volumeConversion./timeConversion*60;
% 
% % Perform linear regression
% RPMtoFluxLM = fitlm(rpmConversion,fluxConversion,'Intercept',false);
% rpmMean = mean(rpmConversion);
% rpmN = length(rpmConversion);
% 
% % Plot linear regression
% RPMtoFluxFigure = figure();
% plot(RPMtoFluxLM);
% hold on
% ylabel('Flux [ml/h]')
% xlabel('RPM Pump')
% title('Conversion from RPM to Pump Flux')
% 
% % Calculate Media Flux based on RPM
% rpmMeasurements = RPMMeasurementsTable{:,2:end}';
% [mediaFluxMean,mediaFluxCI] = predict(RPMtoFluxLM,rpmMeasurements);
% mediaFluxSE = (mediaFluxMean - mediaFluxCI(:,1))./tinv(0.975,rpmN-1);
% 
% % Add prediction intervals to the plot
% plotConversionRPM = (0:0.1:100)';
% [plotConversionFluxMean,plotConversionFluxCI] = predict(RPMtoFluxLM,plotConversionRPM);
% plot(plotConversionRPM,plotConversionFluxCI,'k--','DisplayName','Prediction Interval 95%')
% 
% % Calculate Dilution Rate (Assuming Volume of 500ml)
% dilutionRatesMean = mediaFluxMean/500;
% dilutionRatesSE = mediaFluxSE/500;
% dilutionRatesN = repmat(rpmN,strainConditionsN,1);
% dilutionRatesDF = dilutionRatesN-1;
% 
% % Export Table
% dilutionRatesTable = table(strainConditions,dilutionRatesMean,dilutionRatesSE,dilutionRatesN,dilutionRatesDF);
% dilutionRatesTable = renamevars(dilutionRatesTable,["strainConditions","dilutionRatesMean","dilutionRatesSE","dilutionRatesN","dilutionRatesDF"],["strainConditions","Mean [1/h]","SE [1/h]","N Measurements","Freedom Degrees"]);
% writetable(dilutionRatesTable,'output\DilutionRatesStatistics.xlsx','Sheet','Statistics');


%% Process Dilution Rates (Version Calculating the intercept):

% Calculate flux in ml/h
rpmConversion = RPMFluxConversionTable{1:8,2};
volumeConversion = RPMFluxConversionTable{1:8,3};
timeConversion = RPMFluxConversionTable{1:8,4};
fluxConversion = volumeConversion./timeConversion*60;

% Perform linear regression
RPMtoFluxLM = fitlm(rpmConversion,fluxConversion);
rpmMean = mean(rpmConversion);
rpmN = length(rpmConversion);

% Plot linear regression
RPMtoFluxFigure = figure();
plot(RPMtoFluxLM);
hold on
ylabel('Flux [ml/h]')
xlabel('RPM Pump')
title('Conversion from RPM to Pump Flux')

% Calculate Media Flux based on RPM
rpmMeasurements = RPMMeasurementsTable{:,2:end}';
[mediaFluxMean,mediaFluxCI] = predict(RPMtoFluxLM,rpmMeasurements);
mediaFluxSE = (mediaFluxMean - mediaFluxCI(:,1))./tinv(0.975,rpmN-2);

% Add prediction intervals to the plot
plotConversionRPM = (0:0.1:20)';
[plotConversionFluxMean,plotConversionFluxCI] = predict(RPMtoFluxLM,plotConversionRPM);
plot(plotConversionRPM,plotConversionFluxCI,'k--','DisplayName','Prediction Interval 95%')

% Calculate Dilution Rate (Assuming Volume of 500ml)
dilutionRatesMean = mediaFluxMean/500;
dilutionRatesSE = mediaFluxSE/500;
dilutionRatesN = repmat(rpmN,strainConditionsN,1);
dilutionRatesDF = dilutionRatesN-2;

% Export Table
dilutionRatesTable = table(strainConditions,dilutionRatesMean,dilutionRatesSE,dilutionRatesN,dilutionRatesDF);
dilutionRatesTable = renamevars(dilutionRatesTable,["strainConditions","dilutionRatesMean","dilutionRatesSE","dilutionRatesN","dilutionRatesDF"],["strainConditions","Mean [1/h]","SE [1/h]","N Measurements","Freedom Degrees"]);
writetable(dilutionRatesTable,'output\DilutionRatesStatistics.xlsx','Sheet','Statistics');

%% Process Metabolites Measurements

metabolitesCultureN = length(metabolitesCulture);
metabolitesMediaN = length(metabolitesMedia);
MetabolitesMeasurementsCultureTable = cell(metabolitesCultureN,1);
MetabolitesMeasurementsMediaTable = cell(metabolitesMediaN,1);

for i = 1:metabolitesCultureN
    MetabolitesMeasurementsCultureTable{i,1} = readtable(MetabolitesMeasurementsFile,'Sheet',metabolitesCulture{i,1});
end
for i = 1:metabolitesMediaN
    MetabolitesMeasurementsMediaTable{i,1} = readtable(MetabolitesMeasurementsFile,'Sheet',metabolitesMedia{i,1});
end

metabolitesCultureOutliers = cell(metabolitesCultureN,1);
metabolitesMediaOutliers = cell(metabolitesMediaN,1);

metabolitesProductionStatistics = cell(metabolitesCultureN,1);

metabolitesProductionMean = zeros(strainConditionsN,metabolitesCultureN);
metabolitesProductionSD = zeros(strainConditionsN,metabolitesCultureN);
metabolitesProductionN = zeros(strainConditionsN,metabolitesCultureN);
metabolitesProductionSE = zeros(strainConditionsN,metabolitesCultureN);
metabolitesProductionDF = zeros(strainConditionsN,metabolitesCultureN);

for i = 1:metabolitesMediaN
    measurementsMedia = MetabolitesMeasurementsMediaTable{i,1}{:,2:end};
    measurementsCulture = MetabolitesMeasurementsCultureTable{i,1}{:,2:end};
    metabolitesMediaOutliers{i,1} = isoutlier(measurementsMedia,'median',1);
    metabolitesCultureOutliers{i,1} = isoutlier(measurementsCulture,'median',1);

    if i ==1
        metabolitesMediaOutliers{i,1}(1,4) = true;
        metabolitesCultureOutliers{i,1}(3,4) = true;
    end

    measurementsMediaFiltered = measurementsMedia;
    measurementsCultureFiltered = measurementsCulture;

    measurementsMediaFiltered(metabolitesMediaOutliers{i,1}) = nan;
    measurementsCultureFiltered(metabolitesCultureOutliers{i,1}) = nan;
    
    measurementsMean = (mean(measurementsCultureFiltered,1,'omitnan')-mean(measurementsMediaFiltered,1,'omitnan'))';
    [~,~,~,stats] = ttest2(measurementsCultureFiltered,measurementsMediaFiltered,'Vartype','equal');
    measurementsSD = stats.sd';
    measurementsN = (sum(~metabolitesMediaOutliers{i,1})+sum(~metabolitesCultureOutliers{i,1}))';
    measurementsDF = stats.df';
    measurementsSE = (stats.sd.*sqrt(1./sum(~metabolitesMediaOutliers{i,1})+1./sum(~metabolitesCultureOutliers{i,1})))';
    
    metabolitesProductionMean(:,i) = measurementsMean;
    metabolitesProductionSD(:,i) = measurementsSD;
    metabolitesProductionN(:,i) = measurementsN;
    metabolitesProductionSE(:,i) = measurementsSE;  
    metabolitesProductionDF(:,i) = measurementsDF;

    measurementsTable = table(strainConditions,measurementsMean,measurementsSD,measurementsSE,measurementsN,measurementsDF);
    measurementsTable = renamevars(measurementsTable,["strainConditions","measurementsMean","measurementsSD","measurementsSE","measurementsN","measurementsDF"],["StrainCondition","Mean [g/L]","SD [g/L]","SE [g/L]","N Measurements","Freedom Degrees"]);
    metabolitesProductionStatistics{i,1} = measurementsTable;
    writetable(measurementsTable,"output\MetabolitesProductionStatistics.xlsx","Sheet",metabolitesCulture{i,1})
end

for i = (metabolitesMediaN+1):metabolitesCultureN
    measurements = MetabolitesMeasurementsCultureTable{i,1}{:,2:end};
    metabolitesCultureOutliers{i,1} = isoutlier(measurements,'median',1);

    measurementsFiltered = measurements;
    measurementsFiltered(metabolitesCultureOutliers{i,1}) = nan;
    
    measurementsMean = mean(measurementsFiltered,1,'omitnan')';
    measurementsSD = std(measurementsFiltered,0,1,'omitnan')';
    measurementsN = sum(~metabolitesCultureOutliers{i,1},1)';
    measurementsSE = measurementsSD./sqrt(measurementsN);
    measurementsDF = measurementsN-1;
    
    metabolitesProductionMean(:,i) = measurementsMean;
    metabolitesProductionSD(:,i) = measurementsSD;
    metabolitesProductionN(:,i) = measurementsN;
    metabolitesProductionSE(:,i) = measurementsSE;
    metabolitesProductionDF(:,i) = measurementsDF;

    measurementsTable = table(strainConditions,measurementsMean,measurementsSD,measurementsSE,measurementsN,measurementsDF);
    if i < 5
        measurementsTable = renamevars(measurementsTable,["strainConditions","measurementsMean","measurementsSD","measurementsSE","measurementsN","measurementsDF"],["StrainCondition","Mean [g/L]","SD [g/L]","SE [g/L]","N Measurements","Freedom Degrees"]);
    else
        measurementsTable = renamevars(measurementsTable,["strainConditions","measurementsMean","measurementsSD","measurementsSE","measurementsN","measurementsDF"],["StrainCondition","Mean [mg/L]","SD [mg/L]","SE [mg/L]","N Measurements","Freedom Degrees"]);
    end
    metabolitesProductionStatistics{i,1} = measurementsTable;
    writetable(measurementsTable,"output\MetabolitesProductionStatistics.xlsx","Sheet",metabolitesCulture{i,1})
end

%% Calculate metabolic Fluxes using Simulation

metabolitesFluxesDistributions = cell(metabolitesCultureN,1);

% Extract Molecular Weights
metabolitesMolWeightNames = MolWeightsTable{:,1};
metabolitesMolWeights = MolWeightsTable{:,2};

% Obtain parameters and their distributions for volume fluxes
densityMean = DensityTable{1,1};
densitySD = DensityTable{1,2};
moistWeightConversionMean = MoistWeightConversionTable{1,1};
moistWeightConversionSD = MoistWeightConversionTable{1,2};

densityDistribution = normrnd(densityMean,densitySD,1,simulationN);
moistWeightConversionDistribution = normrnd(moistWeightConversionMean,moistWeightConversionSD,1,simulationN);

% Define arrays to store values
metabolitesBiomassFluxesMean = zeros(strainConditionsN,metabolitesCultureN);
metabolitesBiomassFluxesSE = zeros(strainConditionsN,metabolitesCultureN);
metabolitesBiomassFluxesLCI = zeros(strainConditionsN,metabolitesCultureN);
metabolitesBiomassFluxesUCI = zeros(strainConditionsN,metabolitesCultureN);

metabolitesVolumeFluxesMean = zeros(strainConditionsN,metabolitesCultureN);
metabolitesVolumeFluxesSE = zeros(strainConditionsN,metabolitesCultureN);
metabolitesVolumeFluxesLCI = zeros(strainConditionsN,metabolitesCultureN);
metabolitesVolumeFluxesUCI = zeros(strainConditionsN,metabolitesCultureN);
metabolitesVolumeFluxesMeanGRASP = zeros(strainConditionsN,metabolitesCultureN);
metabolitesVolumeFluxesStdGRASP = zeros(strainConditionsN,metabolitesCultureN);

% Calculate the distribution of fluxes and their statistics, then store
% them
for i = 1:strainConditionsN
    biomassMeasurementDistribution = biomassMeasurementsMean(i,:) + trnd(biomassMeasurementsDF(i,:),1,simulationN).*biomassMeasurementsSE(i,:);
    biomassConversionDistribution = biomassConversionMean(i,:) + trnd(biomassConversionDF(i,:),1,simulationN).*biomassConversionSE(i,:);
    dilutionRateDistribution = dilutionRatesMean(i,:) + trnd(dilutionRatesDF(i,:),1,simulationN).*dilutionRatesSE(i,:);
    for j = 1:metabolitesCultureN
        metaboliteMolWeight = metabolitesMolWeights(strcmp(metabolitesCulture{j,:},metabolitesMolWeightNames),:);
        metaboliteProductionDistribution = metabolitesProductionMean(i,j)+trnd(metabolitesProductionDF(i,j),1,simulationN).*metabolitesProductionSE(i,j);
        
        metaboliteFluxDistribution = metaboliteProductionDistribution*1000./metaboliteMolWeight./biomassMeasurementDistribution./biomassConversionDistribution.*dilutionRateDistribution;
        metabolitesBiomassFluxesMean(i,j) = mean(metaboliteFluxDistribution);
        metabolitesBiomassFluxesSE(i,j) = std(metaboliteFluxDistribution,1);
        metabolitesBiomassFluxesLCI(i,j) = prctile(metaboliteFluxDistribution,100*alpha/2);
        metabolitesBiomassFluxesUCI(i,j) = prctile(metaboliteFluxDistribution,100*(1-alpha/2));

        metaboliteVolumeFluxDistribution = metaboliteFluxDistribution.*densityDistribution*1000./moistWeightConversionDistribution;
        metabolitesVolumeFluxesMean(i,j) = mean(metaboliteVolumeFluxDistribution);
        metabolitesVolumeFluxesSE(i,j) = std(metaboliteVolumeFluxDistribution,1);
        metabolitesVolumeFluxesLCI(i,j) = prctile(metaboliteVolumeFluxDistribution,100*alpha/2);
        metabolitesVolumeFluxesUCI(i,j) = prctile(metaboliteVolumeFluxDistribution,100*(1-alpha/2));
        metabolitesVolumeFluxesMeanGRASP(i,j) = (metabolitesVolumeFluxesUCI(i,j)+metabolitesVolumeFluxesLCI(i,j))/2;
        metabolitesVolumeFluxesStdGRASP(i,j) = (metabolitesVolumeFluxesUCI(i,j)-metabolitesVolumeFluxesLCI(i,j))/4;
    end
end

% Export data to tables
for i = 1:metabolitesCultureN
    fluxBiomassMean = metabolitesBiomassFluxesMean(:,i);
    fluxBiomassSE = metabolitesBiomassFluxesSE(:,i);
    fluxBiomassLCI = metabolitesBiomassFluxesLCI(:,i);
    fluxBiomassUCI = metabolitesBiomassFluxesUCI(:,i);

    fluxVolumeMean = metabolitesVolumeFluxesMean(:,i);
    fluxVolumeSE = metabolitesVolumeFluxesSE(:,i);
    fluxVolumeLCI = metabolitesVolumeFluxesLCI(:,i);
    fluxVolumeUCI = metabolitesVolumeFluxesUCI(:,i);
    fluxVolumeMeanGRASP = metabolitesVolumeFluxesMeanGRASP(:,i);
    fluxVolumeStdGRASP = metabolitesVolumeFluxesStdGRASP(:,i);

    fluxBiomassTable = table(strainConditions,fluxBiomassMean,fluxBiomassSE,fluxBiomassLCI,fluxBiomassUCI);
    fluxBiomassTable = renamevars(fluxBiomassTable,["strainConditions","fluxBiomassMean","fluxBiomassSE","fluxBiomassLCI","fluxBiomassUCI"],["StrainCondition","Mean [mmol/gDCWh]","SE [mmol/gDCWh]","Lower CI 95% [mmol/gDCWh]","Upper CI 95% [mmol/gDCWh]"]);
    writetable(fluxBiomassTable,'output\FluxesBiomassStatistics.xlsx','Sheet',metabolitesCulture{i,1});

    fluxVolumeTable = table(strainConditions,fluxVolumeMean,fluxVolumeSE,fluxVolumeLCI,fluxVolumeUCI,fluxVolumeMeanGRASP,fluxVolumeStdGRASP);
    fluxVolumeTable = renamevars(fluxVolumeTable,["strainConditions","fluxVolumeMean","fluxVolumeSE","fluxVolumeLCI","fluxVolumeUCI","fluxVolumeMeanGRASP","fluxVolumeStdGRASP"],["StrainCondition","Mean [mmol/Lh]","SE [mmol/Lh]","Lower 95% CI [mmol/Lh]","Upper 95% CI [mmol/Lh]","Mean for GRASP  [mmol/Lh]","SD for GRASP [mmol/Lh]"]);
    writetable(fluxVolumeTable,'output\FluxesVolumeStatistics.xlsx','Sheet',metabolitesCulture{i,1});
end
















