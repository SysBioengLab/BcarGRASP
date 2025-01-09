function tableCqs = transformTable(filename,dataFolder)
% This function reorganizes the data in an Aria report to a matrix with the
% values of the Cqs obtained

% INPUTS
% - filename: name of the Aria report
% - dataFolder: name of the folder with the Aria report

% OUPUTS
% - tableCqs: matrix with the values of the Cqs


%% Reorganize data

fullfilename = fullfile(dataFolder,filename);

originalTable = readtable(fullfilename,'Sheet','Tabular Results');

rowNames = ['A','B','C','D','E','F','G','H'];

% wellNames = originalTable{:,{'Well'}};
% wellCqs = originalTable{:,{'Cq__R__T__'}};

wellNames = originalTable{:,1};
wellCqs = originalTable{:,13};

tableCqs = NaN(8,12);

for i = 1:length(wellNames)
    wellRow = find(rowNames == wellNames{i,:}(1));
    wellColumn = str2num(wellNames{i,:}(2:end));
    tableCqs(wellRow,wellColumn) = wellCqs(i);
end

end
    



