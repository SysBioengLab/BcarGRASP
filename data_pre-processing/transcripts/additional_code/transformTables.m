function tablesCqs = transformTables(filenames,dataFolder)
% This function uses the function transformTable in a 'for' loop to return
% multiple table Cqs in a cell array

% INPUT
% - filenames: names of the Aria reports
% - dataFolder: name of the folder with the Aria reports

% OUTPUT
% - tablesCqs: cell array with multiple tables with Cqs.

for i = 1:length(filenames)
    tablesCqs{i,:} = transformTable(filenames{i,:},dataFolder);
end

end