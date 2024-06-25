%% batchProcessStats.m
%% Open up recordings of breathing and vocalization and save summary statistics
%Create metadata file in excel prior to running
%% Load metadata and extract filenames
clear all
close all
pathToData = 'C:\Users\Alastair\Box\USV Behavior\Urine recordings\analyze\';
metaDataFile = [pathToData,'metadata.xlsx'];
metaData = readtable (metaDataFile,'ReadVariableNames',true);
metaData = table2struct(metaData);
for i = 1:length(metaData)
    tempFile = [metaData(i).bl '_callstats.mat'];
    if isfile (tempFile) == 0
        writeUSVvals(metaData(i).bl, pathToData);
    end
    tempFile = [metaData(i).u1 '_callstats.mat'];
    if isfile (tempFile) == 0
        writeUSVvals(metaData(i).u1, pathToData);
    end
    tempFile = [metaData(i).u2 '_callstats.mat'];
    if isfile (tempFile) == 0
        writeUSVvals(metaData(i).u2, pathToData);
    end
    tempFile = [metaData(i).u3 '_callstats.mat'];
    if isfile (tempFile) == 0
        writeUSVvals(metaData(i).u3, pathToData);
    end
end