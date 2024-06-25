%% writeOptoVarsBatch
% calculates descriptive statistics of opto stims and vocalization for many
% recordings contained in a meta data file
clear all
pathToData = '';
metaDataFile = '';
cd(pathToData)
metaData = readtable (metaDataFile);
metaData = table2struct (metaData);

%% Set variables that will not change with loop iterations
breathTraceSampleRate = 1000; %Hz
optoSampRate = 1000;
thresh = 3; %for finding opto peaks
%% Loop through files and calculate variables for selected recording
for i = 1:length(metaData)
    fileName = metaData(i).File;
    txtFile = append (fileName,'.txt');
    breathData = readtable (txtFile,'ReadVariableNames', false);
    breathTrace = table2array (breathData);
    breathTraceSampleRate = 1000; %Hz
    time = 0:1/breathTraceSampleRate:length(breathTrace)/breathTraceSampleRate;
    time = time(1:length(time)-1);
    filtBreathTrace = bandpass (breathTrace,[2,100],1000);
    optoFile = strcat (pathToData,fileName,'_opto.txt');
    optoTrace=readtable(optoFile, 'ReadVariableNames', false);
    optoTrace=table2array(optoTrace);
    optoTrace = downsample (optoTrace, 10);
    optoSampRate = 1000;
    optoTime = 0:1/optoSampRate:length(optoTrace)/optoSampRate;
    optoTime = optoTime(1:length(optoTime)-1);
    optoTrace(optoTrace<5) = 0;
    optoTrace(optoTrace>5) = 5;
    [~, optoStarts] = findpeaks(optoTrace, 'MinPeakHeight', thresh);
    optoStarts = optoStarts';
    [~, stimTrainLocs] = findpeaks(optoTrace, 'MinPeakDistance', 3001, 'MinPeakHeight', thresh);
    UsvData = readtable(strcat(fileName,'_dat'));
    UsvData = table2array(UsvData);
    vocStartTime = UsvData (:,2);
    vocStartTime = vocStartTime';
    noVoc = isempty(vocStartTime);
    vocStartTime = vocStartTime*breathTraceSampleRate;
    vocEndTime = UsvData (:,3);
    vocEndTime = vocEndTime';
    vocEndTime = vocEndTime*breathTraceSampleRate;
    if noVoc == 0
        trimmedVocStartTime = vocStartTime;
        for j = 1:length(vocStartTime)
            if vocStartTime(j)<2000
                trimmedVocStartTime = trimmedVocStartTime (2:end);
                vocEndTime = vocEndTime (2:end);
            elseif vocStartTime(j)>298000
                trimmedVocStartTime = trimmedVocStartTime (1:(length(trimmedVocStartTime)-1));
                vocEndTime = vocEndTime (1:(length(trimmedVocStartTime)));
            end
        end
        vocStartTime = trimmedVocStartTime;
        noVoc = isempty(vocStartTime);
    end
    if noVoc == 0
        if vocStartTime(1)<optoStarts(1)
            vocStartTime = vocStartTime (2:end);
        vocEndTime = vocEndTime (2:end);
        end
    end
    if noVoc == 1
    numcalls = 0;
    interVocInterval = NaN;
    probConVocStims = 0;
    probStimVoc = 0;
    delayOnset = NaN;
elseif noVoc ==0
%identify stims with vocalization
%pull out indices of optoStarts with vocalizations
vocOptoStartInd = interp1 (optoStarts,1:numel(optoStarts), vocStartTime,'previous','extrap');
%use these indices to find stats of vocalization breaths (i.e. value of
%that column in the vector)
vocOptoStarts = optoStarts(vocOptoStartInd);
delayOnset = vocStartTime - vocOptoStarts;
delayOffset = vocEndTime-vocOptoStarts;

%identify multisyllabic vocalizations (same vocExpStart)
%multisyllabindex = or ([diff(vocExpStart) 1] == 0, [1 diff(vocExpStart)] == 0);
%singlesyllabindex = ~multisyllabindex;
[uniqcalls, ~, rnk] = unique(vocOptoStarts);
rnk = transpose (rnk);
numcalls = length (uniqcalls); 
interVocInterval = (vocStartTime(2:end))-(vocEndTime(1:end-1));
%find indicies for opto without vocalization
optoStartInd = [1:length(optoStarts)];
noVocLog = ismember(optoStartInd, vocOptoStartInd);
noVocNum = 1;
for k = 1:length(noVocLog)
    if noVocLog (1,i) == 0
        noVocOptoStartInd(noVocNum)= optoStartInd(k);
        noVocNum = noVocNum +1;
    end
end
%calculate the probability that a voc stim will be followed by a second
%voc stim
consecutiveVocStims = 0;
for q = 1:length(vocOptoStartInd)
    testInd = vocOptoStartInd(q)+1;
    testlog = ismember(vocOptoStartInd,testInd);
    testlog = sum (testlog);
    if testlog == 1
        consecutiveVocStims = consecutiveVocStims+1;
    end
end
probConVocStims = (consecutiveVocStims/length(vocOptoStartInd))*100;
probStimVoc = (length(vocStartTime)/length(optoStarts))*100;
end
save (fullfile(pathToData,strcat(fileName,'optostats.mat')),'numcalls','probStimVoc','probConVocStims','delayOnset','interVocInterval');
end