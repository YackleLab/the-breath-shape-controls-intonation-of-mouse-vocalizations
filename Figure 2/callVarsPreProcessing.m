%% callVarsPreProcessing.m
% Prepare a recording values file to be compiled into the large dataset
% prior to this run vocal mat to detect and classify USVs and
% peakFreqWriter.m to vectorize the USV spectrogram
%% Set up file names and locations
fileName = ['']; %specify file name
pathToData = ''; %specify path to data
cd(pathToData)
breathFile = [fileName '.txt'];
vocFile = [fileName '.xlsx'];
freqFile = [fileName 'peakFreq.mat'];
%% Load in the breathing data, filter, create time vec
breathData = readtable(breathFile,'ReadVariableNames',false);
breathTrace = table2array(breathData);
sampRate = 1000; %Hz
time = 0:1/sampRate:length(breathTrace)/sampRate;
time = time(1:length(time)-1);;
filtBreathTrace = bandpass(breathTrace,[2 35],sampRate);
%% Calculate breath parameters
% This section depends on functions from Bachmutsky et al 2020
% available on https://github.com/YackleLab/Opioids-depress-breathing-through-two-small-brainstem-sites
durThresh = 0.00002; %minimum duration of insp
inspAmpThresh = -0.015; %minimum insp amplitude
expAmpThresh = 0; %minimum exp amplitude
%find starts and trim out filtering artefact at start & end of recording
breathStarts = getbreathstarts(filtBreathTrace, durThresh, inspAmpThresh, expAmpThresh);
breathStarts = breathStarts(2:end);
breathTrimInd = breathStarts > 1500;
breathStarts = breathStarts(breathTrimInd);
breathTrimInd = breathStarts < 298500;
breathStarts = breathStarts(breathTrimInd);
%pull out individual waveforms
breathMat = cell(1,length(breathStarts)-1);
for breathInd = 1:length(breathStarts)-1;
    tempStart = breathStarts(breathInd);
    tempEnd = breathStarts(breathInd+1);
    tempBreath = filtBreathTrace(tempStart:tempEnd);
    breathMat{1,breathInd} = tempBreath;
end
[~,~,inspDur,expDur,~,~,~] = getbreathvals(breathMat);
inspStarts = breathStarts(1:length(breathStarts)-1);
inspDur = transpose(inspDur);
expDur = transpose(expDur);
expStarts = inspStarts +inspDur;

%% Open USV data
vocQuery = isfile(vocFile);
if vocQuery == 1
    vocData = readtable(vocFile,'VariableNamingRule','preserve');
    vocStarts = (transpose(vocData.Start_time))*1000;
    vocEnds = (transpose(vocData.End_time))*1000;
    vocClass = transpose(vocData.Class);
    noVm = isempty(vocStarts);
else noVm = 1;
    disp ('No VocalMat file found')
end
vocClass = string(vocClass);
load(freqFile);
%% trim vocs that occur during filtering artefact
vocTrimInd = vocStarts > 1500;
vocStarts = vocStarts(vocTrimInd);
vocEnds = vocEnds(vocTrimInd);
vocClass = vocClass(vocTrimInd);
vocTrimInd = vocStarts < 298500;
vocStarts = vocStarts(vocTrimInd);
vocEnds = vocEnds(vocTrimInd);
vocClass = vocClass(vocTrimInd);
% remove starts,ends and class for vocs classified as 'noise'
for vocNum = 1:length(vocClass)
    tempClass = vocClass {1,vocNum};
    if contains(tempClass,'noise_dist') == 1
        vocInd(vocNum) = 0;
    else vocInd(vocNum) =1;
    end
end
vocInd = logical(vocInd);
vocClass = vocClass(vocInd);
vocStarts = vocStarts(vocInd);
vocEnds = vocEnds(vocInd);
%% Aligning breath parameters and vocalization times
% check vocalization doesn't occur prior to first inspiration and
% discard if so
if vocStarts(1)<inspStarts(1)
    vocStarts = vocStarts(2:end);
    vocEnds = vocEnds(2:end);
    vocClass = vocClass(1,2:end);
end
%find closest inspStarts value for vocStarts and index
vocInspStartInd = interp1(inspStarts,1:numel(inspStarts),vocStarts,'previous','extrap');
%find features of these breaths
vocExpStarts = expStarts(vocInspStartInd);
vocInspStarts = inspStarts(vocInspStartInd);
vocExpDur = expDur(vocInspStartInd);
vocMat = breathMat(vocInspStartInd);
delayOnsetExp = vocStarts- vocExpStarts;
delayOffsetExp = vocEnds - vocExpStarts;
normDelayOnsetExp = delayOnsetExp./vocExpDur;
normDelayOffsetExp =  delayOffsetExp./vocExpDur;

%% Calculate correlation coefficient of airflow and USV frequency
spectSampRate = 2000; % 1/spectrogram bin width (0.5 ms)
upSampleFactor = round(spectSampRate/sampRate);
%interpolate breathtrace values to match sample rate of spect
upSampBT = interp(smoothBreathTrace,upSampleFactor);
%trim first and last values to match times with spect
upSampBT = upSampBT(2:end-1);
%create time vec
upTime = 0:1/spectSampRate:length(upSampBT)/spectSampRate;
upTime = upTime(2:end-1);
%extract breath and peak frequency waveforms for vocal times
spectVocStart = vocStarts*upSampleFactor;
spectVocEnd = vocEnds*upSampleFactor;
for vocNum = 1:length(spectVocStart)
    tempStart = round(spectVocStart(vocNum));
    tempEnd = round(spectVocEnd(vocNum));
    tempBreath = upSampBT(tempStart:tempEnd);
    tempFreq = threshFreq(tempStart:tempEnd);
    negBreathTrough = min(tempBreath) < 0;
    if negBreathTrough == 0
        normTempBreath = tempBreath-min(tempBreath);
        normTempBreath = normTempBreath/max(normTempBreath);
    else
        normTempBreath = tempBreath+(-min(tempBreath));
        normTempBreath = normTempBreath/max(normTempBreath);
    end
    normTempFreq = tempFreq-min(tempFreq);
    normTempFreq =normTempFreq/max(normTempFreq);
    normFreqs{vocNum,1} = normTempFreq;
    absFreqs{vocNum,1} = tempFreq;
    normUpExps{vocNum,1} = normTempBreath;
    absUpExps{vocNum,1} = tempBreath;
    r = corrcoef(normTempBreath,normTempFreq,'Rows','complete');
    corrs(vocNum) = r(2);
end
%% Identify multisyllabic vocalizations (same vocExpStart)
multiSyllInd = or([diff(vocExpStarts) 1] ==0, [1 diff(vocExpStarts)] ==0);
uniSyllInd = ~(multiSyllInd);
[uniqVoc,~,rnk] = unique(vocInspStarts);
rnk = transpose(rnk);
numVocBreaths = length(uniqVoc);
numMSCalls = sum(histcounts(vocInspStarts,unique(vocInspStarts))>1);
%number of syllables within a given vocalization breath
[vocCounts,~] = histc(vocInspStarts,uniqVoc);
biVocStarts = uniqVoc(vocCounts ==2);
triVocStarts = uniqVoc(vocCounts ==3);
[biInd,~] = ismember(vocInspStarts,biVocStarts);
[triInd,~] = ismember(vocInspStarts,triVocStarts);

%% Assemble unisyllabic data
cellNum = 1;
%vocClass = table2array(vocClass);
for uniNum = 1:length(uniSyllInd);
    if uniSyllInd(uniNum) == 1;
        uniCell{cellNum,1} = cell2mat(vocMat(1,uniNum));
        tempStart = (vocExpStarts(uniNum))-1000;
        tempEnd = (vocExpStarts(uniNum))+1000;
        uniCell{cellNum,2} = filtBreathTrace(tempStart:tempEnd);
        uniCell{cellNum,3} = vocInspStarts(uniNum);
        uniCell{cellNum,4} = vocExpStarts(uniNum);
        uniCell{cellNum,5} = vocStarts(uniNum);
        uniCell{cellNum,6} = vocEnds(uniNum);
        uniCell{cellNum,7} = vocClass(1,uniNum);
        uniCell{cellNum,8} = corrs(1,uniNum);
        uniCell{cellNum,9} = cell2mat(normFreqs(uniNum,1));
        uniCell{cellNum,10} = cell2mat(absFreqs(uniNum,1));
        uniCell{cellNum,11} = cell2mat(normUpExps(uniNum,1));
        uniCell{cellNum,12} = cell2mat(absUpExps(uniNum,1));
        cellNum = cellNum+1;
    end
end

%% Assemble bisyllabic data
cellNum = 1;
firstSyllNum = 1;
for biNum = 1:length(biInd);
    if biInd(biNum) == 1;
        if firstSyllNum == 1;
            biCell{cellNum,1} = cell2mat(vocMat(1,biNum));
            tempStart = (vocExpStarts(biNum))-1000;
            tempEnd = (vocExpStarts(biNum))+1000;
            biCell{cellNum,2} = filtBreathTrace(tempStart:tempEnd);
            biCell{cellNum,3} = vocInspStarts(biNum);
            biCell{cellNum,4} = vocExpStarts(biNum);
            biCell{cellNum,5} = vocStarts(biNum);
            biCell{cellNum,6} = vocEnds(biNum);
            biCell{cellNum,7} = vocClass(1,biNum);
            biCell{cellNum,8} = corrs(1,biNum);
            biCell{cellNum,9} = cell2mat(normFreqs(biNum,1));
            firstSyllNum = 0;
        elseif firstSyllNum == 0;
            biCell{cellNum,10} = vocStarts(biNum);
            biCell{cellNum,11} = vocEnds(biNum);
            biCell{cellNum,12} = vocClass(1,biNum);
            biCell{cellNum,13} = corrs(1,biNum);
            biCell{cellNum,14} = cell2mat(normFreqs(biNum,1));
            cellNum = cellNum+1;
            firstSyllNum = 1;
        end
    end
end

%% Assemble trisyllabic data
cellNum = 1;
firstSyllNum = 1;
secondSyllNum = 0;
thirdSyllNum = 0;
for triNum = 1:length(triInd);
    if triInd(triNum) == 1;
        if firstSyllNum == 1;
            triCell{cellNum,1} = cell2mat(vocMat(1,triNum));
            tempStart = (vocExpStarts(triNum))-1000;
            tempEnd = (vocExpStarts(triNum))+1000;
            triCell{cellNum,2} = filtBreathTrace(tempStart:tempEnd);
            triCell{cellNum,3} = vocInspStarts(triNum);
            triCell{cellNum,4} = vocExpStarts(triNum);
            triCell{cellNum,5} = vocStarts(triNum);
            triCell{cellNum,6} = vocEnds(triNum);
            triCell{cellNum,7} = vocClass(1,triNum);
            triCell{cellNum,8} = corrs(1,triNum);
            triCell{cellNum,9} = cell2mat(normFreqs(triNum,1));
            firstSyllNum = 0;
            secondSyllNum = 1;
        elseif secondSyllNum == 1;
            triCell{cellNum,10} = vocStarts(triNum);
            triCell{cellNum,11} = vocEnds(triNum);
            triCell{cellNum,12} = vocClass(1,triNum);
            triCell{cellNum,13} = corrs(1,triNum);
            triCell{cellNum,14} = cell2mat(normFreqs(triNum,1));
            secondSyllNum = 0;
            thirdSyllNum = 1;
        elseif thirdSyllNum == 1;
            triCell{cellNum,15} = vocStarts(triNum);
            triCell{cellNum,16} = vocEnds(triNum);
            triCell{cellNum,17} = vocClass(1,triNum);
            triCell{cellNum,18} = corrs(1,triNum);
            triCell{cellNum,19} = cell2mat(normFreqs(triNum,1));
            thirdSyllNum = 0;
            firstSyllNum = 1;
            cellNum = cellNum+1;
        end
    end
end
triVar = exist ('triCell');
biVar = exist('biCell');
uniVar = exist('uniCell');
numVars = sum([triVar biVar uniVar]);
% if numVars == 3
%     save (fullfile(pathToData,[fileName 'processed.mat']), 'uniCell','biCell','triCell');
% elseif numVars == 2
%     save (fullfile(pathToData,[fileName 'processed.mat']), 'uniCell','biCell');
% elseif numVars == 1
%     save (fullfile(pathToData,[fileName 'processed.mat']), 'uniCell');
% end

%% color coded scatters of breaths & freqs
vocIndNum = 1;
tempFreq = normFreqs{vocIndNum};
tempExp = normUpExps{vocIndNum};
time = (1:length(tempExp))/2;
c1 = [0.02 0.1 0.2];
c2 = [0.1 0.3 0.53];
c3 = [0.83 0.85 0.92];
c = generateCmap(c1,c2,c3);
steps = length(c)/length(tempExp);
colorInd = 1:steps:length(c);
colorInd = round(colorInd);
tempC = c(colorInd,:);
colormap(tempC);
figure
subplot (2,1,1)
scatter (time,tempExp,[],tempC,'_');
hold on
plot (time,tempFreq,'Color','r');
xlabel('time(ms)');
subplot (2,1,2)
scatter (tempExp,tempFreq,[],tempC,'filled');
hold on
%h = colorbar(tempC);
nans = isnan(tempFreq);
notNans = ~nans;
tempFreq = tempFreq(notNans);
tempExp = tempExp(notNans);
degree = 1;
coefficients = polyfit(tempExp,tempFreq,degree);
xFit = min(tempExp):0.1:max(tempExp);
yFit = polyval(coefficients,xFit);
plot (xFit,yFit);
%%
r = corrcoef (tempExp,tempFreq);