function [saveDest,numcalls] = writeUSVvals(fileName,pathToData)
%writeUSVvals.m
%   This function outputs a mat file containing descriptive statistics of
%   breathing and vocalization from a single recording
%% Load in breathing data and filter
txtFile = append (fileName,'.txt');
cd(pathToData)
breathData = readtable (txtFile,'ReadVariableNames', false);
breathTrace = table2array (breathData);

breathTraceSampleRate = 1000; %Hz
time = 0:1/breathTraceSampleRate:length(breathTrace)/breathTraceSampleRate;
time = time(1:length(time)-1);

%Bandpass filter data
filtBreathTrace = bandpass (breathTrace,[2,35],1000);
%% Get Breath Parameters
% This section depends on functions from Bachmutsky et al 2020
% available on https://github.com/YackleLab/Opioids-depress-breathing-through-two-small-brainstem-sites
durThresh = 0.00002; %minimum duration of insp
inspAmpThresh = -0.015; %minimum insp amplitude
expAmpThresh = 0; %minimum exp amplitude
breathStarts = getbreathstarts(filtBreathTrace, durThresh, inspAmpThresh, expAmpThresh);
breathStarts = breathStarts(2:end);
breathStartsY = zeros(1,length(breathStarts));
breathStartsSec = breathStarts/breathTraceSampleRate;
breathEnds = breathStarts (2:length(breathStarts));
%pull out breath waveforms
breathmat = cell (1,length(breathStarts)-1);
for breathindx = 1:(length(breathStarts)-1);
    tempbreathstart=breathStarts(breathindx);
    tempbreathend = breathStarts(breathindx+1);
    tempbreath = filtBreathTrace(tempbreathstart:tempbreathend);
    breathmat{1,breathindx} = tempbreath;
end
[inspPeak, expPeak, inspDur, expDur, ~, ~, ~] = getbreathvals(breathmat);
%find insp and exp starts
inspStarts = breathStarts (1:length(breathStarts)-1);
inspDur = inspDur';
expStarts = inspStarts + inspDur;

%% USV detection with USVseg algorithm
%run USVdetection with modified USVseg algorithm
segFile = [fileName '_dat.csv'];
segQuery = isfile(segFile);
if segQuery == 0
    usvFile = USVdetect(fileName,pathToData);
    disp ('USVseg file not found. Running detection')
else
    usvFile = segFile;
    disp ('USVseg file found. Retrieving data')
end
usvData = readtable(usvFile,"VariableNamingRule","preserve");
vocStart = transpose(usvData.start);
vocEnd = transpose(usvData.end);

vocStartTime = vocStart*1000;
vocEndTime = vocEnd*1000;
%test if there are any vocalisations
noVoc = isempty(vocStartTime);

%% Aligning breath parameters and vocalization times
% discard the first and last 2s of breathing to remove filtering artefact
trimmedInspStarts = inspStarts;
trimmedFromStart = 0;
trimmedFromEnd = 0;
for i = 1:length(inspStarts)
    if inspStarts(i)<1500
        trimmedInspStarts = trimmedInspStarts (2:end);
        expStarts = expStarts (2:end);
        trimmedFromStart = trimmedFromStart+1;
    elseif inspStarts(i)>298500
        trimmedInspStarts = trimmedInspStarts (1:(length(trimmedInspStarts)-1));
        expStarts = expStarts (1:(length(trimmedInspStarts)-1));
        trimmedFromEnd = trimmedFromEnd+1;
    end
end
inspStarts = trimmedInspStarts;
inspPeak = inspPeak(trimmedFromStart+1:end-trimmedFromEnd);
expPeak = expPeak(trimmedFromStart+1:end-trimmedFromEnd);
%% recalculate some params ensure expStarts expEnds and inspStarts have the
%same number of elements
if inspStarts(1)>expStarts(1)
    inspStarts = inspStarts (2:end);
    inspPeak = inspPeak(2:end);
    expPeak = expPeak(2:end);
end
if expStarts(end)<inspStarts(end)
    inspStarts = inspStarts(1:(length(inspStarts)-1));
    inspPeak = inspPeak(2:end);
    expPeak = expPeak(2:end);
end
expEnds = inspStarts (2:end);
inspStarts = inspStarts(1:end-1);
inspPeak = inspPeak(1:end-1);
expPeak = expPeak(1:end-1);
if expStarts(end)>expEnds(end)
    expStarts = expStarts(1:end-1);
end
inspDur = expStarts-inspStarts;
expDur = expEnds-expStarts;
breathDur = expEnds - inspStarts; 
%% trim early and late vocalizations
if noVoc == 0
    trimmedVocStartTime = vocStartTime;
    for i = 1:length(vocStartTime)
        if vocStartTime(i)<1500
            trimmedVocStartTime = trimmedVocStartTime (2:end);
            vocEndTime = vocEndTime (2:end);
        elseif vocStartTime(i)>298500
            trimmedVocStartTime = trimmedVocStartTime (1:(length(trimmedVocStartTime)-1));
            vocEndTime = vocEndTime (1:(length(trimmedVocStartTime)));
        end
    end
    vocStartTime = trimmedVocStartTime;
    % check vocalization doesn't occur prior to first inspiration and discard
    % it if so
    if vocStartTime(1)<inspStarts(1)
        vocStartTime = vocStartTime (2:end);
        vocEndTime = vocEndTime (2:end);
    end
end
%% Calculate descriptive statistics of vocalization relative to breathing

%if there are none then set all parameters to zero
if noVoc == 1
    numcalls = 0;
    numMScalls = 0;
    numUniCalls = 0;
    numBiCalls = 0;
    numTriCalls = 0;
    interVocInterval = NaN;
    probConVocBreaths = 0;
    delayOnsetExp = NaN;
    normDelayOnsetExp = NaN;
    unisyllabicdelayOnsetExp = NaN;
    unisyllabicNormDelayOnsetExp = NaN;
    multiDelayOnsetExp = NaN;
    multiNormDelayOnsetExp = NaN;
    biDelayOnsetExp = NaN;
    biNormDelayOnsetExp = NaN;
    triDelayOnsetExp = NaN;
    triNormDelayOnsetExp = NaN;
    vocCum = zeros (1,10);
    vocHist = zeros (1,10);
    vocWinInspDur = NaN;
    vocWinExpDur = NaN;
    vocWinInspPeak = NaN;
    vocWinExpPeak = NaN;
elseif noVoc ==0
%identify breaths with vocalization
%pull out indices of inspStarts with vocalizations
vocBreathInspStartInd = interp1 (inspStarts,1:numel(inspStarts), vocStartTime,'previous','extrap');
%use these indices to find stats of vocalization breaths (i.e. value of
%that column in the vector)
vocExpStart = expStarts(vocBreathInspStartInd);
vocInspStart = inspStarts(vocBreathInspStartInd);
vocBreathDur = breathDur(vocBreathInspStartInd);
vocExpDur = expDur(vocBreathInspStartInd);
delayOnsetExp = vocStartTime - vocExpStart;
delayOffsetExp = vocEndTime-vocExpStart;
delayOnsetInsp = vocStartTime - vocInspStart;
delayOffsetInsp = vocEndTime-vocInspStart;
normDelayOnsetExp = delayOnsetExp./vocExpDur;
normDelayOffsetExp =  delayOffsetExp./vocExpDur;
normDelayOnsetInsp = delayOnsetInsp./vocBreathDur;
normDelayOffsetInsp = delayOffsetExp./vocBreathDur;
%% identify multisyllabic vocalizations (same vocExpStart)
multisyllabindex = or ([diff(vocExpStart) 1] == 0, [1 diff(vocExpStart)] == 0);
singlesyllabindex = ~multisyllabindex;
[uniqcalls, ~, rnk] = unique(vocInspStart);
rnk = transpose (rnk);
numcalls = length (uniqcalls); 
numMScalls = sum(histcounts(vocInspStart, unique(vocInspStart)) > 1);
%number of calls within a given vocalisation breath
[callscounts, ~] = histc(vocInspStart, uniqcalls);
bisyllabvocstarts = uniqcalls(callscounts ==2);
trisyllabvocstarts = uniqcalls(callscounts >=3);
[bisyllabindex,~] = ismember(vocInspStart,bisyllabvocstarts);
[trisyllabindex, ~] = ismember(vocInspStart,trisyllabvocstarts);
unisyllabicdelayOnsetExp = delayOnsetExp(singlesyllabindex);
unisyllabicNormDelayOnsetExp = normDelayOnsetExp(singlesyllabindex);
multiDelayOnsetExp = delayOnsetExp(multisyllabindex);
multiNormDelayOnsetExp = normDelayOnsetExp(multisyllabindex);
biDelayOnsetExp = delayOnsetExp(bisyllabindex);
biNormDelayOnsetExp = normDelayOnsetExp(bisyllabindex);
triDelayOnsetExp = delayOnsetExp(trisyllabindex);
triNormDelayOnsetExp = normDelayOnsetExp(trisyllabindex);
numUniCalls = length (unisyllabicNormDelayOnsetExp);
numBiCalls = length (bisyllabvocstarts);
numTriCalls = length (trisyllabvocstarts);
%calculate inter-Vocalization interval
interVocInterval = (vocStartTime(2:end))-(vocEndTime(1:end-1));
%find indicies for breaths without vocalization
inspStartInd = [1:length(inspStarts)];
noVocLog = ismember(inspStartInd, vocBreathInspStartInd);
noVocNum = 1;
for i = 1:length(noVocLog)
    if noVocLog (1,i) == 0
        noVocBreathInspStartInd(noVocNum)= inspStartInd(i);
        noVocNum = noVocNum +1;
    end
end
%calculate the probability that a voc breath will be followed by a second
%voc breath
consecutiveVocBreaths = 0;
for i = 1:length(vocBreathInspStartInd)
    testInd = vocBreathInspStartInd(i)+1;
    testlog = ismember(vocBreathInspStartInd,testInd);
    testlog = sum (testlog);
    if testlog == 1
        consecutiveVocBreaths = consecutiveVocBreaths+1;
    end
end
probConVocBreaths = (consecutiveVocBreaths/length(vocBreathInspStartInd))*100;
%calculate cumulative USV number
vocStartTimeSec = vocStartTime/1000;
timeBins = [0:30:300];
vocHist = histcounts (vocStartTimeSec,timeBins);
vocCum = cumsum (vocHist);
end
%calculate instantaneous frequency
breathDurSec = breathDur/1000;
for tempBreathNum = 1:length(breathDurSec)
    instFreq(tempBreathNum) = 1/breathDurSec(tempBreathNum);
end
if noVoc == 0
    vocBreathFreq = instFreq(vocBreathInspStartInd);
    noVocBreathFreq = instFreq (noVocBreathInspStartInd);
else
    vocBreathFreq = NaN;
    noVocBreathFreq = instFreq;
end
%pull out only in correct area (200s after urine presentation)
if contains (fileName,'urine1') == 1
    windowInds = breathStartsSec<200;
    windowFreqs = instFreq(windowInds);
    windowTi = inspDur(windowInds);
    windowTe = expDur(windowInds);
    windowPi = inspPeak(windowInds);
    windowPe = expPeak(windowInds);
    vocWindowInds = vocStartTimeSec<200;
    vocBreathWindowInds = vocBreathInspStartInd(vocWindowInds);
    vocWindowFreq = windowFreqs(vocBreathWindowInds);
    vocWinTi = windowTi(vocBreathWindowInds);
    vocWinTe = windowTe(vocBreathWindowInds);
    vocWinPi = windowPi(vocBreathWindowInds);
    vocWinPe = windowPe(vocBreathWindowInds);
    noVocBreathStartsSec = breathStartsSec(noVocBreathInspStartInd);
    noVocWindowInds = noVocBreathStartsSec<200;
    noVocBreathWindowInds = noVocBreathInspStartInd(noVocWindowInds);
    noVocWindowFreq = windowFreqs(noVocBreathWindowInds);
    noVocWinTi = windowTi(noVocBreathWindowInds);
    noVocWinTe = windowTe(noVocBreathWindowInds);
    noVocWinPi = windowPi(noVocBreathWindowInds);
    noVocWinPe = windowPe(noVocBreathWindowInds);
elseif contains (fileName,'urine_1') == 1
    windowInds = breathStartsSec<200;
    windowFreqs = instFreq(windowInds);
    windowTi = inspDur(windowInds);
    windowTe = expDur(windowInds);
    windowPi = inspPeak(windowInds);
    windowPe = expPeak(windowInds);
    vocWindowInds = vocStartTimeSec<200;
    vocBreathWindowInds = vocBreathInspStartInd(vocWindowInds);
    vocWindowFreq = windowFreqs(vocBreathWindowInds);
    vocWinTi = windowTi(vocBreathWindowInds);
    vocWinTe = windowTe(vocBreathWindowInds);
    vocWinPi = windowPi(vocBreathWindowInds);
    vocWinPe = windowPe(vocBreathWindowInds);
    noVocBreathStartsSec = breathStartsSec(noVocBreathInspStartInd);
    noVocWindowInds = noVocBreathStartsSec<200;
    noVocBreathWindowInds = noVocBreathInspStartInd(noVocWindowInds);
    noVocWindowFreq = windowFreqs(noVocBreathWindowInds);
    noVocWinTi = windowTi(noVocBreathWindowInds);
    noVocWinTe = windowTe(noVocBreathWindowInds);
    noVocWinPi = windowPi(noVocBreathWindowInds);
    noVocWinPe = windowPe(noVocBreathWindowInds);
else vocWindowFreq = NaN;
    noVocWindowFreq = NaN;
    vocWinTi = NaN;
    vocWinTe = NaN;
    vocWinPi = NaN;
    vocWinPe = NaN;
    noVocWinTi = NaN;
    noVocWinTe = NaN;
    noVocWinPi = NaN;
    noVocWinPe = NaN;
end
%create histogram of breathing rate
[breathRateHist, ~] = histcounts(breathStartsSec,0:5:300);
if noVoc == 0
    [vocRateHist,~] = histcounts(vocStartTimeSec,0:5:300);
else vocRateHist = zeros(1,60);
end
saveDest = [pathToData fileName '_callstats.mat'];
save (fullfile(saveDest),'numcalls','numMScalls','numUniCalls', ...
    'numBiCalls','numTriCalls', 'delayOnsetExp','normDelayOnsetExp', ...
    'unisyllabicdelayOnsetExp','unisyllabicNormDelayOnsetExp', ...
    'multiDelayOnsetExp','multiNormDelayOnsetExp','biDelayOnsetExp', ...
    'biNormDelayOnsetExp','triDelayOnsetExp','triNormDelayOnsetExp', ...
    'interVocInterval', 'probConVocBreaths','vocCum','vocHist', ...
    'instFreq','vocBreathFreq','noVocBreathFreq','vocWinTi',"vocWinTe",'vocWinPi', ...
    'vocWinPe','noVocWinTi','noVocWinTe','noVocWinPi','noVocWinPe', ...
    'breathRateHist','vocRateHist','vocWindowFreq','noVocWindowFreq');
end