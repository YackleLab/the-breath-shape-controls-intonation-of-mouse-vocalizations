%% writeOptoVarsBatch
%Get Opto Stats outputs 5 x Nbreaths matrices describing Pi, Pe, Ti, Te &
%instFreq for opto, baseline,pre and post laser epochs
function [optoStats,blStats,preStats,postStats] = getOptoAmplitude(fileName,stimLength)
breathFile = strcat(fileName,'.txt');
breathData = readtable(breathFile);
breathTrace = table2array(breathData);
sampleRate = 1000;
time = 0:1/sampleRate:length(breathTrace)/sampleRate;
time = time(1:length(time)-1);
filteredBreathTrace = highpass(breathTrace,2,1000);
smoothFilteredBreathTrace = smoothdata (filteredBreathTrace, 'movmean',15);
stimLength = 1*1000;
% Import opto trace
optoFile = strcat (fileName,'_opto.txt');
optoData=readtable(optoFile, 'ReadVariableNames', false);
optoTrace=table2array(optoData);
%see if opto oversampled or not
if length(optoTrace)> length(breathTrace)
    downSampFactor = length(optoTrace)/length(breathTrace);
    optoTrace = downsample(optoTrace,downSampFactor);
end
%binarizeOptoTrace
optoTrace(optoTrace<5) = 0;
optoTrace(optoTrace>5) = 5;
[~, optoStarts] = findpeaks(optoTrace, 'MinPeakDistance', stimLength+1, 'MinPeakHeight', 2.5);
optoEnds = optoStarts + stimLength;
% Get Breath Parameters
durThresh = 0.00002; %minimum duration of insp
inspAmpThresh = -0.015; %minimum insp amplitude
expAmpThresh = 0; %minimum exp amplitude
breathStarts = getbreathstarts(smoothFilteredBreathTrace, durThresh, inspAmpThresh, expAmpThresh);
breathStarts = breathStarts(2:end);
breathStartsY = zeros(1,length(breathStarts));
breathStartsSec = breathStarts/sampleRate;
breathEnds = breathStarts (2:length(breathStarts));
%pull out breath waveforms
breathmat = cell (1,length(breathStarts)-1);
for breathindx = 1:(length(breathStarts)-1);
    tempbreathstart=breathStarts(breathindx);
    tempbreathend = breathStarts(breathindx+1);
    tempbreath = smoothFilteredBreathTrace(tempbreathstart:tempbreathend);
    breathmat{1,breathindx} = tempbreath;
end
[inspPeak, expPeak, inspDur, expDur, ~, ~, ~] = getbreathvals(breathmat);
inspPeak=inspPeak';expPeak=expPeak';inspDur=inspDur';expDur=expDur';
breathLength = (inspDur+expDur)/1000;
instFreq = 1./breathLength;
%find insp and exp starts
inspStarts = breathStarts (1:length(breathStarts)-1);
expStarts = inspStarts + inspDur;
optoIndStart = 1;
%find breaths during opto periods
for i = 1:length(optoStarts);
    optoInd(optoIndStart:optoIndStart+stimLength) = optoStarts(i,1):optoEnds(i,1);
    optoIndStart = length(optoInd)+1;
end
optoBreathInd = ismember(breathStarts,optoInd);
%get stats of optobreaths
optoPi = inspPeak(optoBreathInd);
optoPe = expPeak(optoBreathInd);
optoTi = inspDur(optoBreathInd);
optoTe = expDur(optoBreathInd);
optoFreq = instFreq(optoBreathInd);
optoStats(1,:) = optoPi;optoStats(2,:)=optoPe;optoStats(3,:)=optoTi;optoStats(4,:)=optoTe;optoStats(5,:)=optoFreq;
%get stats of bl breaths
blBreathInd = breathStarts < optoStarts(1,1);
blPi = inspPeak(blBreathInd);
blPe = expPeak(blBreathInd);
blTi = inspDur(blBreathInd);
blTe = expDur(blBreathInd);
blFreq = instFreq(blBreathInd);
blStats(1,:) = blPi;blStats(2,:)=blPe;blStats(3,:)=blTi;blStats(4,:)=blTe;blStats(5,:)=blFreq;
%get stats of pre breaths
preInd = optoInd-(1*sampleRate);
preBreathInd = ismember(breathStarts,preInd);
prePi = inspPeak(preBreathInd);
prePe = expPeak(preBreathInd);
preTi = inspDur(preBreathInd);
preTe = expDur(preBreathInd);
preFreq = instFreq(preBreathInd);
preStats(1,:) = prePi;preStats(2,:)=prePe;preStats(3,:)=preTi;preStats(4,:)=preTe;preStats(5,:)=preFreq;
%get stats of post breaths
postInd = optoInd+stimLength;
postBreathInd = ismember(breathStarts,postInd);
postPi = inspPeak(postBreathInd);
postPe = expPeak(postBreathInd);
postTi = inspDur(postBreathInd);
postTe = expDur(postBreathInd);
postFreq = instFreq(postBreathInd);
postStats(1,:) = postPi;postStats(2,:)=postPe;postStats(3,:)=postTi;postStats(4,:)=postTe;postStats(5,:)=postFreq;
end