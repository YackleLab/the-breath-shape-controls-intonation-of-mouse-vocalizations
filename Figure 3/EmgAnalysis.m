%% EmgAnalysis.m
% For plotting and analysis of audio, airflow and EMG signals
% first compute frequency vector with peakFreqWriter.m
%% load in breathing data and filtering
pathToData = '';
fileName = '';
spectFile = '';
freqFile = '';
diaFile = '';
txtfile = strcat (fileName,'.txt');
emgfile = strcat (fileName,'_emg.txt');
emg2file = strcat (fileName,'_emg2.txt');
cd(pathToData)
breathdata=readtable(txtfile, 'ReadVariableNames', false);
emgData = readtable(emgfile, 'ReadVariableNames', false);
emg2Data = readtable(emg2file, 'ReadVariableNames', false);
emgTrace = table2array(emgData);
emg2Trace = table2array(emg2Data);
breathTrace=table2array(breathdata);
breathSampRate = 1000.; %hz
emgSamprate = 10000; %hz
time = 0:1/breathSampRate:length(breathTrace)/breathSampRate;
time = time(1:length(time)-1);
emgTime = 0:1/emgSamprate:length(emgTrace)/emgSamprate;
emgTime = emgTime(1:length(emgTime)-1);
filtBreathTrace=bandpass(breathTrace,[2,35],1000);% 
%% Filtering EMG
[B,A] = paynter (.02,1000,'mod'); %time constant 20 ms
intLarynx = filtfilt(B,A, abs(emgTrace));
intDia = filtfilt(B,A, abs(emg2Trace));
%% Get Breath Parameters
% This section depends on functions from Bachmutsky et al 2020
% available on https://github.com/YackleLab/Opioids-depress-breathing-through-two-small-brainstem-sites
durThresh = 0.00002; %minimum duration of insp
inspAmpThresh = -0.01; %minimum insp amplitude
expAmpThresh = 0; %minimum exp amplitude
breathStarts = getbreathstarts(filtBreathTrace, durThresh, inspAmpThresh, expAmpThresh);
breathStarts = breathStarts(2:end);
breathStartsY = zeros(1,length(breathStarts));
breathStartsSec = breathStarts/breathSampRate;
breathEnds = breathStarts (2:length(breathStarts));
%pull out breath waveforms
breathmat = cell (1,length(breathStarts)-1);
for breathindx = 1:(length(breathStarts)-1)
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
%% recalculate some params ensure expStarts and inspStarts have the
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
%% USV detection with USVseg algorithm
%run USVdetection with modified USVseg algorithm
% segFile = [fileName '_dat.csv'];
% segQuery = isfile(segFile);
% if segQuery == 0
%     usvFile = USVdetect(fileName,pathToData);
%     disp ('USVseg file not found. Running detection')
% else
%     usvFile = segFile;
%     disp ('USVseg file found. Retrieving data')
% end
% usvData = readtable(usvFile,"VariableNamingRule","preserve");
% vocStart = transpose(usvData.start);
% vocEnd = transpose(usvData.end);
%% find VocalMat analyzed file and import times
vmFile = strcat(fileName,'.xlsx');
vmQuery = isfile(vmFile);
if vmQuery == 1
    vmData = readtable(vmFile,'VariableNamingRule','preserve');
    vocStart = transpose(vmData.Start_time);
    vocEnd = transpose(vmData.End_time);
    vocClass = transpose(vmData.Class);
    noVm = isempty(vocStart);
else noVm = 1;
    disp ('No VocalMat file found')
end
%% trim USVs
%remove USVs identified as noise
for i = 1:length(vocClass)
    if contains(vocClass{1,i},'noise_dist')
        noiseInds(i) = true;
    else noiseInds(i) = false;
    end
end
vocStart = vocStart(~noiseInds);
vocEnd = vocEnd(~noiseInds);
vocClass = vocClass(~noiseInds);
%get rid of usv that occur before first breath
vocStartTime = vocStart*1000;
vocEndTime = vocEnd*1000;
vocEndTime = vocEndTime(vocStartTime>inspStarts(1));
vocStartTime = vocStartTime(vocStartTime>inspStarts(1));
%test if there are any vocalisations
noVoc = isempty(vocStartTime);
%% Calculate descriptive statistics of vocalization relative to breathing
%identify breaths with vocalization
%pull out indices of inspStarts with vocalizations
vocBreathInspStartInd = interp1 (inspStarts,1:numel(inspStarts), vocStartTime,'previous','extrap');
%use these indices to find stats of vocalization breaths (i.e. value of
%that column in the vector)
vocExpStart = expStarts(vocBreathInspStartInd);
vocInspStart = inspStarts(vocBreathInspStartInd);
vocBreathDur = breathDur(vocBreathInspStartInd);
vocExpDur = expDur(vocBreathInspStartInd);
%ignore double
multisyllabindex = or ([diff(vocExpStart) 1] == 0, [1 diff(vocExpStart)] == 0);
singlesyllabindex = ~multisyllabindex;
vocInspStart = vocInspStart(singlesyllabindex);
vocExpStart = vocExpStart(singlesyllabindex);
vocBreathDur = vocBreathDur(singlesyllabindex);
vocExpDur = vocExpDur(singlesyllabindex);
uniVocStart = vocStartTime(singlesyllabindex);
uniVocEnd = vocEndTime(singlesyllabindex);
uniVocClass = vocClass(singlesyllabindex);
%% Load in spectrogram and peak freq
load(spectFile);
load(freqFile);

%% pull out breathing, emg, and usv vector for each voc
close all
for i = 1:length(uniVocClass)
    if contains(uniVocClass{1,i},'complex')
        complexInds(i) = true;
    else complexInds(i) = false;
    end
end
complexStart = uniVocStart(complexInds);
complexEnd = uniVocEnd(complexInds);
complexLength = complexEnd - complexStart;
validComplex = complexLength> 40;
complexStart = complexStart(validComplex);
complexEnd = complexEnd(validComplex);
dsLarynx = downsample(intLarynx,10);
dsDia = downsample(intDia,10);
for i =1:length(uniVocStart)
    tempStart = round(complexStart(i));
    tempEnd = round(complexEnd(i));
    tempAir = filtBreathTrace(tempStart:tempEnd);
    tempLarynx = dsLarynx(tempStart:tempEnd);
    tempDia = dsDia(tempStart:tempEnd);
    tempFreq = threshFreq((tempStart*2)-1:(tempEnd*2)-1);
end
figure
%plot (1:length(tempAir),tempAir);
hold on
plot (1:length(tempAir),normDia);
plot (1:length(tempAir),normAir);
plot (1:length(tempAir),normLarynx);
plot((1:length(tempFreq))/2,normFreq);
legend('dia','airflow','larynx','freq');
%% plot spectrogram
figure;
ax1=subplot (4,1,1);
%imagesc(T, F, 10*log10(P+eps)); % add eps like pspectrogram does
imagesc(T,F,P);
axis xy;
colormap(flipud(gray));
ylabel('Frequency (kHz)');
xlabel('Time (s)');
caxis([5 10]);
ylim ([0 2e+5])
% h = colorbar;
% h.Label.String = 'Power/frequency (dB/Hz)';
%c = gray;
%c = flipud(c);
%colormap(c);
%clim([-102 -80]);
hold on
plot (T, threshFreq,'color','m');

%filtered_breathtrace = breathtrace;
ax2=subplot(4,1,2);
yline(0,':'); 
hold on;
plot (time, filtBreathTrace, 'k');
plot (time, highBreaths,'Color', [0.5 0.5 0.5]);
plot([vocStartTime/1000;vocEndTime/1000], [0,0],'LineWidth',5, 'Color', [0.4660 0.6740 0.1880]);
ylabel('Airflow') ;
xlabel('Time (sec)');


ax3=subplot(4,1,3); 
plot (emgTime, emgTrace, 'Color', [0.5 0.5 0.5]);
hold on;
plot (emgTime, intLarynx, 'k');
ylabel('int Larynx')

ax4=subplot(4,1,4);
plot (emgTime, emg2Trace,'Color', [0.5 0.5 0.5]);
hold on;
plot (emgTime, intDia, 'k');
ylabel('int Dia')
linkaxes ([ax1, ax2, ax3, ax4], 'x')