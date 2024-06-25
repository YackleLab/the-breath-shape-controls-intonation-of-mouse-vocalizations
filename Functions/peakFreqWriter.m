function threshFreq = peakFreqWriter(fileName,pathToData)
%peakFreqWriter - function version of write peak freq
%creates spectrogram, filters & thresholds then saves peak frequency
%as a vector
%find the spectrogram and either load or create
[spectData,fs] = audioread([fileName '.wav']);
spectFile = [fileName '_spect.mat'];
spectQuery = isfile(spectFile);
if spectQuery ==1
    load (spectFile);
else
    [~,F,T,P] = USVdetect (fileName,pathToData);
end
%cut out frequencies below ~40 kHz and ~120 kHZ
freqLut = F(52:155); %if endogenous
%freqLut = F (83:247); %if opto
powerSpec = P(52:155); %if endogenous
%powerSpec = P(83:247,:); %if opto
% cycle through time bins, pull out max power freq bin and store in a vector
for i = 1:length(T);
    [tempVal tempInd] = max(powerSpec(:,i));
    freqVals(i) = tempVal;
    freqInds(i) = tempInd;
    peakFreq(i) = freqLut(tempInd);
end
% threshold for values above 3SD
thresh = 3*std(powerSpec);
threshLog = gt(freqVals,thresh);
for i = 1:length(peakFreq)
    if threshLog (i) == 1
        threshFreq(i) = peakFreq(i);
    else threshFreq(i) = NaN;
    end
end
save (fullfile(pathToData,[fileName 'peakFreq.mat']), 'threshFreq');
end