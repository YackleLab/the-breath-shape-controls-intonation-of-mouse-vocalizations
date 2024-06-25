%% Load in the USV variables of mice & perform group analysis
%Analyze files with batchProcessStats.m
%Create metadata file in excel
%% Load metadata and extract filenames
clear all
close all
pathToData = 'C:\Users\Alastair\Box\USV Behavior\Urine recordings\analyze\';
metaDataFile = [pathToData,'metadata.xlsx'];
metaData = readtable (metaDataFile,'ReadVariableNames',true);
metaData = table2struct(metaData);
%% metaData = table2struct (metaData);
for i = 1:size(metaData,1)
    blFile {i,1} =  [(metaData(i).bl) '_callstats'];
    u1File {i,1} = strcat((metaData(i).u1),'_callstats');
    u2File {i,1} = strcat((metaData(i).u2),'_callstats');
    u3File {i,1} = strcat((metaData(i).u3),'_callstats');
end
%% Import call vars and store
for i = 1:length (blFile)
    fileName = strcat (pathToData, blFile{i,1},'.mat');
    load (fileName);
    BLnumCalls (i,:) = numcalls;
    BLnumMSCalls (i,:) = numMScalls;
    BLvocHist (i,:) = vocHist;
    blInstFreq {i,1} = instFreq;
    blNoVocFreq {i,1} = noVocBreathFreq;
    blBreathHist(i,:) = breathRateHist;
    blVocHist(i,:) = vocRateHist;
end
%%
for i = 1:length (u1File)
    fileName = strcat (pathToData, u1File{i,1},'.mat');
    load (fileName)
    u1numCalls(i,1) = numcalls;
    u1vocHist (i,:) = vocHist;
    u1instFreq {i,1} = instFreq;
    u1VocTi{i,:} = vocWinTi;
    u1VocTe{i,:} = vocWinTe;
    u1VocPi{i,:} = vocWinPi;
    u1VocPe{i,:} = vocWinPe;
    u1SniffTi{i,:} = noVocWinTi;
    u1SniffTe{i,:} = noVocWinTe;
    u1SniffPi{i,:} = noVocWinPi;
    u1SniffPe{i,:} = noVocWinPe;
    u1BreathHist(i,:) = breathRateHist;
    u1VocHist(i,:) = vocRateHist;
    u1VocWindowFreq{i,1} = vocWindowFreq;
    u1NoVocWindowFreq{i,1} = noVocWindowFreq;
end
for i = 1:length (u2File)
    fileName = strcat (pathToData, u2File{i,1},'.mat');
    load (fileName)
    u2numCalls(i,1) = numcalls;
    u2vocHist (i,:) = vocHist;
    u2instFreq {i,1} = instFreq;
    u2noVocFreq {i,1} = noVocBreathFreq;
    u2VocFreq {i,1} = vocBreathFreq;
    u2BreathHist(i,:) = breathRateHist;
    u2VocHist(i,:) = vocRateHist;
end
for i = 1:length (u3File)
    fileName = strcat (pathToData, u3File{i,1},'.mat');
    load (fileName)
    u3numCalls(i,1) = numcalls;
    u3vocHist (i,:) = vocHist;
    u3instFreq {i,1} = instFreq;
    u3noVocFreq {i,1} = noVocBreathFreq;
    u3VocFreq {i,1} = vocBreathFreq;
    u3BreathHist(i,:) = breathRateHist;
    u3VocHist(i,:) = vocRateHist;
end

%% rate histograms
breathRate = [blBreathHist u1BreathHist u2BreathHist u3BreathHist];
vocRate = [blVocHist u1VocHist u2VocHist u3VocHist];
%% Cumulative USV number
figure (3)
cumhist = horzcat (BLvocHist,u1vocHist,u2vocHist,u3vocHist);
cumUSV = cumsum (cumhist,2);
cumTime = [-4.5:0.5:15];
err = std(cumUSV);
err = err/sqrt(3);
meanCumUSV = mean (cumUSV,1);
H = shadedErrorBar (cumTime,meanCumUSV,err,'lineProps','b','patchSaturation',0.15);
H.mainLine.Color = 'r';
H.patch.FaceColor = 'r';
xline(0,':'); 
xlabel('Time (min)');
ylabel('Cumulative # USVs');
ylim ([0 500]);

%% Vocal breath params
allVocTi = [0]; allVocTe = [0]; allVocPi = [0]; allVocPe = [0];
allSniffTi = [0]; allSniffTe = [0]; allSniffPi = [0]; allSniffPe = [0];
for i = 1:length(u1VocTi);
    tempVocTi = u1VocTi{i,1};allVocTi = [allVocTi tempVocTi];meanVocTi(i) = mean(tempVocTi);
    tempVocTe = u1VocTe{i,1};allVocTe = [allVocTe tempVocTe];meanVocTe(i) = mean(tempVocTe);
    tempVocRatio = tempVocTi/tempVocTe; meanVocRatio(i) = mean(tempVocRatio);
    tempVocPi = transpose(u1VocPi{i,1});allVocPi = [allVocPi tempVocPi];meanVocPi(i) = mean(tempVocPi);
    tempVocPe = transpose(u1VocPe{i,1});allVocPe = [allVocPe tempVocPe];meanVocPe(i) = mean(tempVocPe);
    tempSniffTi = u1SniffTi{i,1};allSniffTi = [allSniffTi tempSniffTi];meanSniffTi(i) = mean(tempSniffTi);
    tempSniffTe = u1SniffTe{i,1};allSniffTe = [allSniffTe tempSniffTe];meanSniffTe(i) = mean(tempSniffTe);
    tempSniffRatio = tempSniffTi/tempSniffTe; meanSniffRatio(i) = mean(tempSniffRatio);
    tempSniffPi = transpose(u1SniffPi{i,1});allSniffPi = [allSniffPi tempSniffPi];meanSniffPi(i) = mean(tempSniffPi);
    tempSniffPe = transpose(u1SniffPe{i,1});allSniffPe = [allSniffPe tempSniffPe];meanSniffPe(i) = mean(tempSniffPe);
end
allVocTi = allVocTi(2:end);,allVocTe = allVocTe(2:end);,allVocPi = allVocPi(2:end);,allVocPe = allVocPe(2:end);
allSniffTi = allSniffTi(2:end);,allSniffTe = allSniffTe(2:end);,allSniffPi = allSniffPi(2:end);,allSniffPe = allSniffPe(2:end);
%ti
fig1chart2vars(vertcat(meanVocTi,meanSniffTi));
ylabel ('Ti (ms)');
title ('Inspiration Duration');
ylim([0 50]);
%te
fig1chart2vars(vertcat(meanVocTe,meanSniffTe));
ylabel ('Te (ms)');
title ('Expiration Duration');
ylim([0 140]);
%pi
fig1chart2vars(vertcat(meanVocPi,meanSniffPi));
ylabel ('Pi (au)');
title ('Inspiration Amplitude');
ylim([-0.06 0]);
%pe
fig1chart2vars(vertcat(meanVocPe,meanSniffPe));
ylabel ('Pe (au)');
title ('Expiration Amplitude');
ylim([0 0.04]);
%ti/te
fig1chart2vars(vertcat(meanVocRatio,meanSniffRatio));
ylabel('Ti/Te');
title('insp:exp ratio');
%% rate histograms
figure
breathRate = breathRate/5;
meanBreathRate = mean(breathRate);
semBreathRate = (std(breathRate))/sqrt(5);
h = shadedErrorBar(1:5:1200, meanBreathRate,semBreathRate,'lineProps','b','patchSaturation',0.15);
h.mainLine.Color = 'k';
h.patch.FaceColor = 'k';
hold on
vocRate = vocRate/5;
meanVocRate = mean(vocRate);
semVocRate = (std(vocRate))/sqrt(5);
j = shadedErrorBar(1:5:1200, meanVocRate, semVocRate,'lineProps','b','patchSaturation',0.15);
j.mainLine.Color = 'r';
j.patch.FaceColor = 'r';
%% freq histograms
edges = 0:0.5:20;
for i = 1:6
tempBl = blNoVocFreq{i,1};
tempBlHist = histcounts(tempBl,edges);
tempBlHist = tempBlHist/length(tempBl);
tempVoc = u1VocWindowFreq{i,1};
tempVocHist = histcounts(tempVoc,edges);
tempVocHist = tempVocHist/length(tempVoc);
tempNoVoc = u1NoVocWindowFreq{i,1};
tempNoVocHist = histcounts(tempNoVoc,edges);
tempNoVocHist = tempNoVocHist/length(tempNoVoc);
blWindowHist(i,:) = tempBlHist;
vocalWindowHist(i,:) = tempVocHist;
sniffWindowHist(i,:) = tempNoVocHist;
end

%% plot for individual mice
figure
%plot (edges(2:end),blWindowHist, 'Color', [0.5 0.5 0.5],'LineWidth',0.25);
hold on
plot (edges(2:end),vocalWindowHist,'Color','r','LineWidth',0.25);
plot (edges(2:end), sniffWindowHist,'Color','k','LineWidth',0.25);
%plot (edges(2:end),mean(blWindowHist), 'Color', [0.5 0.5 0.5],'LineWidth',3);
plot (edges(2:end),mean(vocalWindowHist),'Color','r','LineWidth',3);
plot (edges(2:end), mean(sniffWindowHist),'Color','k','LineWidth',3);
%%
figure
k = shadedErrorBar(edges(2:end),mean(sniffWindowHist),(std(sniffWindowHist))/sqrt(5),'lineProps','b','patchSaturation',0.15);
k.mainLine.Color = 'k';
k.patch.FaceColor = 'k';
h = shadedErrorBar(edges(2:end),mean(vocalWindowHist),(std(vocalWindowHist))/sqrt(5),'lineProps','b','patchSaturation',0.15);
h.mainLine.Color = 'r';
h.patch.FaceColor = 'r';
%scatters
for i = 6
    vocTi = u1VocTi{i,:};
    vocTe = u1VocTe{i,:};
    vocPi = u1VocPi{i,:};
    vocPe = u1VocPe{i,:};
    sniffTi = u1SniffTi{i,:};
    sniffTe = u1SniffTe{i,:};
    sniffPi = u1SniffPi{i,:};
    sniffPe = u1SniffPe{i,:};
    %randomly pull out sniff breaths so Ns are equal
    randInd = rand(1,length(vocTi));
    randInd = randInd*length(sniffTi);
    randInd = round(randInd);
    sniffTi = sniffTi(randInd);
    sniffTe = sniffTe(randInd);
    sniffPi = sniffPi(randInd);
    sniffPe = sniffPe(randInd);
    figure
    scatter(sniffTi,sniffTe);
    hold on
    scatter(vocTi,vocTe);
    xlim([0 140])
    ylim([0 350])
    figure
    scatter(sniffPe,sniffPi);
    hold on
    scatter(vocPe,vocPi);
    xlim([0 0.2])
    ylim([-0.2 0])
end

%% get inst freq means
for i = 1:length(u1VocWindowFreq)
meanVocWindowFreq(i) = mean(u1VocWindowFreq{i,:});
end
for i = 1:length(u1NoVocWindowFreq)
meanNoVocWindowFreq(i) = mean(u1NoVocWindowFreq{i,:});
end
fig1chart2vars(vertcat(meanVocWindowFreq,meanNoVocWindowFreq));
ylim([0 15])
ylabel('inst freq (Hz)');
title('Instantaneous Frequency');