%% Compile all calls type
% Compile call info from multiple recordings (analyzed by callVarsPreProcessing.m)
% into one dataset and calculate statistics sorted by call type
clear all
close all
%% Load metadata and extract filenames
pathToData = ''; %set filepath
cd(pathToData)
metaDataFile = [pathToData,'vocdata.xlsx'];
metaData = readtable (metaDataFile,'ReadVariableNames',true);
metaData = table2struct(metaData);

%% Import call vars and store
unis = {'breath trace','long trace','inspStart','expStart','vocStart','vocEnd','vocClass','corr','norm freq trace','abs freq trace','norm upsampled exp', 'abs upsampled exp'};
bis = {'breath trace','long trace','inspStart','expStart','vocStart1','vocEnd1','vocClass1','corr1','freq trace 1','vocStart2','vocEnd2','vocClass2','corr2','freq trace 2'};
tris = {'breath trace','long trace','inspStart','expStart','vocStart1','vocEnd1','vocClass1','corr1','freq trace 1','vocStart2','vocEnd2','vocClass2','corr2','freq trace 2','vocStart3','vocEnd3','vocClass3','corr3','freq trace 3'};
for i = 1:length(metaData)
    if isfile ([metaData(i).file 'processed.mat'])
        tempFile = [metaData(i).file 'processed'];
        load (tempFile)
        unis = [unis; uniCell];
        if exist('biCell') == 1
            bis = [bis; biCell];
        end
        if exist('triCell') == 1
            tris = [tris;triCell];
        end
    end
end
%clear uniCell; clear biCell; clear triCell;
unis = unis(2:end,:); bis = bis(2:end,:); tris = tris(2:end,:);
%get uni stats
expStarts = cell2mat(unis(:,4));
inspStarts = cell2mat(unis(:,3));
inspDur = expStarts-inspStarts;
breaths = unis(:,1);
vocStarts = cell2mat(unis(:,5));
vocEnds = cell2mat(unis(:,6));
vocOnsets = vocStarts-expStarts;
vocOffsets = vocEnds-expStarts;
fullBreaths =transpose(unis(:,2));
fullBreaths = cell2mat(fullBreaths);
allAbsFreq = [0];
allAbsExp = [0];
for i = 1:length(unis)
    tempBreath = breaths{i,1};
    breathLength(i,1) = length(tempBreath);
    tempFull = fullBreaths(:,i);
    vocStartBreathInd = round(vocOnsets(i,1)+1001);
    vocEndBreathInd = round(vocOffsets(i,1)+1001);
    vocExps{i,1} = tempFull(vocStartBreathInd:vocEndBreathInd);
    allAbsFreq = [allAbsFreq unis{i,10}];
    allAbsExp = [allAbsExp; unis{i,12}];
end
expDur = breathLength-inspDur;
corrs = cell2mat(unis(:,8));
%% Plot calltypes unisyllabic
flatNum=1;ufmNum=1;chevNum=1;dfmNum=1;revChevNum=1;
shortNum=1;compNum=1;stepUpNum=1;stepDownNum=1;twoStepNum=1;multiNum = 1;
for i = 1:length(unis)
    tempVal = unis{i,7};
    if contains(tempVal,'flat') == 1
        flatInd(flatNum) = i;
        flatNum = flatNum +1;
    elseif contains(tempVal,'up_fm') == 1
        ufmInd(ufmNum) = i;
        ufmNum = ufmNum +1;
    elseif contains(tempVal,'chevron') == 1
        chevInd(chevNum) = i;
        chevNum = chevNum +1;
    elseif contains(tempVal,'down_fm') == 1
        dfmInd(dfmNum) = i;
        dfmNum = dfmNum +1;
    elseif contains(tempVal,'rev_chevron') == 1
        revChevInd(revChevNum) = i;
        revChevNum = revChevNum +1;
    elseif contains(tempVal,'short') == 1
        shortInd(shortNum) = i;
        shortNum = shortNum +1;
    elseif contains(tempVal,'complex') == 1
        compInd(compNum) = i;
        compNum = compNum +1;
    elseif contains(tempVal,'step_up') == 1
        stepUpInd(stepUpNum) = i;
        stepUpNum = stepUpNum +1;
    elseif contains(tempVal,'step_down') == 1
        stepDownInd(stepDownNum) = i;
        stepDownNum = stepDownNum +1;
    elseif contains(tempVal,'two_steps') == 1
        twoStepInd(twoStepNum) = i;
        twoStepNum = twoStepNum +1;
    elseif contains(tempVal,'mult_steps') == 1
        multiInd(multiNum) = i;
        multiNum = multiNum +1;
    end
end
flatNum = flatNum-1;ufmNum = ufmNum-1;chevNum =chevNum-1;dfmNum =dfmNum-1;
revChevNum = revChevNum-1;shortNum = shortNum-1;compNum = compNum-1;stepUpNum = stepUpNum-1;
stepDownNum = stepDownNum-1;twoStepNum =twoStepNum-1;multiNum =multiNum-1;
allNums = [flatNum,ufmNum,chevNum,dfmNum,revChevNum,shortNum,compNum,stepUpNum,stepDownNum,twoStepNum,multiNum];
sortNums = sort(allNums,'descend');
propNums = allNums/sum(allNums);
figure (1)
labels = categorical({'up fm','step down','flat','short','complex','chevron','two step','multi step','step up','down fm','rev chev'});
labels = reordercats(labels,{'up fm','step down','flat','short','complex','chevron','two step','multi step','step up','down fm','rev chev'});
bar (labels,sortNums);
ufmInspDur=inspDur(ufmInd);ufmExpDur=expDur(ufmInd);ufmBreaths=breaths(ufmInd); ufmVocOnsets=vocOnsets(ufmInd); ufmVocOffsets=vocOffsets(ufmInd); ufmCorr = corrs(ufmInd);
chevInspDur=inspDur(chevInd);chevExpDur=expDur(chevInd);chevBreaths=breaths(chevInd);chevVocOnsets=vocOnsets(chevInd);chevVocOffsets=vocOffsets(chevInd); chevCorr = corrs(chevInd);
twoStepInspDur=inspDur(twoStepInd);twoStepExpDur=expDur(twoStepInd);twoStepBreaths=breaths(twoStepInd);twoStepVocOnsets=vocOnsets(twoStepInd);twoStepVocOffsets=vocOffsets(twoStepInd); twoStepCorr = corrs(twoStepInd);
stepDownInspDur=inspDur(stepDownInd);stepDownExpDur=expDur(stepDownInd);stepDownBreaths=breaths(stepDownInd);stepDownVocOnsets=vocOnsets(stepDownInd);stepDownVocOffsets=vocOffsets(stepDownInd); stepDownCorr = corrs(stepDownInd);
flatInspDur=inspDur(flatInd);flatExpDur=expDur(flatInd);flatBreaths=breaths(flatInd);flatVocOnsets=vocOnsets(flatInd);flatVocOffsets=vocOffsets(flatInd); flatCorr = corrs(flatInd);
dfmInspDur=inspDur(dfmInd);dfmExpDur=expDur(dfmInd);dfmBreaths=breaths(dfmInd);dfmVocOnsets=vocOnsets(dfmInd);dfmVocOffsets=vocOffsets(dfmInd); dfmCorr = corrs(dfmInd);
multiInspDur=inspDur(multiInd);multiExpDur=expDur(multiInd);multiBreaths=breaths(multiInd);multiVocOnsets=vocOnsets(multiInd);multiVocOffsets=vocOffsets(multiInd); multiCorr = corrs(multiInd);
stepUpInspDur = inspDur(stepUpInd);stepUpExpDur = expDur(stepUpInd);stepUpBreaths=breaths(stepUpInd);stepUpVocOnsets=vocOnsets(stepUpInd);stepUpVocOffsets=vocOffsets(stepUpInd); stepUpCorr = corrs(stepUpInd);
compInspDur=inspDur(compInd);compExpDur=expDur(compInd);compBreaths=breaths(compInd);compVocOnsets=vocOnsets(compInd);compVocOffsets=vocOffsets(compInd); compCorr = corrs (compInd);
shortInspDur=inspDur(shortInd);shortExpDur=expDur(shortInd);shortBreaths=breaths(shortInd);shortVocOnsets=vocOnsets(shortInd);shortVocOffsets=vocOffsets(shortInd); shortCorr = corrs(shortInd);
edges = [-100:10:300];
for i = 1:length(edges)-1
    centers(i) = (edges(i)+(edges(i+1)))/2;
end
%onset hists
[ufmHist ~] = histcounts(ufmVocOnsets,edges);ufmHist = ufmHist/ufmNum;
[chevHist ~] = histcounts(chevVocOnsets,edges);chevHist = chevHist/chevNum;
[twoStepHist ~] = histcounts(twoStepVocOnsets,edges);twoStepHist = twoStepHist/twoStepNum;
[stepDownHist ~] = histcounts(stepDownVocOnsets,edges);stepDownHist = stepDownHist/stepDownNum;
[flatHist ~] = histcounts(flatVocOnsets,edges);flatHist = flatHist/flatNum;
[dfmHist ~] = histcounts(dfmVocOnsets,edges);dfmHist = dfmHist/dfmNum;
[multiHist ~] = histcounts(multiVocOnsets,edges);multiHist = multiHist/multiNum;
[stepUpHist ~] = histcounts(stepUpVocOnsets,edges);stepUpHist = stepUpHist/stepUpNum;
[compHist ~] = histcounts(compVocOnsets,edges);compHist = compHist/compNum;
[shortHist ~] = histcounts(shortVocOnsets,edges);shortHist = shortHist/shortNum;
%offset hists
[ufmOffHist ~] = histcounts(ufmVocOffsets,edges);ufmOffHist = ufmOffHist/ufmNum;
[chevOffHist ~] = histcounts(chevVocOffsets,edges);chevOffHist = chevOffHist/chevNum;
[twoStepOffHist ~] = histcounts(twoStepVocOffsets,edges);twoStepOffHist = twoStepOffHist/twoStepNum;
[stepDownOffHist ~] = histcounts(stepDownVocOffsets,edges);stepDownOffHist = stepDownOffHist/stepDownNum;
[flatOffHist ~] = histcounts(flatVocOffsets,edges);flatOffHist = flatOffHist/flatNum;
[dfmOffHist ~] = histcounts(dfmVocOffsets,edges);dfmOffHist = dfmOffHist/dfmNum;
[multiOffHist ~] = histcounts(multiVocOffsets,edges);multiOffHist = multiOffHist/multiNum;
[stepUpOffHist ~] = histcounts(stepUpVocOffsets,edges);stepUpOffHist = stepUpOffHist/stepUpNum;
[compOffHist ~] = histcounts(compVocOffsets,edges);compOffHist = compOffHist/compNum;
[shortOffHist ~] = histcounts(shortVocOffsets,edges);shortOffHist = shortOffHist/shortNum;
%% Pull out waveforms
figure (2)
ax1 = subplot (2,5,1);
c = turbo;
colorInd = transpose(1:25:255);
hold on
yline(0,':');
xline(0,':'); 
tempCol = [0.5 0.5 0.5];
for i=1:ufmNum
    tempBreath = ufmBreaths{i,1};
    tempX = 1:length(tempBreath);
    tempX = tempX-ufmInspDur(i);
    tempCol = c(colorInd(2),:);
    plot (tempX,tempBreath,'Color',tempCol);
    ylim ([-0.14 0.08])
    xlim ([-100 300])
    ylabel('Airflow (au)');
    title ('up fm (467)');
    hold on
end
yyaxis right
plot (centers,ufmHist,'Color',tempCol,'LineWidth',3);
hold on
plot (centers,ufmOffHist,'Color',tempCol,'LineWidth',3,'LineStyle',':');
ylim ([0 1])
yticklabels([]);
ax2 = subplot (2,5,2);
yline(0,':'); 
xline(0,':'); 
hold on
for i=1:stepDownNum
    tempBreath = stepDownBreaths{i,1};
    tempX = 1:length(tempBreath);
    tempX = tempX-stepDownInspDur(i);
    tempCol = c(colorInd(3),:);
    plot (tempX,tempBreath,'Color',tempCol);
    ylim ([-0.14 0.08])
    xlim ([-100 300])
    title ('step down (355)');
end
yticklabels([]);
yyaxis right
plot (centers,stepDownHist,'Color',tempCol,'LineWidth',3);
hold on
plot (centers,stepDownOffHist,'Color',tempCol,'LineWidth',3,'LineStyle',':');
ylim ([0 1])
yticklabels([]);
ax3 = subplot (2,5,3);
yline(0,':'); 
xline(0,':'); 
hold on
for i=1:flatNum
    tempBreath = flatBreaths{i,1};
    tempX = 1:length(tempBreath);
    tempX = tempX-flatInspDur(i);
    tempCol = c(colorInd(4),:);
    plot (tempX,tempBreath,'Color',tempCol);
    ylim ([-0.14 0.08])
    xlim ([-100 300])
    title ('flat (272)');
end
yticklabels([]);
yyaxis right
plot (centers,flatHist,'Color',tempCol,'LineWidth',3);
hold on
plot (centers,flatOffHist,'Color',tempCol,'LineWidth',3,'LineStyle',':');
ylim ([0 1])
yticklabels([]);
ax4 = subplot (2,5,4);
yline(0,':'); 
xline(0,':'); 
yticklabels([]);
hold on
for i=1:shortNum
    tempBreath = shortBreaths{i,1};
    tempX = 1:length(tempBreath);
    tempX = tempX-shortInspDur(i);
    tempCol = c(colorInd(5),:);
    plot (tempX,tempBreath,'Color',tempCol);
    ylim ([-0.14 0.08])
    xlim ([-100 300])
    title ('short (158)');
end
yticklabels([]);
yyaxis right
plot (centers,shortHist,'Color',tempCol,'LineWidth',3);
hold on
plot (centers,shortOffHist,'Color',tempCol,'LineWidth',3,'LineStyle',':');
ylim ([0 1])
yticklabels([]);
ax5 = subplot (2,5,5);
yline(0,':'); 
xline(0,':'); 
hold on
for i=1:compNum
    tempBreath = compBreaths{i,1};
    tempX = 1:length(tempBreath);
    tempX = tempX-compInspDur(i);
    tempCol = c(colorInd(6),:);
    plot (tempX,tempBreath,'Color',tempCol);
    ylim ([-0.14 0.08])
    xlim ([-100 300])
    title ('complex (126)');
end
yticklabels([]);
yyaxis right
plot (centers,compHist,'Color',tempCol,'LineWidth',3);
hold on
plot (centers,compOffHist,'Color',tempCol,'LineWidth',3,'LineStyle',':');
ylim ([0 1])
ylabel('proportion of USVs');
ax6 = subplot (2,5,6);
yline(0,':'); 
xline(0,':'); 
hold on
for i=1:chevNum
    tempBreath = chevBreaths{i,1};
    tempX = 1:length(tempBreath);
    tempX = tempX-chevInspDur(i);
    tempCol = c(colorInd(7),:);
    plot (tempX,tempBreath,'Color',tempCol);
    ylim ([-0.14 0.06])
    xlim ([-100 300])
    title ('chevron (82)');
    ylabel('Airflow (au)');
    xlabel('Time (ms)');
end
yyaxis right
plot (centers,chevHist,'Color',tempCol,'LineWidth',3);
hold on
plot (centers,chevOffHist,'Color',tempCol,'LineWidth',3,'LineStyle',':');
ylim ([0 1])
yticklabels([]);
ax7 = subplot (2,5,7);
yline(0,':'); 
xline(0,':'); 
hold on
for i=1:twoStepNum
    tempBreath = twoStepBreaths{i,1};
    tempX = 1:length(tempBreath);
    tempX = tempX-twoStepInspDur(i);
    tempCol = c(colorInd(8),:);
    plot (tempX,tempBreath,'Color',tempCol);
    ylim ([-0.14 0.08])
    xlim ([-100 300])
    title ('two step (63)');
    xlabel('Time (ms)');
end
yticklabels([]);
yyaxis right
plot (centers,twoStepHist,'Color',tempCol,'LineWidth',3);
hold on
plot (centers,twoStepOffHist,'Color',tempCol,'LineWidth',3,'LineStyle',':');
ylim ([0 1])
yticklabels([]);
ax8 = subplot (2,5,8);
yline(0,':'); 
xline(0,':'); 
hold on
for i=1:multiNum
    tempBreath = multiBreaths{i,1};
    tempX = 1:length(tempBreath);
    tempX = tempX-multiInspDur(i);
    tempCol = c(colorInd(9),:);
    plot (tempX,tempBreath,'Color',tempCol);
    ylim ([-0.14 0.08])
    xlim ([-100 300])
    title ('multi (46)');
    xlabel('Time (ms)');
end
yticklabels([]);
yyaxis right
plot (centers,multiHist,'Color',tempCol,'LineWidth',3);
hold on
plot (centers,multiOffHist,'Color',tempCol,'LineWidth',3,'LineStyle',':');
ylim ([0 1])
yticklabels([]);
ax9 = subplot (2,5,9);
yline(0,':'); 
xline(0,':'); 
hold on
for i=1:stepUpNum
    tempBreath = stepUpBreaths{i,1};
    tempX = 1:length(tempBreath);
    tempX = tempX-stepUpInspDur(i);
    tempCol = c(colorInd(10),:);
    plot (tempX,tempBreath,'Color',tempCol);
    ylim ([-0.14 0.08])
    xlim ([-100 300])
    title ('step up (42)');
    xlabel('Time (ms)');
end
yticklabels([]);
yyaxis right
plot (centers,stepUpHist,'Color',tempCol,'LineWidth',3);
hold on
plot (centers,stepUpOffHist,'Color',tempCol,'LineWidth',3,'LineStyle',':');
ylim ([0 1])
yticklabels([]);
ax10 = subplot (2,5,10);
yline(0,':'); 
xline(0,':'); 
hold on
for i=1:dfmNum
    tempBreath = dfmBreaths{i,1};
    tempX = 1:length(tempBreath);
    tempX = tempX-dfmInspDur(i);
    tempCol = c(colorInd(11),:);
    plot (tempX,tempBreath,'Color',tempCol);
    ylim ([-0.14 0.08])
    xlim ([-100 300])
    title ('down fm (39)');
    xlabel('Time (ms)');
end
yticklabels([]);
yyaxis right
plot (centers,dfmHist,'Color',tempCol,'LineWidth',3);
hold on
plot (centers,dfmOffHist,'Color',tempCol,'LineWidth',3,'LineStyle',':');
ylim ([0 1])
ylabel('proportion of USV');

%% breath features
%Ti
boxInspDur = nan(10,length(ufmInspDur));
boxInspDur(1,:)=ufmInspDur;boxInspDur(2,1:length(stepDownInspDur))=stepDownInspDur;boxInspDur(3,1:length(flatInspDur))=flatInspDur;
boxInspDur(4,1:length(shortInspDur))=shortInspDur;boxInspDur(5,1:length(compInspDur))=compInspDur;
boxInspDur(6,1:length(chevInspDur))=chevInspDur;boxInspDur(7,1:length(twoStepInspDur))=twoStepInspDur;
boxInspDur(8,1:length(multiInspDur))=multiInspDur;boxInspDur(9,1:length(stepUpInspDur))=stepUpInspDur;boxInspDur(10,1:length(dfmInspDur))=dfmInspDur;
boxInspDur = transpose (boxInspDur);
figure
boxLabels = categorical({'up fm','step down','flat','short','complex','chevron','two step','multi step','step up','down fm'});
boxLabels = reordercats(boxLabels,{'up fm','step down','flat','short','complex','chevron','two step','multi step','step up','down fm'});
boxplot (boxInspDur,boxLabels,'symbol','');
ylim([0 100]);
ylabel('Ti(ms)');
%Te
boxExpDur = nan(10,length(ufmExpDur));
boxExpDur(1,:)=ufmExpDur;boxExpDur(2,1:length(stepDownExpDur))=stepDownExpDur;boxExpDur(3,1:length(flatExpDur))=flatExpDur;
boxExpDur(4,1:length(shortExpDur))=shortExpDur;boxExpDur(5,1:length(compExpDur))=compExpDur;
boxExpDur(6,1:length(chevExpDur))=chevExpDur;boxExpDur(7,1:length(twoStepExpDur))=twoStepExpDur;
boxExpDur(8,1:length(multiExpDur))=multiExpDur;boxExpDur(9,1:length(stepUpExpDur))=stepUpExpDur;boxExpDur(10,1:length(dfmExpDur))=dfmExpDur;
boxExpDur = transpose (boxExpDur);
figure
boxplot (boxExpDur,boxLabels,'symbol','');
ylim([0 300]);
ylabel('Te(ms)');
%Corr
boxCorr = nan(10,length(ufmCorr));
boxCorr(1,:)=ufmCorr;boxCorr(2,1:length(stepDownCorr))=stepDownCorr;boxCorr(3,1:length(flatCorr))=flatCorr;
boxCorr(4,1:length(shortCorr))=shortCorr;boxCorr(5,1:length(compCorr))=compCorr;
boxCorr(6,1:length(chevCorr))=chevCorr;boxCorr(7,1:length(twoStepCorr))=twoStepCorr;
boxCorr(8,1:length(multiCorr))=multiCorr;boxCorr(9,1:length(stepUpCorr))=stepUpCorr;boxCorr(10,1:length(dfmCorr))=dfmCorr;
boxCorr = transpose (boxCorr);
figure
boxplot (boxCorr,boxLabels,'symbol','');
ylim([-1 1]);
ylabel('r');
exportCorr = nan(11,length(corrs));
exportCorr(1,:) = corrs;
exportCorr(2:11,1:length(boxCorr)) = transpose(boxCorr);
exportCorr = transpose(exportCorr);
%% pull out expirations and peak frequencies for USV types
freqs = unis(:,9);
compExps = vocExps(compInd); compFreqs = freqs(compInd);
ufmExps = vocExps(ufmInd); ufmFreqs = freqs(ufmInd);
dfmExps = vocExps(dfmInd); dfmFreqs = freqs(dfmInd);
shortExps = vocExps(shortInd); shortFreqs = freqs(shortInd);
stepDownExps = vocExps(stepDownInd); stepDownFreqs = freqs(stepDownInd);
stepUpExps = vocExps(stepUpInd); stepUpFreqs = freqs(stepUpInd);
flatExps = vocExps(flatInd); flatFreqs = freqs(flatInd);
chevExps = vocExps(chevInd); chevFreqs = freqs(chevInd);
twoStepExps = vocExps(twoStepInd); twoStepFreqs = freqs(twoStepInd);
multiExps = vocExps(multiInd); multiFreqs = freqs(multiInd);

%% get bi stats
corrsBi1 = cell2mat(bis(:,8));
corrsBi2 = cell2mat(bis(:,13));
corrsBi = vertcat(corrsBi1, corrsBi2);
typesBi1 = bis(:,7);
typesBi2 = bis(:,12);
typesBi = vertcat(typesBi1,typesBi2);
%% Extract call types bi syllabic
flatNumBi=1;ufmNumBi=1;chevNumBi=1;dfmNumBi=1;revChevNumBi=1;
shortNumBi=1;compNumBi=1;stepUpNumBi=1;stepDownNumBi=1;twoStepNumBi=1;multiNumBi=1;
for i = 1:length(typesBi)
   tempVal = typesBi{i,1};
   if contains(tempVal,'flat') == 1
       flatIndBi(flatNumBi) = i;
       flatNumBi = flatNumBi +1;
   elseif contains(tempVal,'up_fm') == 1
       ufmIndBi(ufmNumBi) = i;
       ufmNumBi = ufmNumBi +1;
   elseif contains(tempVal,'chevron') == 1
       chevIndBi(chevNumBi) = i;
       chevNumBi = chevNumBi +1;
   elseif contains(tempVal,'down_fm') == 1
       dfmIndBi(dfmNumBi) = i;
       dfmNumBi = dfmNumBi +1;
   elseif contains(tempVal,'rev_chevron') == 1
       revChevIndBi(revChevNumBi) = i;
       revChevNumBi = revChevNumBi +1;
   elseif contains(tempVal,'short') == 1
       shortIndBi(shortNumBi) = i;
       shortNumBi = shortNumBi +1;
   elseif contains(tempVal,'complex') == 1
       compIndBi(compNumBi) = i;
       compNumBi = compNumBi +1;
   elseif contains(tempVal,'step_up') == 1
       stepUpIndBi(stepUpNumBi) = i;
       stepUpNumBi = stepUpNumBi +1;
   elseif contains(tempVal,'step_down') == 1
       stepDownIndBi(stepDownNumBi) = i;
       stepDownNumBi = stepDownNumBi +1;
   elseif contains(tempVal,'two_steps') == 1
       twoStepIndBi(twoStepNumBi) = i;
       twoStepNumBi = twoStepNumBi +1;
   elseif contains(tempVal,'mult_steps') == 1
       multiIndBi(multiNumBi) = i;
       multiNumBi = multiNumBi +1;
   end
end
flatNumBi = flatNumBi-1;ufmNumBi=ufmNumBi-1;chevNumBi=chevNumBi-1;dfmNumBi=dfmNumBi-1;
revChevNumBi = revChevNumBi-1;shortNumBi = shortNumBi-1;compNumBi = compNumBi-1;stepUpNumBi = stepUpNumBi-1;
stepDownNumBi = stepDownNumBi-1;twoStepNumBi =twoStepNumBi-1;multiNumBi =multiNumBi-1;
allNumsBi = [flatNumBi,ufmNumBi,chevNumBi,dfmNumBi,revChevNumBi,shortNumBi,compNumBi,stepUpNumBi,stepDownNumBi,twoStepNumBi,multiNumBi];
flatCorrBi = corrsBi(flatIndBi);ufmCorrBi = corrsBi(ufmIndBi);chevCorrBi = corrsBi(chevIndBi);dfmCorrBi = corrsBi(dfmIndBi);
revChevCorrBi = corrsBi(revChevIndBi);
shortCorrBi = corrsBi(shortIndBi);compCorrBi = corrsBi(compIndBi);stepUpCorrBi = corrsBi(stepUpIndBi);
stepDownCorrBi = corrsBi(stepDownIndBi);twoStepCorrBi = corrsBi(twoStepIndBi);multiCorrBi = corrsBi(multiIndBi);
%% get tri stats
corrsTri1 = cell2mat(tris(:,8));
corrsTri2 = cell2mat(tris(:,13));
corrsTri3 = cell2mat(tris(:,18));
corrsTri = vertcat(corrsTri1, corrsTri2,corrsTri3);
typesTri1 = tris(:,7);
typesTri2 = tris(:,12);
typesTri3 = tris(:,17);
typesTri = vertcat(typesTri1,typesTri2,typesTri3);
%% Extract call types bi syllabic
flatNumTri=1;ufmNumTri=1;chevNumTri=1;dfmNumTri=1;revChevNumTri=1;
shortNumTri=1;compNumTri=1;stepUpNumTri=1;stepDownNumTri=1;twoStepNumTri=1;multiNumTri=1;
for i = 1:length(typesTri)
   tempVal = typesTri{i,1};
   if contains(tempVal,'flat') == 1
       flatIndTri(flatNumTri) = i;
       flatNumTri = flatNumTri +1;
   elseif contains(tempVal,'up_fm') == 1
       ufmIndTri(ufmNumTri) = i;
       ufmNumTri = ufmNumTri +1;
   elseif contains(tempVal,'chevron') == 1
       chevIndTri(chevNumTri) = i;
       chevNumTri = chevNumTri +1;
   elseif contains(tempVal,'down_fm') == 1
       dfmIndTri(dfmNumTri) = i;
       dfmNumTri = dfmNumTri +1;
   elseif contains(tempVal,'rev_chevron') == 1
       revChevIndTri(revChevNumTri) = i;
       revChevNumTri = revChevNumTri +1;
   elseif contains(tempVal,'short') == 1
       shortIndTri(shortNumTri) = i;
       shortNumTri = shortNumTri +1;
   elseif contains(tempVal,'complex') == 1
       compIndTri(compNumTri) = i;
       compNumTri = compNumTri +1;
   elseif contains(tempVal,'step_up') == 1
       stepUpIndTri(stepUpNumTri) = i;
       stepUpNumTri = stepUpNumTri +1;
   elseif contains(tempVal,'step_down') == 1
       stepDownIndTri(stepDownNumTri) = i;
       stepDownNumTri = stepDownNumTri +1;
   elseif contains(tempVal,'two_steps') == 1
       twoStepIndTri(twoStepNumTri) = i;
       twoStepNumTri = twoStepNumTri +1;
   elseif contains(tempVal,'mult_steps') == 1
       multiIndTri(multiNumTri) = i;
       multiNumTri = multiNumTri +1;
   end
end
flatNumTri = flatNumTri-1;ufmNumTri=ufmNumTri-1;chevNumTri=chevNumTri-1;dfmNumTri=dfmNumTri-1;
revChevNumTri = revChevNumTri-1;shortNumTri = shortNumTri-1;compNumTri = compNumTri-1;stepUpNumTri = stepUpNumTri-1;
stepDownNumTri = stepDownNumTri-1;twoStepNumTri =twoStepNumTri-1;multiNumTri =multiNumTri-1;
allNumsTri = [flatNumTri,ufmNumTri,chevNumTri,dfmNumTri,revChevNumTri,shortNumTri,compNumTri,stepUpNumTri,stepDownNumTri,twoStepNumTri,multiNumTri];
flatCorrTri = corrsTri(flatIndTri);ufmCorrTri = corrsTri(ufmIndTri);chevCorrTri = corrsTri(chevIndTri);dfmCorrTri = corrsTri(dfmIndTri);
revChevCorrTri = corrsTri(revChevIndTri);
shortCorrTri = corrsTri(shortIndTri);compCorrTri = corrsTri(compIndTri);stepUpCorrTri = corrsTri(stepUpIndTri);
stepDownCorrTri = corrsTri(stepDownIndTri);twoStepCorrTri = corrsTri(twoStepIndTri);multiCorrTri = corrsTri(multiIndTri);%

%% Compile totals
totalNums = vertcat(allNums,allNumsBi,allNumsTri);
totalNums = sum(totalNums,1);
orderNums = sort(totalNums,'descend');
orderLabels = categorical({'up fm','flat','step down','short','complex','chevron','two step','down fm','step up','multi','rev chev'});
orderLabels = reordercats(orderLabels,{'up fm','flat','step down','short','complex','chevron','two step','down fm','step up','multi','rev chev'});
figure
bar (orderLabels,orderNums);
flatCorrs=vertcat(flatCorr,flatCorrBi,flatCorrTri);ufmCorrs=vertcat(ufmCorr,ufmCorrBi,ufmCorrTri);stepDownCorrs=vertcat(stepDownCorr,stepDownCorrBi,stepDownCorrTri);
shortCorrs=vertcat(shortCorr,shortCorrBi,shortCorrTri);compCorrs=vertcat(compCorr,compCorrBi,compCorrTri);chevCorrs=vertcat(chevCorr,chevCorrBi,chevCorrTri);
twoStepCorrs=vertcat(twoStepCorr,twoStepCorrBi,twoStepCorrTri);dfmCorrs=vertcat(dfmCorr,dfmCorrBi,dfmCorrTri);stepUpCorrs=vertcat(stepUpCorr,stepUpCorrBi,stepUpCorrTri);
multiCorrs=vertcat(multiCorr,multiCorrBi,multiCorrTri);
%Corr
boxCorrs = nan(10,length(ufmCorrs));
boxCorrs(1,:)=ufmCorrs;boxCorrs(2,1:length(flatCorrs))=flatCorrs;boxCorrs(3,1:length(stepDownCorrs))=stepDownCorrs;
boxCorrs(4,1:length(shortCorrs))=shortCorrs;boxCorrs(5,1:length(compCorrs))=compCorrs;
boxCorrs(6,1:length(chevCorrs))=chevCorrs;boxCorrs(7,1:length(twoStepCorrs))=twoStepCorrs;
boxCorrs(8,1:length(dfmCorrs))=dfmCorrs;boxCorrs(9,1:length(stepUpCorrs))=stepUpCorrs;boxCorrs(10,1:length(multiCorrs))=multiCorrs;
boxCorrs = transpose (boxCorrs);
boxLabels2 = categorical({'up fm','flat','step down','short','complex','chevron','two step','down fm','step up','multi step'});
boxLabels2 = reordercats(boxLabels2,{'up fm','flat','step down','short','complex','chevron','two step','down fm','step up','multi step'});
figure
boxplot (boxCorrs,boxLabels2,'symbol','');
ylim([-1 1]);
ylabel('r');
%% shuffle each set of pitch values and compare correlation
numPermutations = 1; %number of shuffled pitch vectors to simulate
normExps = unis(:,11); %extract upsampled exps
for vocNum = 1:length(freqs) %cycle through frequency vectors
    tempFreq = freqs{vocNum}; 
    tempExp = normExps{vocNum};
    tempCorr = corrs(vocNum);
    permMat = zeros(numPermutations,length(tempFreq)); %make somewhere to store the simulated vectors
    for i = 1:numPermutations %create shuffled vectors
        randInds = randperm(length(tempFreq));
        permMat(i,:) = tempFreq(randInds);
    end
    for i = 1:numPermutations %for each vector calculate R and store
        tempPerm = permMat(i,:);
        r = corrcoef(tempExp,tempPerm,'Rows','complete');
        randCorrs(i) = r(2);
    end
    %find out how many SDs real value is from mean of simulated
    %distribution
    zCorrs(vocNum) = (tempCorr-mean(randCorrs))/std(randCorrs);
    simCorrs(vocNum) = randCorrs(1);
end
upZ = zCorrs(ufmInd); stepDownZ = zCorrs(stepDownInd); flatZ = zCorrs(flatInd);
shortZ = zCorrs(shortInd); compZ = zCorrs(compInd); chevZ = zCorrs(chevInd);
twoStepZ = zCorrs(twoStepInd);multiZ = zCorrs(multiInd);stepUpZ=zCorrs(stepUpInd);
dfmZ = zCorrs(dfmInd);
boxZ = nan(10,length(upZ));
boxZ(1,:) = upZ;boxZ(2,1:length(stepDownZ)) = stepDownZ;boxZ(3,1:length(flatZ)) = flatZ;
boxZ(4,1:length(shortZ)) = shortZ; boxZ(5,1:length(compZ)) = compZ; boxZ(6,1:length(chevZ)) = chevZ;
boxZ(7,1:length(twoStepZ)) = twoStepZ;boxZ(8,1:length(multiZ)) = multiZ;
boxZ(9,1:length(stepUpZ)) = stepUpZ; boxZ(10,1:length(dfmZ)) = dfmZ;
boxZ = transpose(boxZ);
figure
boxplot (boxZ,boxLabels,'symbol','');
hold on
yline(0,':');
yline(-1.96,':');
yline(1.96,':');
yline(-2.58,':');
yline(2.58,':');
ylabel('z-score');
title('pitch shuffled');
% make boxplot bars for shuffled data
upSim = simCorrs(ufmInd); stepDownSim = simCorrs(stepDownInd); flatSim = simCorrs(flatInd);
shortSim = simCorrs(shortInd); compSim = simCorrs(compInd); chevSim = simCorrs(chevInd);
twoStepSim = simCorrs(twoStepInd);multiSim = simCorrs(multiInd);stepUpSim=simCorrs(stepUpInd);
dfmSim = simCorrs(dfmInd);
boxSim = nan(10,length(upSim));
boxSim(1,:) = upSim;boxSim(2,1:length(stepDownSim)) = stepDownSim;boxSim(3,1:length(flatSim)) = flatSim;
boxSim(4,1:length(shortSim)) = shortSim; boxSim(5,1:length(compSim)) = compSim; boxSim(6,1:length(chevSim)) = chevSim;
boxSim(7,1:length(twoStepSim)) = twoStepSim;boxSim(8,1:length(multiSim)) = multiSim;
boxSim(9,1:length(stepUpSim)) = stepUpSim; boxSim(10,1:length(dfmSim)) = dfmSim;
boxSim = transpose(boxSim);
figure
boxplot (boxSim,boxLabels,'symbol','');
ylim([-1 1]);
hold on
yline(0,':');
ylabel('r');
title('pitch shuffled');
%% shuffle each set of breath values and compare correlation
numPermutations = 1; %number of shuffled pitch vectors to simulate
normExps = unis(:,11); %extract upsampled exps
for vocNum = 1:length(freqs) %cycle through frequency vectors
    tempFreq = freqs{vocNum}; 
    tempExp = normExps{vocNum};
    tempCorr = corrs(vocNum);
    permMat = zeros(numPermutations,length(tempExp)); %make somewhere to store the simulated vectors
    for i = 1:numPermutations %create shuffled vectors
        randInds = randperm(length(tempExp));
        permMat(i,:) = tempExp(randInds);
    end
    for i = 1:numPermutations %for each vector calculate R and store
        tempPerm = permMat(i,:);
        r = corrcoef(tempPerm,tempExp,'Rows','complete');
        randCorrs(i) = r(2);
    end
    %find out how many SDs real value is from mean of simulated
    %distribution
    zCorrs(vocNum) = (tempCorr-mean(randCorrs))/std(randCorrs);
    simCorrs(vocNum) = randCorrs(1);
end
upZ = zCorrs(ufmInd); stepDownZ = zCorrs(stepDownInd); flatZ = zCorrs(flatInd);
shortZ = zCorrs(shortInd); compZ = zCorrs(compInd); chevZ = zCorrs(chevInd);
twoStepZ = zCorrs(twoStepInd);multiZ = zCorrs(multiInd);stepUpZ=zCorrs(stepUpInd);
dfmZ = zCorrs(dfmInd);
boxZ = nan(10,length(upZ));
boxZ(1,:) = upZ;boxZ(2,1:length(stepDownZ)) = stepDownZ;boxZ(3,1:length(flatZ)) = flatZ;
boxZ(4,1:length(shortZ)) = shortZ; boxZ(5,1:length(compZ)) = compZ; boxZ(6,1:length(chevZ)) = chevZ;
boxZ(7,1:length(twoStepZ)) = twoStepZ;boxZ(8,1:length(multiZ)) = multiZ;
boxZ(9,1:length(stepUpZ)) = stepUpZ; boxZ(10,1:length(dfmZ)) = dfmZ;
boxZ = transpose(boxZ);
figure
boxplot (boxZ,boxLabels,'symbol','');
hold on
yline(0,':');
yline(-1.96,':');
yline(1.96,':');
yline(-2.58,':');
yline(2.58,':');
ylabel('z-score');
title('airflow shuffled');
% make boxplot bars for shuffled data
upSim = simCorrs(ufmInd); stepDownSim = simCorrs(stepDownInd); flatSim = simCorrs(flatInd);
shortSim = simCorrs(shortInd); compSim = simCorrs(compInd); chevSim = simCorrs(chevInd);
twoStepSim = simCorrs(twoStepInd);multiSim = simCorrs(multiInd);stepUpSim=simCorrs(stepUpInd);
dfmSim = simCorrs(dfmInd);
boxSim = nan(10,length(upSim));
boxSim(1,:) = upSim;boxSim(2,1:length(stepDownSim)) = stepDownSim;boxSim(3,1:length(flatSim)) = flatSim;
boxSim(4,1:length(shortSim)) = shortSim; boxSim(5,1:length(compSim)) = compSim; boxSim(6,1:length(chevSim)) = chevSim;
boxSim(7,1:length(twoStepSim)) = twoStepSim;boxSim(8,1:length(multiSim)) = multiSim;
boxSim(9,1:length(stepUpSim)) = stepUpSim; boxSim(10,1:length(dfmSim)) = dfmSim;
boxSim = transpose(boxSim);
figure
boxplot (boxSim,boxLabels,'symbol','');
hold on
yline(0,':');
ylim([-1 1]);
hold on
ylabel('r');
title('airflow shuffled');