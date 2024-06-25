%% Compile all calls raster type
% Compile call info from multiple recordings (analyzed by callVarsPreProcessing.m)
% into one dataset and calculate statistics sorted by call type Compute raster and onset offset histograms
%% Load metadata and extract filenames
pathToData = '';
cd(pathToData)
metaDataFile = [pathToData,'compdata.xlsx'];
metaData = readtable (metaDataFile,'ReadVariableNames',true);
metaData = table2struct(metaData);

%% Import call vars and store
unis = {'breath trace','long trace','inspStart','expStart','vocStart','vocEnd','vocClass','corr','norm freq trace','abs freq trace','norm upsampled exp', 'abs upsampled exp'};
%unis = {'breath trace','long trace','inspStart','expStart','vocStart','vocEnd','vocClass','corr','freq trace'};
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
clear uniCell; clear biCell; clear triCell;
unis = unis(2:end,1:9); bis = bis(2:end,:); tris = tris(2:end,:);
%get uni stats
expStarts = cell2mat(unis(:,4));
inspStarts = cell2mat(unis(:,3));
vocStarts = cell2mat(unis(:,5));
vocEnds = cell2mat(unis(:,6));
vocClass = unis(:,7);
vocOnsets = vocStarts-expStarts;
vocOffsets = vocEnds-expStarts;
breaths =transpose(unis(:,2));
breaths = cell2mat(breaths);
breaths = transpose(breaths);
for i = 1:length(unis)
    breathLength(i,1) = length(unis{i,1});
end
%breathLength=transpose(breathLength);
inspDur = expStarts-inspStarts;
expDur = breathLength-inspDur;
%% get bi stats
biBreaths = transpose(bis(:,2));
biBreaths = cell2mat(biBreaths);
biBreaths = transpose(biBreaths);
biInspStarts = cell2mat(bis(:,3));
biExpStarts = cell2mat(bis(:,4));
vocStarts1 = cell2mat(bis(:,5));
vocEnds1 =  cell2mat(bis(:,6));
vocStarts2 = cell2mat(bis(:,10));
vocEnds2 = cell2mat(bis(:,11));
vocOnsets1 = vocStarts1-biExpStarts;
vocOffsets1 = vocEnds1-biExpStarts;
vocOnsets2 = vocStarts2-biExpStarts;
vocOffsets2 = vocEnds2-biExpStarts;
for i = 1:length(bis)
    biBreathLength(i,1) = length(bis{i,1});
end
biInspDur = biExpStarts-biInspStarts;
biExpDur = biBreathLength-biInspDur;
%% find call types for hist
flatNum=1;ufmNum=1;chevNum=1;dfmNum=1;revChevNum=1;
shortNum=1;compNum=1;stepUpNum=1;stepDownNum=1;twoStepNum=1;multiNum = 1;
for i = 1:length(vocClass)
    tempVal = vocClass{i,1};
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
edges = [-100:10:300];
for i = 1:length(edges)-1
    centers(i) = (edges(i)+(edges(i+1)))/2;
end
normEdges = [0:0.05:1];
allInd = 1:length(unis);
%hists and ordered breaths and stats by type
[allOnHist,allOffHist,allOn,allOff,allBreaths,allNormOnHist,allNormOffHist,allNormOn,allNormOff] = getOnsetOffset(breaths,vocOnsets,vocOffsets,edges,allInd,expDur,normEdges);
[ufmOnHist,ufmOffHist,ufmOn,ufmOff,ufmBreaths,ufmNormOnHist,ufmNormOffHist,ufmNormOn,ufmNormOff] = getOnsetOffset(breaths,vocOnsets,vocOffsets,edges,ufmInd,expDur,normEdges);
[stepDownOnHist,stepDownOffHist,stepDownOn,stepDownOff,stepDownBreaths,stepDownNormOnHist,stepDownNormOffHist,stepDownNormOn,stepDownNormOff] = getOnsetOffset(breaths,vocOnsets,vocOffsets,edges,stepDownInd,expDur,normEdges);
[shortOnHist,shortOffHist,shortOn,shortOff,shortBreaths,shortNormOnHist,shortNormOffHist,shortNormOn,shortNormOff] = getOnsetOffset(breaths,vocOnsets,vocOffsets,edges,shortInd,expDur,normEdges);
[flatOnHist,flatOffHist,flatOn,flatOff,flatBreaths,flatNormOnHist,flatNormOffHist,flatNormOn,flatNormOff] = getOnsetOffset(breaths,vocOnsets,vocOffsets,edges,flatInd,expDur,normEdges);
[compOnHist,compOffHist,compOn,compOff,compBreaths,compNormOnHist,compNormOffHist,compNormOn,compNormOff] = getOnsetOffset(breaths,vocOnsets,vocOffsets,edges,compInd,expDur,normEdges);
[chevOnHist,chevOffHist,chevOn,chevOff,chevBreaths,chevNormOnHist,chevNormOffHist,chevNormOn,chevNormOff] = getOnsetOffset(breaths,vocOnsets,vocOffsets,edges,chevInd,expDur,normEdges);
[twoStepOnHist,twoStepOffHist,twoStepOn,twoStepOff,twoStepBreaths,twoStepNormOnHist,twoStepNormOffHist,twoStepNormOn,twoStepNormOff] = getOnsetOffset(breaths,vocOnsets,vocOffsets,edges,twoStepInd,expDur,normEdges);
[multiOnHist,multiOffHist,multiOn,multiOff,multiBreaths,multiNormOnHist,multiNormOffHist,multiNormOn,multiNormOff] = getOnsetOffset(breaths,vocOnsets,vocOffsets,edges,multiInd,expDur,normEdges);
[stepUpOnHist,stepUpOffHist,stepUpOn,stepUpOff,stepUpBreaths,stepUpNormOnHist,stepUpNormOffHist,stepUpNormOn,stepUpNormOff] = getOnsetOffset(breaths,vocOnsets,vocOffsets,edges,stepUpInd,expDur,normEdges);
[dfmOnHist,dfmOffHist,dfmOn,dfmOff,dfmBreaths,dfmNormOnHist,dfmNormOffHist,dfmNormOn,dfmNormOff] = getOnsetOffset(breaths,vocOnsets,vocOffsets,edges,dfmInd,expDur,normEdges);

%calculate for bis
%remove USVs misattributed to the preceeding breath (onset<0)
preVocInd = vocOnsets1>0;
vocOnsets1 = vocOnsets1(preVocInd);vocOnsets2 = vocOnsets2(preVocInd);
vocOffsets1 = vocOffsets1(preVocInd);vocOffsets2 = vocOffsets2(preVocInd);
biExpDur = biExpDur(preVocInd);
biBreaths = biBreaths(preVocInd,:);
biNormOnsets1 = vocOnsets1./biExpDur;biNormOnsets2 = vocOnsets2./biExpDur;
biNormOffsets1 = vocOffsets1./biExpDur;biNormOffsets2 = vocOffsets2./biExpDur;
%order breaths in increasing voc1 Onset
[vocOnsets1, orderInd] = sort(vocOnsets1);
vocOnsets2 = vocOnsets2(orderInd);
vocOffsets1 = vocOffsets1(orderInd); vocOffsets2 = vocOffsets2(orderInd);
biBreaths = biBreaths(orderInd,:);
[bi1Hist ~] = histcounts(vocOnsets1, edges);[bi2Hist ~] = histcounts(vocOnsets2, edges);
[bi1OffHist ~] = histcounts(vocOffsets1, edges);[bi2OffHist ~] = histcounts(vocOffsets2, edges);
[bi1NormHist ~] = histcounts(biNormOnsets1, normEdges);[bi2NormHist ~] = histcounts(biNormOnsets2, normEdges);
[bi1NormOffHist ~] = histcounts(biNormOffsets1, normEdges);[bi2NormOffHist ~] = histcounts(biNormOffsets2, normEdges);
bi1NormHist = bi1NormHist/length(vocOnsets1); bi1NormOffHist = bi1NormOffHist/length(vocOnsets1);
bi2NormHist = bi2NormHist/length(vocOnsets1); bi2NormOffHist = bi2NormOffHist/length(vocOnsets1);
%figure time
%figure
%subplot(3,4,1)
%plot (edges(2:end),ufmOnHist,'Color','k','LineWidth',2);
%hold on
%plot (edges(2:end),ufmOffHist, 'Color','r','LineWidth',2);
%ylabel('#USVs')
%xlabel('Time (ms)')
%title('up fm')
%subplot(3,4,2)
%plot (edges(2:end),stepDownOnHist,'Color','k','LineWidth',2);
%hold on
%plot (edges(2:end),stepDownOffHist, 'Color','r','LineWidth',2);
%ylabel('#USVs')
%xlabel('Time (ms)')
%title('step down')
%subplot(3,4,3)
%plot (edges(2:end),flatOnHist,'Color','k','LineWidth',2);
%hold on
%plot (edges(2:end),flatOffHist, 'Color','r','LineWidth',2);
%ylabel('#USVs')
%xlabel('Time (ms)')
%title('flat')
%subplot(3,4,4)
%plot (edges(2:end),shortOnHist,'Color','k','LineWidth',2);
%hold on
%plot (edges(2:end),shortOffHist, 'Color','r','LineWidth',2);
%ylabel('#USVs')
%xlabel('Time (ms)')
%title('short')
%subplot(3,4,5)
%plot (edges(2:end),compOnHist,'Color','k','LineWidth',2);
%hold on
%plot (edges(2:end),compOffHist, 'Color','r','LineWidth',2);
%ylabel('#USVs')
%xlabel('Time (ms)')
%title('complex')
%subplot(3,4,6)
%plot (edges(2:end),chevOnHist,'Color','k','LineWidth',2);
%hold on
%plot (edges(2:end),chevOffHist, 'Color','r','LineWidth',2);
%ylabel('#USVs')
%xlabel('Time (ms)')
%title('chevron')
%subplot(3,4,7)
%plot (edges(2:end),twoStepOnHist,'Color','r','LineWidth',2);
%hold on
%plot (edges(2:end),twoStepOffHist, 'Color','r','LineWidth',2);
%ylabel('#USVs')
%xlabel('Time (ms)')
%title('two step')
%subplot(3,4,8)
%plot (edges(2:end),multiOnHist,'Color','k','LineWidth',2);
%hold on
%plot (edges(2:end),multiOffHist, 'Color','r','LineWidth',2);
%ylabel('#USVs')
%xlabel('Time (ms)')
%title('multi')
%subplot(3,4,9)
%plot (edges(2:end),stepUpOnHist,'Color','k','LineWidth',2);
%hold on
%plot (edges(2:end),stepUpOffHist, 'Color','r','LineWidth',2);
%ylabel('#USVs')
%xlabel('Time (ms)')
%title('step up')
%subplot(3,4,10)
%plot (edges(2:end),dfmOnHist,'Color','k','LineWidth',2);
%hold on
%plot (edges(2:end),dfmOffHist, 'Color','r','LineWidth',2);
%ylabel('#USVs')
%xlabel('Time (ms)')
%title('dfm')
%subplot(3,4,11)
%plot (edges(2:end),bi1Hist,'Color','k','LineWidth',2);
%hold on
%plot (edges(2:end),bi1OffHist, 'Color','r','LineWidth',2);
%plot (edges(2:end),bi2Hist,'Color','k','LineWidth',2);
%plot (edges(2:end),bi2OffHist, 'Color','r','LineWidth',2);
%ylabel('#USVs')
%xlabel('Time (ms)')
%title('bisyllabic')
%%
c1 = [0.0 0.2 0.5];
c2= [1 1 1];
c3 = [0.5 0.5 0.5];
c = generateCmap(c1,c2,c3);
makeHeatRaster(ufmBreaths,ufmOn,ufmOff,c,true);
title('up fm');
makeHeatRaster(stepDownBreaths,stepDownOn,stepDownOff,c,true);
title('step down');
makeHeatRaster(flatBreaths,flatOn,flatOff,c,true);
title('flat');
makeHeatRaster(shortBreaths,shortOn,shortOff,c,true);
title('short');
makeHeatRaster(compBreaths,compOn,compOff,c,true);
title('complex')
makeHeatRaster(chevBreaths,chevOn,chevOff,c,true);
title('chevron');
makeHeatRaster(twoStepBreaths,twoStepOn,twoStepOff,c,true);
title('two step');
makeHeatRaster(multiBreaths,multiOn,multiOff,c,true);
title('multi step')
makeHeatRaster(stepUpBreaths,stepUpOn,stepUpOff,c,true);
title('step up');
makeHeatRaster(dfmBreaths,dfmOn,dfmOff,c,true);
title('down fm');
makeHeatRaster(allBreaths,allOn,allOff,c,true);
title('all');
figure
h=imagesc(biBreaths);
colormap(c);
j = colorbar;
xlim([900 1300]);
%find values for color axis
minval = min(biBreaths,[],'all');
maxval = -minval;
caxis([minval maxval]);
plotVocOn = vocOnsets1 +1000;
plotVocOff = vocOffsets1 +1000;
hold on
scatter (plotVocOn, 1:length(plotVocOn),0.16,"_","MarkerEdgeColor",'k');
scatter (plotVocOff, 1:length(plotVocOff),0.16,"_","MarkerEdgeColor",'r');
plotVocOn = vocOnsets2 +1000;
plotVocOff = vocOffsets2 +1000;
scatter (plotVocOn, 1:length(plotVocOn),0.16,"_","MarkerEdgeColor",'k');
scatter (plotVocOff, 1:length(plotVocOff),0.16,"_","MarkerEdgeColor",'r');
title('bisyllabic');
close all
%figure time
figure
subplot(3,4,1)
plot (normEdges,[0 ufmNormOnHist],'Color','k','LineWidth',2);
hold on
plot (normEdges,[0 ufmNormOffHist], 'Color','r','LineWidth',2);
ylabel('prop. USVs')
xlabel('exp duration (norm)')
title('up fm')
ylim([0 0.4]);
subplot(3,4,2)
plot (normEdges,[0 stepDownNormOnHist],'Color','k','LineWidth',2);
hold on
plot (normEdges,[0 stepDownNormOffHist], 'Color','r','LineWidth',2);
ylabel('prop. USVs')
xlabel('exp duration (norm)')
title('step down')
ylim([0 0.4]);
subplot(3,4,3)
plot (normEdges,[0 flatNormOnHist],'Color','k','LineWidth',2);
hold on
plot (normEdges,[0 flatNormOffHist], 'Color','r','LineWidth',2);
ylabel('prop. USVs')
xlabel('exp duration (norm)')
title('flat')
ylim([0 0.4]);
subplot(3,4,4)
plot (normEdges,[0 shortNormOnHist],'Color','k','LineWidth',2);
hold on
plot (normEdges,[0 shortNormOffHist], 'Color','r','LineWidth',2);
ylabel('prop. USVs')
xlabel('exp duration (norm)')
title('short')
ylim([0 0.4]);
subplot(3,4,5)
plot (normEdges,[0 compNormOnHist],'Color','k','LineWidth',2);
hold on
plot (normEdges,[0 compNormOffHist], 'Color','r','LineWidth',2);
ylabel('prop. USVs')
xlabel('exp duration (norm)')
title('complex')
ylim([0 0.4]);
subplot(3,4,6)
plot (normEdges,[0 chevNormOnHist],'Color','k','LineWidth',2);
hold on
plot (normEdges,[0 chevNormOffHist], 'Color','r','LineWidth',2);
ylabel('prop. USVs')
xlabel('exp duration (norm)')
title('chevron')
ylim([0 0.4]);
subplot(3,4,7)
plot (normEdges,[0 twoStepNormOnHist],'Color','k','LineWidth',2);
hold on
plot (normEdges,[0 twoStepNormOffHist], 'Color','r','LineWidth',2);
ylabel('prop. USVs')
xlabel('exp duration (norm)')
title('two step')
ylim([0 0.4]);
subplot(3,4,8)
plot (normEdges,[0 multiNormOnHist],'Color','k','LineWidth',2);
hold on
plot (normEdges,[0 multiNormOffHist], 'Color','r','LineWidth',2);
ylabel('prop. USVs')
xlabel('exp duration (norm)')
title('multi')
ylim([0 0.4]);
subplot(3,4,9)
plot (normEdges,[0 stepUpNormOnHist],'Color','k','LineWidth',2);
hold on
plot (normEdges,[0 stepUpNormOffHist], 'Color','r','LineWidth',2);
ylabel('prop. USVs')
xlabel('exp duration (norm)')
title('step up')
ylim([0 0.4]);
subplot(3,4,10)
plot (normEdges,[0 dfmNormOnHist],'Color','k','LineWidth',2);
hold on
plot (normEdges,[0 dfmNormOffHist], 'Color','r','LineWidth',2);
ylabel('prop. USVs')
xlabel('exp duration (norm)')
title('dfm')
ylim([0 0.4]);
subplot(3,4,11)
plot (normEdges,[0 bi1NormHist],'Color','k','LineWidth',2);
hold on
plot (normEdges,[0 bi1NormOffHist], 'Color','r','LineWidth',2);
plot (normEdges,[0 bi2NormHist],'Color','k','LineWidth',2);
plot (normEdges,[0 bi2NormOffHist], 'Color','r','LineWidth',2);
ylabel('#USVs')
xlabel('Time (ms)')
title('bisyllabic')
ylim([0 0.4]);
subplot(3,4,12)
plot (normEdges,[0 allNormOnHist],'Color','k','LineWidth',2);
hold on
plot (normEdges,[0 allNormOffHist], 'Color','r','LineWidth',2);
ylabel('prop. USVs')
xlabel('exp duration (norm)')
title('all')
ylim([0 0.4]);
%% remove normalized values >1 (exp End misidentified)
ufmNormOn = ufmNormOn(ufmNormOff <= 1);ufmNormOff = ufmNormOff(ufmNormOff <= 1);
stepDownNormOn = stepDownNormOn(stepDownNormOff <= 1);stepDownNormOff = stepDownNormOff(stepDownNormOff <= 1);
flatNormOn = flatNormOn(flatNormOff <= 1);flatNormOff = flatNormOff(flatNormOff <= 1);
shortNormOn = shortNormOn(shortNormOff <= 1);shortNormOff = shortNormOff(shortNormOff <= 1);
compNormOn = compNormOn(compNormOff <= 1);compNormOff = compNormOff(compNormOff <= 1);
chevNormOn = chevNormOn(chevNormOff <= 1);chevNormOff = chevNormOff(chevNormOff <= 1);
twoStepNormOn = twoStepNormOn(twoStepNormOff <= 1);twoStepNormOff = twoStepNormOff(twoStepNormOff <= 1);
multiNormOn = multiNormOn(multiNormOff <= 1);multiNormOff = multiNormOff(multiNormOff <= 1);
stepUpNormOn = stepUpNormOn(stepUpNormOff <= 1);stepUpNormOff = stepUpNormOff(stepUpNormOff <= 1);
dfmNormOn = dfmNormOn(dfmNormOff <= 1);dfmNormOff = dfmNormOff(dfmNormOff <= 1);
%% box plot
%onsets
boxOn = nan(13,length(allOn));
boxOn(1,:)=allOn;boxOn(2,1:length(ufmOn))=ufmOn;boxOn(3,1:length(stepDownOn))=stepDownOn;boxOn(4,1:length(flatOn))=flatOn;
boxOn(5,1:length(shortOn))=shortOn;boxOn(6,1:length(compOn))=compOn;boxOn(7,1:length(chevOn))=chevOn;
boxOn(8,1:length(twoStepOn))=twoStepOn;boxOn(9,1:length(multiOn))=multiOn;boxOn(10,1:length(stepUpOn))=stepUpOn;
boxOn(11,1:length(dfmOn))=dfmOn;boxOn(12,1:length(vocOnsets1))=vocOnsets1;boxOn(13,1:length(vocOnsets2))=vocOnsets2;
boxOn = transpose(boxOn);
boxLabels = categorical({'all','up fm','step down','flat','short','complex','chevron','two step','multi step','step up','down fm','bi1','bi2'});
boxLabels = reordercats(boxLabels,{'all','up fm','step down','flat','short','complex','chevron','two step','multi step','step up','down fm','bi1','bi2'});
figure
boxplot (boxOn,boxLabels,'symbol','');
title('absolute onset');
ylabel('onset (ms)');
ylim([0 250]);
%offsets
boxOff = nan(13,length(allOff));
boxOff(1,:)=allOff;boxOff(2,1:length(ufmOff))=ufmOff;boxOff(3,1:length(stepDownOff))=stepDownOff;boxOff(4,1:length(flatOff))=flatOff;
boxOff(5,1:length(shortOff))=shortOff;boxOff(6,1:length(compOff))=compOff;boxOff(7,1:length(chevOff))=chevOff;
boxOff(8,1:length(twoStepOff))=twoStepOff;boxOff(9,1:length(multiOff))=multiOff;boxOff(10,1:length(stepUpOff))=stepUpOff;
boxOff(11,1:length(dfmOff))=dfmOff;boxOff(12,1:length(vocOffsets1))=vocOffsets1;boxOff(13,1:length(vocOffsets2))=vocOffsets2;
boxOff = transpose(boxOff);
figure
boxplot (boxOff,boxLabels,'symbol','');
title('absolute offset');
ylabel('offset (ms)');
ylim([0 250]);
%norm onsets
boxNormOn = nan(13,length(allNormOn));
boxNormOn(1,:)=allNormOn;boxNormOn(2,1:length(ufmNormOn))=ufmNormOn;boxNormOn(3,1:length(stepDownNormOn))=stepDownNormOn;boxNormOn(4,1:length(flatNormOn))=flatNormOn;
boxNormOn(5,1:length(shortNormOn))=shortNormOn;boxNormOn(6,1:length(compNormOn))=compNormOn;boxNormOn(7,1:length(chevNormOn))=chevNormOn;
boxNormOn(8,1:length(twoStepNormOn))=twoStepNormOn;boxNormOn(9,1:length(multiNormOn))=multiNormOn;boxNormOn(10,1:length(stepUpNormOn))=stepUpNormOn;
boxNormOn(11,1:length(dfmNormOn))=dfmNormOn;boxNormOn(12,1:length(biNormOnsets1))=biNormOnsets1;boxNormOn(13,1:length(biNormOnsets2))=biNormOnsets2;
boxNormOn = transpose(boxNormOn);
figure
boxplot (boxNormOn,boxLabels,'symbol','');
ylim([0 1]);
title('norm onset');
ylabel('onset (norm)');
%norm offsets
boxNormOff = nan(13,length(allNormOff));
boxNormOff(1,:)=allNormOff;boxNormOff(2,1:length(ufmNormOff))=ufmNormOff;boxNormOff(3,1:length(stepDownNormOff))=stepDownNormOff;boxNormOff(4,1:length(flatNormOff))=flatNormOff;
boxNormOff(5,1:length(shortNormOff))=shortNormOff;boxNormOff(6,1:length(compNormOff))=compNormOff;boxNormOff(7,1:length(chevNormOff))=chevNormOff;
boxNormOff(8,1:length(twoStepNormOff))=twoStepNormOff;boxNormOff(9,1:length(multiNormOff))=multiNormOff;boxNormOff(10,1:length(stepUpNormOff))=stepUpNormOff;
boxNormOff(11,1:length(dfmNormOff))=dfmNormOff;boxNormOff(12,1:length(biNormOffsets1))=biNormOffsets1;boxNormOff(13,1:length(biNormOffsets2))=biNormOffsets2;
boxNormOff = transpose(boxNormOff);
figure
boxplot (boxNormOff,boxLabels,'symbol','');
ylim([0 1]);
title('norm offset');
ylabel('offset (norm)');
%%
figure
scatter(ufmNormOn,ufmNormOff);
hold on
scatter (dfmNormOn,dfmNormOff);
scatter (compNormOn, compNormOff);
scatter (chevNormOn,chevNormOff);
xlim([0 1]);
ylim([0 1]);
xlabel('norm Onset');
ylabel('normOffset');