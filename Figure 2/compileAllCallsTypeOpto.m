%% Compile all calls type
% Compile call info from multiple recordings (analyzed by callVarsPreProcessing.m)
% into one dataset and calculate statistics sorted by call type
%% Load metadata and extract filenames
pathToData = 'C:\Users\Alastair\Box\USV Behavior\Adult USV Opto\Reachr\';
cd(pathToData)
metaDataFile = [pathToData,'metadata.xlsx'];
metaData = readtable (metaDataFile,'ReadVariableNames',true);
metaData = table2struct(metaData);

%% Import call vars and store
unis = {'breath trace','long trace','inspStart','expStart','vocStart','vocEnd','vocClass','corr','freq trace'};
bis = {'breath trace','long trace','inspStart','expStart','vocStart1','vocEnd1','vocClass1','corr1','freq trace 1','vocStart2','vocEnd2','vocClass2','corr2','freq trace 2'};
tris = {'breath trace','long trace','inspStart','expStart','vocStart1','vocEnd1','vocClass1','corr1','freq trace 1','vocStart2','vocEnd2','vocClass2','corr2','freq trace 2','vocStart3','vocEnd3','vocClass3','corr3','freq trace 3'};
for i = 24:26%length(metaData)
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
for i = 1:length(unis)
    tempBreath = breaths{i,1};
    breathLength(i,1) = length(tempBreath);
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
if flatNum == 0
    flatInd = 0;
elseif ufmNum == 0
    ufmInd = 0;
elseif chevNum == 0
    chevInd = 0;
elseif dfmNum == 0
    dfmInd = 0;
elseif revChevNum == 0 
    revChevInd = 0;
elseif shortNum == 0
    shortInd = 0;
elseif compNum == 0
    compInd = 0;
elseif stepUpNum == 0
    stepUpInd = 0;
elseif stepDownNum ==0
    stepDownInd = 0;
elseif twoStepNum == 0
    twoStepInd = 0;
elseif multiNum == 0 
    multiInd = 0
end
flatNum = flatNum-1;ufmNum = ufmNum-1;chevNum =chevNum-1;dfmNum =dfmNum-1;
revChevNum = revChevNum-1;shortNum = shortNum-1;compNum = compNum-1;stepUpNum = stepUpNum-1;
stepDownNum = stepDownNum-1;twoStepNum =twoStepNum-1;multiNum =multiNum-1;
allNums = [flatNum,ufmNum,chevNum,dfmNum,shortNum,compNum,stepUpNum,stepDownNum,twoStepNum];
sortNums = sort(allNums,'descend');
propNums = allNums/sum(allNums);
figure (1)
labels = categorical({'down fm','step down','flat','short','step up','up fm','chevron','two step','complex'});
labels = reordercats(labels,{'down fm','step down','flat','short','step up','up fm','chevron','two step','complex'});
bar (labels,sortNums);
ufmInspDur=inspDur(ufmInd);ufmExpDur=expDur(ufmInd);ufmBreaths=breaths(ufmInd); ufmVocOnsets=vocOnsets(ufmInd); ufmVocOffsets=vocOffsets(ufmInd); ufmCorr = corrs(ufmInd);
chevInspDur=inspDur(chevInd);chevExpDur=expDur(chevInd);chevBreaths=breaths(chevInd);chevVocOnsets=vocOnsets(chevInd);chevVocOffsets=vocOffsets(chevInd); chevCorr = corrs(chevInd);
twoStepInspDur=inspDur(twoStepInd);twoStepExpDur=expDur(twoStepInd);twoStepBreaths=breaths(twoStepInd);twoStepVocOnsets=vocOnsets(twoStepInd);twoStepVocOffsets=vocOffsets(twoStepInd); twoStepCorr = corrs(twoStepInd);
stepDownInspDur=inspDur(stepDownInd);stepDownExpDur=expDur(stepDownInd);stepDownBreaths=breaths(stepDownInd);stepDownVocOnsets=vocOnsets(stepDownInd);stepDownVocOffsets=vocOffsets(stepDownInd); stepDownCorr = corrs(stepDownInd);
flatInspDur=inspDur(flatInd);flatExpDur=expDur(flatInd);flatBreaths=breaths(flatInd);flatVocOnsets=vocOnsets(flatInd);flatVocOffsets=vocOffsets(flatInd); flatCorr = corrs(flatInd);
dfmInspDur=inspDur(dfmInd);dfmExpDur=expDur(dfmInd);dfmBreaths=breaths(dfmInd);dfmVocOnsets=vocOnsets(dfmInd);dfmVocOffsets=vocOffsets(dfmInd); dfmCorr = corrs(dfmInd);
%multiInspDur=inspDur(multiInd);multiExpDur=expDur(multiInd);multiBreaths=breaths(multiInd);multiVocOnsets=vocOnsets(multiInd);multiVocOffsets=vocOffsets(multiInd); multiCorr = corrs(multiInd);
stepUpInspDur = inspDur(stepUpInd);stepUpExpDur = expDur(stepUpInd);stepUpBreaths=breaths(stepUpInd);stepUpVocOnsets=vocOnsets(stepUpInd);stepUpVocOffsets=vocOffsets(stepUpInd); stepUpCorr = corrs(stepUpInd);
compInspDur=inspDur(compInd);compExpDur=expDur(compInd);compBreaths=breaths(compInd);compVocOnsets=vocOnsets(compInd);compVocOffsets=vocOffsets(compInd); compCorr = corrs (compInd);
shortInspDur=inspDur(shortInd);shortExpDur=expDur(shortInd);shortBreaths=breaths(shortInd);shortVocOnsets=vocOnsets(shortInd);shortVocOffsets=vocOffsets(shortInd); shortCorr = corrs(shortInd);

%% breath features
%Ti
boxInspDur = nan(9,length(dfmInspDur));
boxInspDur(1,:)=dfmInspDur;boxInspDur(2,1:length(stepDownInspDur))=stepDownInspDur;boxInspDur(3,1:length(flatInspDur))=flatInspDur;
boxInspDur(4,1:length(shortInspDur))=shortInspDur;boxInspDur(5,1:length(stepUpInspDur))=stepUpInspDur;
boxInspDur(6,1:length(ufmInspDur))=ufmInspDur;boxInspDur(7,1:length(chevInspDur))=chevInspDur;
boxInspDur(8,1:length(twoStepInspDur))=twoStepInspDur;boxInspDur(9,1:length(compInspDur))=compInspDur;
boxInspDur = transpose (boxInspDur);
figure
boxplot (boxInspDur,labels,'symbol','');
ylim([0 100]);
ylabel('Ti(ms)');
%Te
boxExpDur = nan(9,length(dfmExpDur));
boxExpDur(1,:)=dfmExpDur;boxExpDur(2,1:length(stepDownExpDur))=stepDownExpDur;boxExpDur(3,1:length(flatExpDur))=flatExpDur;
boxExpDur(4,1:length(shortExpDur))=shortExpDur;boxExpDur(5,1:length(stepUpExpDur))=stepUpExpDur;
boxExpDur(6,1:length(ufmExpDur))=ufmExpDur;boxExpDur(7,1:length(chevExpDur))=chevExpDur;
boxExpDur(8,1:length(twoStepExpDur))=twoStepExpDur;boxExpDur(9,1:length(compExpDur))=compExpDur;
boxExpDur = transpose (boxExpDur);
figure
boxplot (boxExpDur,labels,'symbol','');
ylim([0 300]);
ylabel('Te(ms)');
%Corr
boxCorr = nan(9,length(dfmCorr));
boxCorr(1,:)=dfmCorr;boxCorr(2,1:length(stepDownCorr))=stepDownCorr;boxCorr(3,1:length(flatCorr))=flatCorr;
boxCorr(4,1:length(shortCorr))=shortCorr;boxCorr(5,1:length(stepUpCorr))=stepUpCorr;
boxCorr(6,1:length(ufmCorr))=ufmCorr;boxCorr(7,1:length(chevCorr))=chevCorr;
boxCorr(8,1:length(twoStepCorr))=twoStepCorr;boxCorr(9,1:length(compCorr))=compCorr;
boxCorr = transpose (boxCorr);
figure
boxplot (boxCorr,labels,'symbol','');
ylim([-1 1]);
ylabel('r');
