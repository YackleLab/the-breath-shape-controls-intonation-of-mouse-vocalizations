function [outputArg1,outputArg2] = expRaster(delayOnsetExp,delayOffsetExp,vocInspStart,vocExpStart,highBreathTrace)
%UNTITLED Generate a heatmap raster of individual USV breaths aligned by expiration offset 
%   And PSTH underneath
% order parameters by delay between expiration and vocalization onset
[orderVocOn, indVocOn] = sort (delayOnsetExp);
orderVocOff = delayOffsetExp(indVocOn);
orderExpStarts = vocExpStart (indVocOn);
%remove USVs that are misattributed to the preceeding (onset < 0) or
%following (onset > 400 ms) breath
preVocInd = orderVocOn > 0;
orderVocOn = orderVocOn(preVocInd);
orderVocOff = orderVocOff(preVocInd);
orderExpStarts = orderExpStarts (preVocInd);
nextVocInd = orderVocOn < 200;
orderVocOn = orderVocOn(nextVocInd);
orderVocOff = orderVocOff(nextVocInd);
orderExpStarts = orderExpStarts (nextVocInd);
%pull out underlying breath traces for each vocal breath
for i = 1:length(orderExpStarts)
    tempTraceStart = (orderExpStarts(i) - 250);
    tempTraceEnd = (orderExpStarts(i) + 250);
    orderBreaths(i,:) = highBreathTrace (tempTraceStart:tempTraceEnd);
end
%plot
figure
h = imagesc (orderBreaths);
c = redblue(255);
colormap (c);
h = colorbar;
% find values for color axis
minval = min(orderBreaths,[],'all');
maxval = -minval;
caxis([minval maxval]);
plotVocOn = orderVocOn +250;
plotVocOff = orderVocOff +250;
hold on
scatter (plotVocOn, 1:length(orderVocOn),0.5,"|","MarkerEdgeColor",[0.466 0.6740 0.188]);
scatter (plotVocOff, 1:length(orderVocOn),0.5,"|","MarkerEdgeColor",[0.7 0. 0.7]);
end