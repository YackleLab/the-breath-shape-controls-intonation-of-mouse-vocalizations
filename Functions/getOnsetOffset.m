function [onHist,offHist,orderOn,orderOff,orderBreaths,normOnHist,normOffHist,normOnsets,normOffsets] = getOnsetOffset(breaths,vocOnset,vocOffset,edges,ind,expDur,normEdges)
%Compile vocal timing stats of a subset of data specified by indices and
%compute histograms
%index breaths, onset and offset
indBreaths = breaths(ind,:);
onsets = vocOnset(ind);
offsets = vocOffset(ind);
expLength = expDur(ind);
%trim out breaths misattributed to the previous breath (onset < 0)
preVocInd = onsets>0;
onsets = onsets(preVocInd);
offsets = offsets(preVocInd);
indBreaths = indBreaths(preVocInd,:);
expLength = expLength(preVocInd,:);
% calculate ons and offs normalized to exp length
normOnsets = onsets./expLength;
normOffsets = offsets./expLength;
% calculate histogram
[onHist ~] = histcounts(onsets,edges);
[offHist ~] = histcounts(offsets,edges);
% calculate normalized histogram
[normOnHist ~] = histcounts(normOnsets,normEdges);
[normOffHist ~] = histcounts(normOffsets, normEdges);
normOnHist = normOnHist/length(onsets);
normOffHist = normOffHist/length(onsets);
% order breaths onsets and offsets for raster
[orderOn orderInd] = sort(onsets);
orderOff = offsets(orderInd);
orderBreaths = indBreaths(orderInd,:);
end