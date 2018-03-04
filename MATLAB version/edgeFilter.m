function cleanFROG = edgeFilter(FROG)
    edgeWidth = 2;
    firstVertical = find(FROG(:,round(end/2)),1,'first');
    lastVertical = find(FROG(:,round(end/2)),1,'last');
    firstHorizontal = find(FROG(round(end/2),:),1,'first');
    lastHorizontal = find(FROG(round(end/2),:),1,'last');
    edgeMean = mean(mean([FROG(firstVertical:firstVertical+edgeWidth, round(end/2)-edgeWidth:round(end/2)+edgeWidth)' ...
                          FROG(lastVertical-edgeWidth:lastVertical, round(end/2)-edgeWidth:round(end/2)+edgeWidth)' ...
                          FROG(firstHorizontal:firstHorizontal+edgeWidth, round(end/2)-edgeWidth:round(end/2)+edgeWidth)' ...
                          FROG(lastHorizontal-edgeWidth:lastHorizontal, round(end/2)-edgeWidth:round(end/2)+edgeWidth)']));
    cleanFROG = FROG - edgeMean;
end

