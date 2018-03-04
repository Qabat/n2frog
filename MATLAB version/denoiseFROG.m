function [cleanFROG, header] = denoiseFROG(FROG)

header = FROG(1:5,:);
FROG(1:5,:) = [];

cleanFROG = FROG;
cleanFROG(cleanFROG<0) = 0;

% edge
for ii = 1:5
    cleanFROG = edgeFilter(cleanFROG);
    cleanFROG(cleanFROG<0) = 0;
end

% full spectrum
for ii = 1:5
    cleanFROG = cleanFROG - mean(cleanFROG(:,[1 2 end-1 end-2]),2);
    cleanFROG(cleanFROG<0) = 0;
end

% lone pixels
cleanFROG(bwareaopen(cleanFROG, 100)==0) = 0;
cleanFROG(cleanFROG<0) = 0;
end
