%   ------------------------------------------------------------
%   denoiseFROG runs three noise cleaning algorithms
%   that remove excess noise while preserving
%   important parts of the trace
%   ------------------------------------------------------------

function [cleanFROG, header] = denoiseFROG(FROG, edgeFiltering)

header = FROG(1:5,1);
FROG(1:5,:) = [];
cleanFROG = FROG;

cleanFROG = cleanFROG/max(max(cleanFROG));

cleanFROG(cleanFROG<0) = 0;

% edge
for ii = 1:edgeFiltering
    cleanFROG = edgeFilter(cleanFROG);
    cleanFROG(cleanFROG<0) = 0;
end

% full spectrum
for ii = 1:10
    cleanFROG = cleanFROG - mean(cleanFROG(:,[1 2 end-1 end-2]),2);
    cleanFROG(cleanFROG<0) = 0;
end

% lone pixels
cleanFROG(bwareaopen(cleanFROG, 10)==0) = 0;
cleanFROG(cleanFROG<0) = 0;

end