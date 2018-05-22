%   ------------------------------------------------------------
%   This function uses two simple noise filtering methods
%   applied similarly to Femtosoft FROG
%   ------------------------------------------------------------

function [cleanFROG, header] = denoise(FROG)

header = FROG(1:5,1);
FROG(1:5,:) = [];
cleanFROG = FROG;

cleanFROG = cleanFROG/max(max(cleanFROG));

cleanFROG(cleanFROG<0) = 0;

% subtract full spectrum
for ii = 1:10
    cleanFROG = cleanFROG - mean(cleanFROG(:,[1 2 end-1 end-2]),2);
    cleanFROG(cleanFROG<0) = 0;
end

% remove lone pixels
cleanFROG(bwareaopen(cleanFROG, 10)==0) = 0;
cleanFROG(cleanFROG<0) = 0;

end