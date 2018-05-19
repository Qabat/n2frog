%   --------------------------------------------------------------
%   makePulse calculates the outputPulse from given electricField
%   using either Power Method or Singular Value Decomposition
%   --------------------------------------------------------------

function outputPulse = makePulse(electricFROG, lastPulse, whichMethod)

N = size(electricFROG, 1);

% undo the line: electricFROG = circshift(fft(electricFROG, [], 1),ceil(N/2)-1);
electricFROG = ifft(circshift(electricFROG, 1-ceil(N/2)), [], 1);

% undo the line: electricFROG = fliplr(fftshift(electricFROG, 2));
electricFROG = ifftshift(fliplr(electricFROG),2);

% undo row rotation
for n = 2:N
	electricFROG(n,:) = circshift(electricFROG(n,:), [0, (n-1)]);
end

% now electricFROG is in the outer product form
if (whichMethod == 0)
    % Power method
	outputPulse = electricFROG*electricFROG'*lastPulse;
else
    % SVD method
    [U, S, V] = svds(electricFROG,1);
	outputPulse = U;
end

% normalize the pulse
outputPulse = (abs(outputPulse)/max(abs(outputPulse))).*exp(1i.*angle(outputPulse));

end