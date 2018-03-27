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

% Anti-alias in time domain. See makeFROG for explanation.
% electricFROG=electricFROG-tril(electricFROG,-ceil(N/2))-triu(electricFROG,ceil(N/2));

% now electricFROG is in the "outer product form", see Kane1999
if (whichMethod == 0)
    % Power method
	outputPulse = electricFROG*(electricFROG'*lastPulse);
else
    % SVD method
    % check if this part is ok
% 	[U, S, V] = svd(electricFROG,'econ');
    [U, S, V] = svds(electricFROG,3);
	outputPulse = U(:,1);
end

% normalize to Euclidean norm 1
outputPulse = outputPulse/norm(outputPulse);

end