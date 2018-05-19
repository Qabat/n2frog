%   ------------------------------------------------------------
%   makeFROG calculates FROG trace from given electricField
%   ------------------------------------------------------------

function [intensityFROG, electricFROG] = makeFROG(inputPulse)

N = length(inputPulse);
electricFROG = inputPulse*(inputPulse.');

% row rotation
for n=2:N
	electricFROG(n,:) = circshift(electricFROG(n,:), [0 1-n]);
end

% permute the columns to the right order
electricFROG = fliplr(fftshift(electricFROG, 2));

% FFT each column and put 0 frequency in the correct place:
electricFROG = circshift(fft(electricFROG, [], 1),ceil(N/2)-1);

% generate FROG trace
intensityFROG = abs(electricFROG).^2;

end