function [intensityFROG, electricFROG] = makeFROG(electricField)

N = length(electricField);
electricFROG = electricField*(electricField.');
	
% row rotation
for n=2:N
	electricFROG(n,:) = circshift(electricFROG(n,:), [0 1-n]);
end

% permute the columns to the right order
electricFROG=fliplr(fftshift(electricFROG,2));

% FFT each column and put 0 frequency in the correct place:
electricFROG=circshift(fft(electricFROG,[],1),ceil(N/2)-1);

% generate FROG trace (= |field|^2)
intensityFROG = abs(electricFROG).^2;
