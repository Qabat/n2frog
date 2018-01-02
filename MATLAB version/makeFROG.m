function [IF, electricField] = makeFROG(Pulse)

N = length(Pulse);
electricField = Pulse*(Pulse.');
	
% row rotation
for n=2:N
	electricField(n,:) = circshift(electricField(n,:), [0 1-n]);
end

% permute the columns to the right order
electricField=fliplr(fftshift(electricField,2));

% FFT each column and put 0 frequency in the correct place:
electricField=circshift(fft(electricField,[],1),ceil(N/2)-1);

% generate FROG trace (= |field|^2)
IF = abs(electricField).^2;
