function outputPulse = makePulse(electricFROG, lastPulse, whichMethod)

N = size(electricFROG, 1);

%Do the exact inverse of the procedure in makeFROG...
%Undo the line: EF=circshift(fft(EF,[],1),ceil(N/2)-1);
electricFROG = ifft(circshift(electricFROG, 1-ceil(N/2)), [], 1);

%Undo the line: EF=fliplr(fftshift(EF,2));
electricFROG = ifftshift(fliplr(electricFROG),2);

%Undo the lines: for n=2:N  EF(n,:) = circshift(EF(n,:), [0 1-n]);  end
for n=2:N
	electricFROG(n,:) = circshift(electricFROG(n,:), [0, (n-1)]);
end

% Now EF is the "outer product form", see Kane1999.
% Anti-alias in time domain. See makeFROG for explanation.
electricFROG=electricFROG-tril(electricFROG,-ceil(N/2))-triu(electricFROG,ceil(N/2));

if (whichMethod == 0) % Power method
	outputPulse=electricFROG*(electricFROG'*lastPulse);
else % SVD method
	[U, S, V] = svds(electricFROG,1);
	outputPulse = U(:,1);
end

outputPulse=outputPulse/norm(outputPulse); %Normalize to Euclidean norm 1