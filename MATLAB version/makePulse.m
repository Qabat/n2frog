function outputPulse = makePulse(electricField, lastPulse, whichMethod)

N = size(electricField, 1);

%Do the exact inverse of the procedure in makeFROG...
%Undo the line: EF=circshift(fft(EF,[],1),ceil(N/2)-1);
electricField = ifft(circshift(electricField, 1-ceil(N/2)), [], 1);

%Undo the line: EF=fliplr(fftshift(EF,2));
electricField = ifftshift(fliplr(electricField),2);

%Undo the lines: for n=2:N  EF(n,:) = circshift(EF(n,:), [0 1-n]);  end
for n=2:N
	electricField(n,:) = circshift(electricField(n,:), [0, (n-1)]);
end

% Now EF is the "outer product form", see Kane1999.
% Anti-alias in time domain. See makeFROG for explanation.
electricField=electricField-tril(electricField,-ceil(N/2))-triu(electricField,ceil(N/2));

if (whichMethod == 0) % Power method
	outputPulse=electricField*(electricField'*lastPulse);
else % SVD method
	[U, S, V] = svds(electricField,1);
	outputPulse = U(:,1);
end

outputPulse=outputPulse/norm(outputPulse); %Normalize to Euclidean norm 1