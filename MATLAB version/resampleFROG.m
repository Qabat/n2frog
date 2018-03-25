function [interpedFROG, header] = resampleFROG(filteredFROG, header, scaleDelay, scaleLambda, N)

% read in parameters of spectrogram
lenDelay = header(1);
lenLambda = header(2);
deltaDelay = header(3);
deltaLambda = header(4);
centralLambda = header(5);

% calculate new vectors
vectDelay = -deltaDelay*lenDelay/2:deltaDelay:deltaDelay*(lenDelay/2-1);
vectLambda = centralLambda + (-deltaLambda*lenLambda/2:deltaLambda:deltaLambda*(lenLambda/2-1));
newDelay = linspace(-deltaDelay*lenDelay/2 * scaleDelay, deltaDelay*(lenDelay/2-1) * scaleDelay, N);
newLambda = centralLambda + linspace(-deltaLambda*lenLambda/2 * scaleLambda, deltaLambda*(lenLambda/2-1) * scaleLambda, N);

% interpolate to new axis
[XIn, YIn] = meshgrid(vectDelay, vectLambda);
[XOut, YOut] = meshgrid(newDelay, newLambda);
interpedFROG = interp2(XIn, YIn, filteredFROG, XOut, YOut,'spline', 0);
interpedFROG(interpedFROG<0) = 0;

% update header values
header(1) = N;
header(2) = N;
header(3) = abs(newDelay(2)-newDelay(1));
header(4) = abs(newLambda(2)-newLambda(1));

end

