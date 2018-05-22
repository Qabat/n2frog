%   --------------------------------------------------------------
%   This function interpolates the FROG trace to new size
%   rescaling delay and lambda spacing.
%   --------------------------------------------------------------

function [interpedFROG, header] = resample(filteredFROG, header, scaleDelay, scaleLambda, N)

% read in parameters of spectrogram
lenDelay = header(1);
lenLambda = header(2);
deltaDelay = header(3);
deltaLambda = header(4);
centralLambda = header(5);

% cut the trace
cutDelay = floor(lenDelay*(1-scaleDelay)/2);
cutLambda = floor(lenLambda*(1-scaleLambda)/2);

filteredFROG(:,(end-cutDelay+1):end) = [];
filteredFROG(:,1:cutDelay) = [];
filteredFROG((end-cutLambda+1):end,:) = [];
filteredFROG(1:cutLambda,:) = [];

vectDelay = (-deltaDelay*(lenDelay-2*cutDelay)/2):deltaDelay:(deltaDelay*((lenDelay-2*cutDelay)/2-1));
vectLambda = centralLambda + ((-deltaLambda*(lenLambda-2*cutLambda)/2):deltaLambda:(deltaLambda*((lenLambda-2*cutLambda)/2-1)));

% resample axes
newdeltaDelay = deltaDelay/scaleLambda;
newdeltaLambda = deltaLambda/scaleDelay;
interDelay = (-deltaDelay*(lenDelay-2*cutDelay)/2):newdeltaDelay:(deltaDelay*((lenDelay-2*cutDelay)/2-1));
interLambda = centralLambda + ((-deltaLambda*(lenLambda-2*cutLambda)/2):newdeltaLambda:(deltaLambda*((lenLambda-2*cutLambda)/2-1)));

[XIn, YIn] = meshgrid(vectDelay, vectLambda);
[XOut, YOut] = meshgrid(interDelay, interLambda);
interpedFROG = interp2(XIn, YIn, filteredFROG, XOut, YOut,'spline');
interpedFROG(interpedFROG<0) = 0;

% stretch to NxN
newDelay = (-newdeltaDelay*(N/2)):newdeltaDelay:(newdeltaDelay*(N/2-1));
newLambda = centralLambda + ((-newdeltaLambda*(N/2)):newdeltaLambda:(newdeltaLambda*(N/2-1)));

[XIn, YIn] = meshgrid(interDelay, interLambda);
[XOut, YOut] = meshgrid(newDelay, newLambda);
interpedFROG = interp2(XIn, YIn, interpedFROG, XOut, YOut,'spline',0);
interpedFROG(interpedFROG<0) = 0;

% update header values
header(1) = N;
header(2) = N;
header(3) = abs(newDelay(2)-newDelay(1));
header(4) = abs(newLambda(2)-newLambda(1));

end