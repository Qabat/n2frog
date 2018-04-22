%   --------------------------------------------------------------
%   switchDomain converts the trace
%   from wavelength domain to frequency domain
%   --------------------------------------------------------------

function [omegaFROG, newDelay, newOmega] = switchDomain(lambdaFROG, header, N)

% speed of light in nm/fs
c = 299.792458;

% read in parameters of spectrogram
lenDelay = header(1);
lenLambda = header(2);
deltaDelay = header(3);
deltaLambda = header(4);
centralLambda = header(5);

% change domain from lambda to omega
vectDelay = -deltaDelay*lenDelay/2:deltaDelay:deltaDelay*(lenDelay/2-1);
vectLambda = centralLambda + (-deltaLambda*lenLambda/2:deltaLambda:deltaLambda*(lenLambda/2-1));
vectOmega = 2*pi*c./vectLambda;
vectOmega = vectOmega - mean(vectOmega);
omegaFROG = lambdaFROG./(2*pi*c) .* (vectLambda'.^2) ;

% set new spacing satisfying FFT 
temporalSpan = max(vectDelay) - min(vectDelay);
newDelay = linspace(0, temporalSpan, N);
newDelay = fftshift(newDelay);
newDelay = fftshift(newDelay - newDelay(1));
newOmega = linspace(2*pi*N./temporalSpan, 0, N);
newOmega = fftshift(newOmega);
newOmega = fftshift(newOmega - newOmega(1));

% interpolate to new data points
[XIn, YIn] = meshgrid(vectDelay,vectOmega);
[XOut, YOut] = meshgrid(newDelay, newOmega);
omegaFROG = interp2(XIn, YIn, omegaFROG, XOut, YOut,'spline', 0);
omegaFROG(omegaFROG<0) = 0;

% shift maximum of autocorrelation to 0 delay
temporalMarginal = sum(omegaFROG, 1);
[~, maxIndex] = max(temporalMarginal);
omegaFROG = circshift(omegaFROG, [0 -abs(N/2-maxIndex)]);

% shift spectral center of mass to 0 frequency via FFT
spectralMarginal = sum(omegaFROG, 2);
cmShift = -sum(newOmega' .* abs(spectralMarginal))/sum(abs(spectralMarginal));
omegaFROG = abs(ifftshift(ifft(ifftshift(fftshift(fft(fftshift(sqrt(omegaFROG), 1), [], 1), 1).*exp(1i.*newDelay'.*cmShift), 1), [], 1), 1)).^2;

% normalize
omegaFROG = omegaFROG/max(max(omegaFROG));

end