function [omegaFROG, vT, vAF] = switchDomain(lambdaFROG, header, N)

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

fExpTempSpan = max(vectDelay) - min(vectDelay);
 
% if 2*pi*N/fExpTempSpan < max(vectOmega) - min(vectOmega)
    vT = linspace(0, fExpTempSpan, N);
    vT = fftshift(vT);
    vT = vT - vT(1);
    vT = fftshift(vT);
    
%     vAF = linspace(0, 2*pi*N./fExpTempSpan, N); 
    
    vAF = linspace(2*pi*N./fExpTempSpan, 0, N);
    vAF = fftshift(vAF);
    vAF = vAF - vAF(1);
    vAF = fftshift(vAF);
% else
%     fExpAFSpan = max(vectOmega) - min(vectOmega);
%     
%     vAF = linspace(0, fExpAFSpan, N);
%     vAF = fftshift(vAF);
%     vAF = vAF - vAF(1);
%     vAF = fftshift(vAF);
%     
%     vT = linspace(0, 2*pi*N/fExpAFSpan, N);
%     vT = fftshift(vT);
%     vT = vT - vT(1);
%     vT = fftshift(vT);
% end

% interpolate to new data points
[XIn, YIn] = meshgrid(vectDelay,vectOmega);
[XOut, YOut] = meshgrid(vT, vAF);
omegaFROG = interp2(XIn, YIn, omegaFROG, XOut, YOut,'spline', 0);
omegaFROG(omegaFROG<0) = 0;

% shift maximum to 0 delay
[Maxrows, Maxcols] = find(omegaFROG == max(omegaFROG(:)));
omegaFROG = circshift(omegaFROG, [-abs(N/2-Maxrows) -abs(N/2-Maxcols)]);

% normalize spectrogram
omegaFROG = omegaFROG/max(max(omegaFROG));
end

