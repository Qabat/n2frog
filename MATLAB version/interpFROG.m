function [omegaFROG, lenDelay, lenOmega, deltaDelay, deltaOmega, centralOmega] = interpFROG(lambdaFROG)
N = 256;

c = 299.792458; % speed of light in nm/fs

% read in parameters and spectrogram
lenDelay = lambdaFROG(1);
lenLambda = lambdaFROG(2);
deltaDelay = lambdaFROG(3);
deltaLambda = lambdaFROG(4);
centralLambda = lambdaFROG(5);
centralOmega = 2*pi*c / centralLambda;
lambdaFROG(1:5,:) = [];

% change domain from lambda to omega
vectDelay = -deltaDelay*lenDelay/2:deltaDelay:deltaDelay*(lenDelay/2-1);
vectLambda = centralLambda + (-deltaLambda*lenLambda/2:deltaLambda:deltaLambda*(lenLambda/2-1));
vectOmega = 2*pi*c./vectLambda; % to get vector of delOmegas -> diff(vectOmega)
lenOmega = length(vectOmega);
omegaFROG = (lambdaFROG)./(2*pi*c) .* (vectLambda'.^2) ;

display(deltaDelay);

% interpolate for FFT-compatible NxN array
%     if (2*pi/deltaDelay < abs(vectOmega(end)-vectOmega(1)))
%         deltaOmega = 2*pi / (deltaDelay*lenDelay);
%         display(deltaOmega);
%         N = lenDelay;
%     else
%         deltaOmega = abs(vectOmega(1) - vectOmega(end)) / 255;
%         deltaDelay = 2*pi / (deltaOmega*lenOmega);
%         display(deltaOmega);
%         display(deltaDelay);
%         N = lenOmega;
%     end
%
%     deltaDelay = 10;
%     deltaOmega = 0.0001;

fExpTempSpan = max(vectDelay) - min(vectDelay);

% if 2*pi*N/fExpTempSpan < max(vectOmega) - min(vectOmega)
    vT = linspace(0, fExpTempSpan, N);
    vT = fftshift(vT);
    vT = vT - vT(1);
    vT = fftshift(vT);
    
    vAF = linspace(0, 2*pi*N./fExpTempSpan, N);
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
vectOmega = vectOmega - mean(vectOmega);

%     newVectDelay = -deltaDelay*N/2:deltaDelay:deltaDelay*(N/2-1);
%     newVectOmega = centralOmega + (deltaOmega*N/2:-deltaOmega:-deltaOmega*(N/2-1));

[XIn, YIn] = meshgrid(vectDelay,vectOmega);
%     [XOut, YOut] = meshgrid(newVectDelay,newVectOmega);
[XOut, YOut] = meshgrid(vT, vAF);

omegaFROG = interp2(XIn, YIn, omegaFROG, XOut, YOut,'spline', 0);
%     omegaFROG = interp2(vectDelay,vectOmega,omegaFROG,newVectDelay,newVectOmega,'spline');

omegaFROG = omegaFROG/max(max(omegaFROG));
end

