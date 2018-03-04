clc;
clear;
close all;

N = 256;

% prepare FROG trace from pulse retrieved by Femtosoft FROG
Pulse = dlmread('..\testfrog\result.txt');
Time = Pulse(:,1);
Intensity = Pulse(:,2);
Phase = Pulse(:,3);

% calculate frog from pulse
% experimentalFROG = makeFROG(sqrt(Intensity).*exp(1i.*Phase));
% deltaDelay = 6.515;
% deltaOmega = 1/(N*deltaDelay);
% deltaOmega = 17; % nawet dla losowej delty zbiega, jak nie spe³nia tego warunku na FFT...

% read measured FROG from file
experimentalFROG = dlmread('..\testfrog\60.txt');
[experimentalFROG, header] = denoiseFROG(experimentalFROG);
[experimentalFROG, delays, omegas] = interpFROG(experimentalFROG, header, N);

% input parameters for FROG algorithm
errorTolerance = 1e-3;
maxIterations = 200;
whichMethod = 0;
hidePlots = 0;
useBootstrap = 0;

% main
howMany = 1;
for n=1:howMany

    [retrievedPulse, retrievedFROG, finalGError, finalIterations] = mainFROG(experimentalFROG, errorTolerance, maxIterations, delays, omegas, whichMethod, hidePlots, useBootstrap);
% 
%     retrievedIntensity = abs(retrievedPulse).^2;
%     retrievedIntensity = retrievedIntensity/max(retrievedIntensity);
%     retrievedPhase = angle(retrievedPulse);
%     retrievedPhase(retrievedIntensity<0.1) = NaN; % phase blanking
% 
%     outputFile = [retrievedIntensity, retrievedPhase];
%     method = 'power';
%     dlmwrite(['.\output_' method '\' num2str(n) '.txt'],outputFile,'\t');

end

%compare retrieved pulses
%     figure()
% for n=1:howMany
%     file = dlmread(['.\output_' method '\' num2str(i) '.txt']);
%     plot(file(:,2))
%     hold on
% end