clc;
clear;
close all;

% prepare FROG trace from pulse retrieved by Femtosoft FROG
Pulse = dlmread('..\testfrog\result.txt');
Time = Pulse(:,1);
Intensity = Pulse(:,2);
Phase = Pulse(:,3);
computedFROG = makeFROG(sqrt(Intensity).*exp(1i.*Phase));

% input parameters for FROG algorithm
errorTolerance = 0.0000001;
maxIterations = 100;
deltaDelay = 6.515;
whichMethod = 0;
hidePlots = 0;
useBootstrap = 0;

% main
howMany = 1;
for n=1:howMany
   
    [retrievedPulse, retrievedFROG, finalGError, finalIterations] = mainFROG(computedFROG, errorTolerance, maxIterations, deltaDelay, whichMethod, hidePlots, useBootstrap);

    retrievedIntensity = abs(retrievedPulse).^2;
    retrievedIntensity = retrievedIntensity/max(retrievedIntensity);
    retrievedPhase = angle(retrievedPulse);
    retrievedPhase(retrievedIntensity<0.1) = 0; % phase blanking

    outputFile = [Time, retrievedIntensity, retrievedPhase];
    method = 'power';
    dlmwrite(['.\output_' method '\' num2str(i) '.txt'],outputFile,'\t');

end

% compare retrieved pulses
% for n=1:howMany
%     file = dlmread(['.\output_' method '\' num2str(i) '.txt']);
%     figure()
%     plot(file(:,2))
%     hold on
% end