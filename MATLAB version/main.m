%   ------------------------------------------------------------
%   Main program for running FROG algorithm multiple times.
%   ------------------------------------------------------------

clc;
clear;
close all; 

% set parameters of a trace
N = 128;
scaleDelay = 1;
scaleLambda = 0.4;

% read measured FROG from file
experimentalFROG = dlmread('../testfrog/60.txt');
[experimentalFROG, header] = denoiseFROG(experimentalFROG);
[experimentalFROG, header] = resampleFROG(experimentalFROG, header, scaleDelay, scaleLambda, N);
[experimentalFROG, delays, omegas] = switchDomain(experimentalFROG, header, N);

% input parameters for FROG algorithm
errorTolerance = 1e-3;
maxIterations = 500;
whichMethod = 0;
hidePlots = 0;
useBootstrap = 0;

% main
howMany = 1;
for n=1:howMany

    % main algorithm
    [retrievedPulse, retrievedFROG, finalGError, finalIterations] = algoFROG(experimentalFROG, errorTolerance, maxIterations, delays, omegas, whichMethod, hidePlots, useBootstrap);

    % smooth out the pulse
    retrievedIntensity = abs(retrievedPulse).^2;
    newDelays = linspace(delays(1),delays(end),1024);
    retrievedIntensity = interp1(delays, retrievedIntensity, newDelays, 'spline');
    retrievedIntensity = abs(retrievedIntensity/max(retrievedIntensity));
    retrievedPhase = angle(retrievedPulse);
    retrievedPhase = interp1(delays, retrievedPhase, newDelays, 'spline');
    retrievedPhase(retrievedIntensity<0.1) = NaN; % phase blanking
    
    % set max intensity and zero phase at zero
    [maxValue, maxIndex] = max(retrievedIntensity);
    retrievedIntensity = circshift(retrievedIntensity, -round(maxIndex-length(newDelays)/2));
    retrievedPhase = circshift(retrievedPhase, -round(maxIndex-length(newDelays)/2));
    retrievedPhase = retrievedPhase - retrievedPhase(length(newDelays)/2);
    
    % save pulse to file
    outputFile = [retrievedIntensity', retrievedPhase'];
    method = 'power';
    dlmwrite(['.\output_' method '\' num2str(n) '.txt'],outputFile,'\t');

end

% compare retrieved pulses
% close all;
% figure()
% for n=1:howMany
%     file = dlmread(['.\output_' method '\' num2str(n) '.txt']);
%     subplot(1,2,1)
%     plot(file(:,1))
%     hold on
%     subplot(1,2,2)
%     plot(file(:,2))
%     hold on
% end