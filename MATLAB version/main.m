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
edgeTimes = 30;
mirror = 'none'; % looks like symmetrizing only makes things worse

% read measured FROG from file
experimentalFROG = dlmread('../testfrog/60.txt');
[experimentalFROG, header] = denoiseFROG(experimentalFROG, edgeTimes);
[experimentalFROG, header] = resampleFROG(experimentalFROG, header, scaleDelay, scaleLambda, N);
[experimentalFROG, delays, omegas] = switchDomain(experimentalFROG, header, N);
experimentalFROG = mirrorFROG(experimentalFROG, mirror);

% interpolate calibration phase for use in n2 reconstruction
calibrationPhase = dlmread('ref pulse.txt');
denseDelays = linspace(delays(1), delays(end), 2^14);
calibrationPhase = interp1(calibrationPhase(:,1), calibrationPhase(:,3), denseDelays, 'spline', 0);
% something may be wrong with interpolating calibration phase?
% using it should change all phases in 100 runs in same way
% but it seems they are more noisy this way

% input parameters for FROG algorithm
errorTolerance = 1e-3;
maxIterations = 500;
whichMethod = 0;
hidePlots = 1;
useBootstrap = 1;

% main
howMany = 100;
for n=1:howMany

    % main algorithm
    [retrievedPulse, retrievedFROG, finalGError, finalIterations] = algoFROG(experimentalFROG, errorTolerance, maxIterations, delays, omegas, whichMethod, hidePlots, useBootstrap);

    % smooth out the pulse
    retrievedIntensity = abs(retrievedPulse).^2;
    retrievedIntensity = interp1(delays, retrievedIntensity, denseDelays, 'spline');
    retrievedIntensity = abs(retrievedIntensity/max(retrievedIntensity));
    retrievedPhase = angle(retrievedPulse);
    retrievedPhase = interp1(delays, retrievedPhase, denseDelays, 'spline');
%     retrievedPhase = retrievedPhase - calibrationPhase;
	retrievedPhase(retrievedIntensity < 0.1) = NaN; % phase blanking
    
    % set max intensity at  and zero phase at max phase
    [~, maxIntensity] = max(retrievedIntensity);
    retrievedIntensity = circshift(retrievedIntensity, length(denseDelays)/2-maxIntensity);
    retrievedPhase = circshift(retrievedPhase, length(denseDelays)/2-maxIntensity);
    denseDelays = denseDelays - denseDelays(length(denseDelays)/2);
    [~, maxPhase] = max(retrievedPhase);
    retrievedPhase = retrievedPhase - retrievedPhase(maxPhase);

    % save pulse to file
    outputFile = [denseDelays', retrievedIntensity', retrievedPhase'];
    method = 'power';
    dlmwrite(['.\output_' method '\' num2str(n) '.txt'],outputFile,'\t');

end

% compare retrieved pulses
close all;
figure()
for n=1:howMany
    file = dlmread(['.\output_' method '\' num2str(n) '.txt']);
    plot(file(:,1),file(:,2)*pi)
    xlim([-300 400]);
    hold on
    plot(file(:,1), file(:,3)+pi/2)
    xlim([-300 400]);
end