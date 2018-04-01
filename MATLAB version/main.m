%   ------------------------------------------------------------
%   Main program for running FROG algorithm multiple times.
%   ------------------------------------------------------------

clc;
clear;
close all; 

% fullRun = 1 runs 100 times without bootstrap to calculate mean pulse and
% then 100 times with bootstrap to calculate errorbars, =0 just one time
fullRun = 1;
experimentalFROG = dlmread('../test data/100.txt');
fileName = 'CALC 100e';

% set parameters of a trace
N = 128;

% 0.5, 0.6 for old, 1, 0.4 for new measurements
scaleDelay = 0.5;
scaleLambda = 0.6;
edgeFiltering = 0;
mirror = 'none'; % symmetrize only when without it traces don't look similar
flipPhase = 1; % for measuring n2 phase is flipped so the n2 sign is correct

% prepare FROG trace for running the algorithm
[experimentalFROG, header] = denoiseFROG(experimentalFROG, edgeFiltering);
[experimentalFROG, header] = resampleFROG(experimentalFROG, header, scaleDelay, scaleLambda, N);
[experimentalFROG, delays, omegas] = switchDomain(experimentalFROG, header, N);
experimentalFROG = mirrorFROG(experimentalFROG, mirror);

% for smooth plots
denseDelays = linspace(delays(1), delays(end), 2^12);
denseOmegas = linspace(omegas(1), omegas(end), 2^12);

% input parameters for FROG algorithm
errorTolerance = 1e-3;
maxIterations = 500;
whichMethod = 0;
howMany = 100;

if ~fullRun
    hidePlots = 0;
    useBootstrap = 0;
	[retrievedPulse, retrievedFROG, finalError, finalIterations] = algoFROG(experimentalFROG, errorTolerance, maxIterations, delays, omegas, flipPhase, whichMethod, hidePlots, useBootstrap);
else
	hidePlots = 1;
    useBootstrap = 0; % when using bootstrap for calculating errors make howMany = 100
    for n = 1:howMany
        
        disp(num2str(n));
        
        % main algorithm
        [retrievedPulse, retrievedFROG, finalError, finalIterations] = algoFROG(experimentalFROG, errorTolerance, maxIterations, delays, omegas, flipPhase, whichMethod, hidePlots, useBootstrap);

        % smooth out the pulse
        retrievedIntensity = abs(retrievedPulse).^2;
        retrievedPhase = angle(retrievedPulse);
        retrievedSPulse = fftshift(fft(fftshift(retrievedPulse)));
        retrievedSpectrum = abs(retrievedSPulse).^2;
        retrievedSPhase = angle(retrievedSPulse);

        retrievedIntensity = interp1(delays, retrievedIntensity, denseDelays, 'spline');
        retrievedPhase = interp1(delays, retrievedPhase, denseDelays, 'spline');
        retrievedSpectrum = interp1(omegas, retrievedSpectrum, denseOmegas, 'spline');
        retrievedSPhase = interp1(omegas, retrievedSPhase, denseOmegas, 'spline');

        retrievedIntensity = abs(retrievedIntensity/max(retrievedIntensity));
        retrievedSpectrum = abs(retrievedSpectrum/max(retrievedSpectrum));

        % set max intensity at at zero and zero phase at max phase
        [~, maxIntensity] = max(retrievedIntensity);
        retrievedIntensity = circshift(retrievedIntensity, length(denseDelays)/2-maxIntensity);
        retrievedPhase = circshift(retrievedPhase, length(denseDelays)/2-maxIntensity);
        denseDelays = denseDelays - denseDelays(length(denseDelays)/2);
        retrievedPhase = unwrap(retrievedPhase);
        width = 0.05*length(retrievedPhase);
        temp = retrievedPhase(round((length(retrievedPhase)-width)/2):round((length(retrievedPhase)+width)/2));
        retrievedPhase = retrievedPhase - max(temp);

        % save pulse to file
        outputFile = [denseDelays' retrievedIntensity' retrievedPhase' 1000*denseOmegas' retrievedSpectrum' retrievedSPhase' finalError*ones(length(denseDelays),1)];
        dlmwrite(['../../output/' num2str(n) '.txt'], outputFile, '\t');

    end

    % compare retrieved pulses
    intensities = [];
    phases = [];
    weights = [];
    figure()
    for n=1:howMany

        % plotting results of multiple FROG runs
        file = dlmread(['../../output/' num2str(n) '.txt']);
        subplot(1,2,1)
        plot(file(:,1),file(:,2)*pi)
        hold on
        plot(file(:,1), file(:,3)+pi/2)
        xlim([-1000 1000]);
        subplot(1,2,2)
        plot(file(:,4),file(:,5)*pi)
        hold on
        plot(file(:,4), file(:,6)+pi/2)
        xlim([-50 50]);

        % collecting data for calculating errors
        intensities = [intensities file(:,2)];
        phases = [phases file(:,3)];
        weights = [weights (1/file(1,7)).^4];

    end

    % calculating pulse as a weighted average of 100 runs without bootstrap
    intensity = sum(weights.*intensities, 2)/sum(weights);
    phase = sum(weights.*phases, 2)/sum(weights);
    phase(intensity < 0.1) = NaN; % phase blanking
    figure()
    plot(denseDelays, intensity);
    hold on
    plot(denseDelays, phase);

    hidePlots = 1;
    useBootstrap = 1; % when using bootstrap for calculating errors make howMany = 100
    for n = 1:howMany
        
        disp([num2str(n) 'b']);
        
        % main algorithm
        [retrievedPulse, retrievedFROG, finalError, finalIterations] = algoFROG(experimentalFROG, errorTolerance, maxIterations, delays, omegas, flipPhase, whichMethod, hidePlots, useBootstrap);

        % smooth out the pulse
        retrievedIntensity = abs(retrievedPulse).^2;
        retrievedPhase = angle(retrievedPulse);
        retrievedSPulse = fftshift(fft(fftshift(retrievedPulse)));
        retrievedSpectrum = abs(retrievedSPulse).^2;
        retrievedSPhase = angle(retrievedSPulse);

        retrievedIntensity = interp1(delays, retrievedIntensity, denseDelays, 'spline');
        retrievedPhase = interp1(delays, retrievedPhase, denseDelays, 'spline');
        retrievedSpectrum = interp1(omegas, retrievedSpectrum, denseOmegas, 'spline');
        retrievedSPhase = interp1(omegas, retrievedSPhase, denseOmegas, 'spline');

        retrievedIntensity = abs(retrievedIntensity/max(retrievedIntensity));
        retrievedSpectrum = abs(retrievedSpectrum/max(retrievedSpectrum));

        % set max intensity at at zero and zero phase at max phase
        [~, maxIntensity] = max(retrievedIntensity);
        retrievedIntensity = circshift(retrievedIntensity, length(denseDelays)/2-maxIntensity);
        retrievedPhase = circshift(retrievedPhase, length(denseDelays)/2-maxIntensity);
        denseDelays = denseDelays - denseDelays(length(denseDelays)/2);
        retrievedPhase = unwrap(retrievedPhase);
        width = 0.05*length(retrievedPhase);
        temp = retrievedPhase(round((length(retrievedPhase)-width)/2):round((length(retrievedPhase)+width)/2));
        retrievedPhase = retrievedPhase - max(temp);

        % save pulse to file
        outputFile = [denseDelays' retrievedIntensity' retrievedPhase' 1000*denseOmegas' retrievedSpectrum' retrievedSPhase' finalError*ones(length(denseDelays),1)];
        dlmwrite(['../../output/' num2str(n) '.txt'], outputFile, '\t');

    end

    % compare retrieved pulses
    intensities = [];
    phases = [];
    weights = [];
    figure()
    for n=1:howMany

        % plotting results of multiple FROG runs
        file = dlmread(['../../output/' num2str(n) '.txt']);
        subplot(1,2,1)
        plot(file(:,1),file(:,2)*pi)
        hold on
        plot(file(:,1), file(:,3)+pi/2)
        xlim([-1000 1000]);
        subplot(1,2,2)
        plot(file(:,4),file(:,5)*pi)
        hold on
        plot(file(:,4), file(:,6)+pi/2)
        xlim([-50 50]);

        % collecting data for calculating errors
        intensities = [intensities file(:,2)];
        phases = [phases file(:,3)];
        weights = [weights (1/file(1,7)).^2];

    end
    
    % plot average on top of bootstrap
    subplot(1,2,1)
    plot(denseDelays, intensity*pi, 'LineWidth', 3);
    hold on
    plot(denseDelays, phase+pi/2, 'LineWidth', 3);
    
    % calculating and plotting errors (bootstrap method)
    intensityError = std(intensities, 0, 2);
    phaseError = std(phases, 0, 2);
    figure()
    subplot(1,2,1)
    errorbar(denseDelays, intensity, intensityError);
    xlim([-1000 1000]);
    subplot(1,2,2)
    errorbar(denseDelays, phase, phaseError);
    xlim([-1000 1000]);
    
    % writing final result to file
    dlmwrite(['../../fullruns/' fileName '.txt'], [denseDelays', intensity, phase, intensityError, phaseError], '\t');
end