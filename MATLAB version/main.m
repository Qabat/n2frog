%   ------------------------------------------------------------
%   Main program for running FROG algorithm multiple times.
%   ------------------------------------------------------------

clc;
clear;
close all; 

% general settings
fullRun = 1;
experimentalFROG = dlmread('../test data/60.txt');
fileName = 'YAG 60 test';

% set parameters of a trace
N = 128;

% 0.5, 0.6 for old, 1, 0.4 for new measurements
scaleDelay = 1;
scaleLambda = 0.4;
edgeFiltering = 0;
mirror = 'both';
flipPhase = 1; % for measuring n2 phase is flipped so the n2 sign is correct

% prepare FROG trace for running the algorithm
[experimentalFROG, header] = denoiseFROG(experimentalFROG, edgeFiltering);
[experimentalFROG, header] = resampleFROG(experimentalFROG, header, scaleDelay, scaleLambda, N);
[experimentalFROG, delays, omegas] = switchDomain(experimentalFROG, header, N);
experimentalFROG = mirrorFROG(experimentalFROG, delays, omegas, mirror);

% input parameters for FROG algorithm
errorTolerance = 1e-3;
maxIterations = 500;
whichMethod = 0; % 0 for power method, 1 for SVD method
normalRuns = 3;
bootstrapRuns = 10;
outlierTreshold = 10;

% main thing
if ~fullRun
    hidePlots = 0;
    useBootstrap = 0;
	[retrievedPulse, retrievedFROG, finalError, finalIterations] = algoFROG(experimentalFROG, errorTolerance, maxIterations, delays, omegas, flipPhase, whichMethod, hidePlots, useBootstrap);
    retrievedIntensity = abs(retrievedPulse).^2;
    retrievedPhase = angle(retrievedPulse);
    retrievedSPulse = fftshift(fft(fftshift(retrievedPulse)));
    retrievedSpectrum = abs(retrievedSPulse).^2;
    retrievedSPhase = angle(retrievedSPulse);
    outputFile = [delays' retrievedIntensity retrievedPhase 1000*omegas' retrievedSpectrum retrievedSPhase finalError*ones(length(delays),1)];
    dlmwrite('../../output/testrun.txt', outputFile, '\t');
else
	hidePlots = 1;
    useBootstrap = 0;
    for n = 1:normalRuns
        
        disp(['Current run: ' num2str(n)]);
        
        % main algorithm
        [retrievedPulse, retrievedFROG, finalError, finalIterations] = algoFROG(experimentalFROG, errorTolerance, maxIterations, delays, omegas, flipPhase, whichMethod, hidePlots, useBootstrap);

        % smooth out the pulse
        retrievedIntensity = abs(retrievedPulse).^2;
        retrievedPhase = angle(retrievedPulse);
        retrievedSPulse = fftshift(fft(fftshift(retrievedPulse)));
        retrievedSpectrum = abs(retrievedSPulse).^2;
        retrievedSPhase = angle(retrievedSPulse);
        retrievedIntensity = abs(retrievedIntensity/max(retrievedIntensity));
        retrievedSpectrum = abs(retrievedSpectrum/max(retrievedSpectrum));

        % set max intensity at at zero and zero phase at max phase
        [~, maxIntensity] = max(retrievedIntensity);
        retrievedIntensity = circshift(retrievedIntensity, length(delays)/2-maxIntensity);
        retrievedPhase = circshift(retrievedPhase, length(delays)/2-maxIntensity);
        delays = delays - delays(length(delays)/2);
        retrievedPhase = unwrap(retrievedPhase);
        width = 0.05*length(retrievedPhase);
        temp = retrievedPhase(round((length(retrievedPhase)-width)/2):round((length(retrievedPhase)+width)/2));
        retrievedPhase = retrievedPhase - max(temp);

        % save pulse to file
        outputFile = [delays' retrievedIntensity retrievedPhase 1000*omegas' retrievedSpectrum retrievedSPhase finalError*ones(length(delays),1)];
        dlmwrite(['../../output/normal/' num2str(n) '.txt'], outputFile, '\t');

    end

    % compare retrieved pulses
    intensities = [];
    phases = [];
    errors = [];
    for n=1:normalRuns
        
        file = dlmread(['../../output/normal/' num2str(n) '.txt']);

%         % plotting results of multiple FROG runs
%         figure()
%         subplot(1,2,1)
%         plot(file(:,1),file(:,2)*pi)
%         hold on
%         plot(file(:,1), file(:,3)+pi/2)
%         xlim([-1000 1000]);
%         subplot(1,2,2)
%         plot(file(:,4),file(:,5)*pi)
%         hold on
%         plot(file(:,4), file(:,6)+pi/2)
%         xlim([-50 50]);

        % collecting data for calculating errors
        intensities = [intensities file(:,2)];
        phases = [phases file(:,3)];
        errors = [errors file(1,7)];

    end
    
    % calculating pulse as a weighted average of 10 runs without bootstrap
    [minError, minIndex] = min(errors);
    intensity = intensities(:, minIndex);
    phase = phases(:, minIndex);    
    phase(intensity < 0.2) = NaN; % phase blanking
    
    useBootstrap = 1; % when using bootstrap for calculating errors make howMany = 100
    for n = 1:bootstrapRuns
        
        disp(['Current run: ' num2str(n) 'b']);
        
        % main algorithm
        [retrievedPulse, retrievedFROG, finalError, finalIterations] = algoFROG(experimentalFROG, errorTolerance, maxIterations, delays, omegas, flipPhase, whichMethod, hidePlots, useBootstrap);                
        
        % smooth out the pulse
        retrievedIntensity = abs(retrievedPulse).^2;
        retrievedPhase = angle(retrievedPulse);
        retrievedSPulse = fftshift(fft(fftshift(retrievedPulse)));
        retrievedSpectrum = abs(retrievedSPulse).^2;
        retrievedSPhase = angle(retrievedSPulse);
        retrievedIntensity = abs(retrievedIntensity/max(retrievedIntensity));
        retrievedSpectrum = abs(retrievedSpectrum/max(retrievedSpectrum));

        % set max intensity at at zero and zero phase at max phase
        [~, maxIntensity] = max(retrievedIntensity);
        retrievedIntensity = circshift(retrievedIntensity, length(delays)/2-maxIntensity);
        retrievedPhase = circshift(retrievedPhase, length(delays)/2-maxIntensity);
        delays = delays - delays(length(delays)/2);
        retrievedPhase = unwrap(retrievedPhase);
        width = 0.05*length(retrievedPhase);
        temp = retrievedPhase(round((length(retrievedPhase)-width)/2):round((length(retrievedPhase)+width)/2));
        retrievedPhase = retrievedPhase - max(temp);

        % save pulse to file
        outputFile = [delays' retrievedIntensity retrievedPhase 1000*omegas' retrievedSpectrum retrievedSPhase finalError*ones(length(delays),1)];
        dlmwrite(['../../output/bootstrap/' num2str(n) '.txt'], outputFile, '\t');

    end

    % compare retrieved pulses
    intensities = [];
    phases = [];
    figure('Position',[150 75 1600 900]);
    for n=1:bootstrapRuns

        % plotting results of multiple FROG runs
        file = dlmread(['../../output/bootstrap/' num2str(n) '.txt']);
        subplot(2,2,1)
        plot(file(:,1),file(:,2)*pi)
        hold on
        plot(file(:,1), file(:,3)+pi/2)
        xlim([-500 500]);
        ylim([-0.2 3.5]);
        title('Bootstrap in time');
        xlabel('time [fs]');
        ylabel('phase [rad]');
        subplot(2,2,2)
        plot(file(:,4),file(:,5)*pi)
        hold on
        plot(file(:,4), file(:,6)+pi/2)
        xlim([-30 30]);
        ylim([-0.2 3.5]);
        title('Bootstrap in frequency');
        xlabel('frequency [THz]');
        ylabel('phase [rad]');

        % collecting data for calculating errors
        intensities = [intensities file(:,2)];
        phases = [phases file(:,3)];

    end
    
    % plot average on top of bootstrap
    subplot(2,2,1)
    plot(delays, intensity*pi, 'r', 'LineWidth', 3);
    hold on
    plot(delays, phase+pi/2, 'r', 'LineWidth', 3);

    % remove outliers
    mask = [];
    for n=1:bootstrapRuns
        if (sqrt(mean((phases(:,n) - mean(phases, 2)).^2)) > outlierTreshold)
            mask = [mask n];
        end
    end
    phases(:,mask) = [];
    intensities(:,mask) = [];
    
    % calculating and plotting errors (bootstrap method)
    intensityError = std(intensities, 0, 2);
    phaseError = std(phases, 0, 2);
    subplot(2,2,3)
    errorbar(delays, intensity, intensityError);
    xlim([-500 500]);
    ylim([-0.2 1.2]);
    title('Intensity with errors');
    xlabel('time [fs]');
    ylabel('intensity');
    subplot(2,2,4)
    errorbar(delays, phase, phaseError);
    xlim([-500 500]);
    title('Phase with errors');
    xlabel('time [fs]');
    ylabel('phase [rad]');
    
    % display best error
    disp(['Minimum error: ' num2str(minError)]);
    
    % writing final result to file
    dlmwrite(['../../fullruns/' fileName '.txt'], [delays', intensity, phase, intensityError, phaseError], '\t');
end