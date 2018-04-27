%   ------------------------------------------------------------
%   Main program for running FROG algorithm multiple times.
%   ------------------------------------------------------------

clc;
clear;
close all; 

% general settings
fullRun = 1;
sample = 'YAG';
power = '95';
experimentalFROG = dlmread(['../../measurements/' sample '/' power '/' power '.txt']);

% set parameters of a trace
N = 128;

% 0.5, 0.7 for old, 1, 0.4 for new measurements
scaleDelay = 1;
scaleLambda = 0.4;
mirror = 'both';

% for measuring n2 phase is flipped so the n2 sign is correct
% 0 no flipping, 1 flip for positive n2, -1 flip for negative n2
flipPhase = 1;

% prepare FROG trace for running the algorithm
[experimentalFROG, header] = denoiseFROG(experimentalFROG);
[experimentalFROG, header] = resampleFROG(experimentalFROG, header, scaleDelay, scaleLambda, N);
[experimentalFROG, delays, omegas] = switchDomain(experimentalFROG, header, N);
experimentalFROG = mirrorFROG(experimentalFROG, delays, omegas, mirror);

% input parameters for FROG algorithm
errorTolerance = 1e-3;
maxIterations = 500;
whichMethod = 0; % 0 for power method, 1 for SVD method
normalRuns = 30;
bootstrapRuns = 100;
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
else % normal runs
	hidePlots = 1;
    useBootstrap = 0;
    intensities = [];
    phases = [];
    errors = [];
    for n = 1:normalRuns
        
        disp(['Current run: ' num2str(n)]);
        
        % main algorithm
        [retrievedPulse, retrievedFROG, finalError, finalIterations] = algoFROG(experimentalFROG, errorTolerance, maxIterations, delays, omegas, flipPhase, whichMethod, hidePlots, useBootstrap);

        % save pulse to file
        retrievedIntensity = abs(retrievedPulse).^2;
        retrievedPhase = angle(retrievedPulse);
        retrievedSPulse = fftshift(fft(fftshift(retrievedPulse)));
        retrievedSpectrum = abs(retrievedSPulse).^2;
        retrievedSPhase = angle(retrievedSPulse);
        retrievedIntensity = abs(retrievedIntensity/max(retrievedIntensity));
        retrievedSpectrum = abs(retrievedSpectrum/max(retrievedSpectrum));      
        outputFile = [delays' retrievedIntensity retrievedPhase 1000*omegas' retrievedSpectrum retrievedSPhase finalError*ones(length(delays),1)];
        dlmwrite(['../../output/normal/' num2str(n) '.txt'], outputFile, '\t');

        intensities = [intensities retrievedIntensity];
        phases = [phases retrievedPhase];
        errors = [errors finalError];
    end
    
    % picking the pulse with lowest error
    [minError, minIndex] = min(errors);
    intensity = intensities(:, minIndex);
    phase = phases(:, minIndex);    

    % bootstrap runs
    useBootstrap = 1;
    for n = 1:bootstrapRuns
        
        disp(['Current run: ' num2str(n) 'b']);
        
        % main algorithm
        [retrievedPulse, retrievedFROG, finalError, finalIterations] = algoFROG(experimentalFROG, errorTolerance, maxIterations, delays, omegas, flipPhase, whichMethod, hidePlots, useBootstrap);                

        % maximizing temporal pulse overlap
        bestOverlap = 1000;
        bestTau = 0;
        for tau = -200:200
            
            overlapPulse = ifftshift(ifft(ifftshift(fftshift(fft(fftshift(retrievedPulse))).*exp(-1i.*omegas'.*tau))));
            [~, centerPulse] = max(abs(overlapPulse).^2);
            overlapPulse = overlapPulse./overlapPulse(centerPulse);
            shiftedIntensity = abs(overlapPulse).^2;
    
            shiftedOverlap = trapz(abs(shiftedIntensity-intensity));
            
            if (shiftedOverlap < bestOverlap)
                bestOverlap = shiftedOverlap;
                bestTau = tau;
            end
            
        end
        
        retrievedPulse = ifftshift(ifft(ifftshift(fftshift(fft(fftshift(retrievedPulse))).*exp(-1i.*omegas'.*bestTau))));

        % save pulse to file
        retrievedIntensity = abs(retrievedPulse).^2;
        retrievedPhase = angle(retrievedPulse);
        retrievedSPulse = fftshift(fft(fftshift(retrievedPulse)));
        retrievedSpectrum = abs(retrievedSPulse).^2;
        retrievedSPhase = angle(retrievedSPulse);
        retrievedIntensity = abs(retrievedIntensity/max(retrievedIntensity));
        retrievedSpectrum = abs(retrievedSpectrum/max(retrievedSpectrum));

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
        
        bootstrapPhase = file(:,3);
        bootstrapPhase(intensity < 0.1) = [];
        mainPhase = phase;
        mainPhase(intensity < 0.1) = [];
        
        if (mean((bootstrapPhase - mainPhase).^4) < outlierTreshold)
            
            % maximizing phase overlap
            bestOverlap = 1000;
            bestConst = 0;
            for const = -1:0.005:1
    
                shiftedPhase = bootstrapPhase + const;
    
                shiftedOverlap = trapz(abs(shiftedPhase - mainPhase));
    
                if (shiftedOverlap < bestOverlap)
                    bestOverlap = shiftedOverlap;
                    bestConst = const;
                end
    
            end
    
            subplot(2,2,1)
            plot(file(:,1), file(:,2)*pi)
            hold on
            plot(file(:,1), file(:,3) + bestConst + pi/2)
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
            phases = [phases file(:,3) + bestConst];
        end
        
    end
    
    % calculating bootstrap errors
    intensityError = std(intensities, 0, 2);
    phaseError = std(phases, 0, 2);
    
    % plots
    subplot(2,2,1)
    plot(delays, intensity*pi, 'r-*', 'LineWidth', 1);
    hold on
    plot(delays, phase+pi/2, 'r-*', 'LineWidth', 1);
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
    
    % display best error and write to file
    disp(['Minimum error: ' num2str(minError)]);
    disp(['Bootstrap runs: ' num2str(size(phases,2))]);
    dlmwrite(['../../fullruns/' sample '/' sample ' ' power '.txt'], [delays', intensity, phase, intensityError, phaseError], '\t');
end