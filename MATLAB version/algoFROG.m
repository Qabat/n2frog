%   ------------------------------------------------------------
%   algoFROG runs PCGPA algorithm for retrieving
%   the ultrashort pulse from measured FROG trace.
%   ------------------------------------------------------------

function [retrievedPulse, retrievedFROG, rmsError, finalIterations] = algoFROG(originalFROG, errorTolerance, maxIterations, delays, omegas, flipPhase, whichMethod, hidePlots, useBootstrap)

% get trace dimensions
N = size(originalFROG, 1);

% generate initial guess
initialIntensity = exp(-2*log(2)*(((0:N-1)'-N/2)/(N/10)).^2);
initialIntensity = initialIntensity/max(initialIntensity) + rand(N,1);
initialPhase = pi*rand(N,1);
retrievedPulse = initialIntensity.*exp(1i.*initialPhase);
[retrievedFROG, retrievedEFROG] = makeFROG(retrievedPulse);

% prepare bootstrap mask if bootstrap is to be used
if (useBootstrap == 1)
    bootstrapMask = zeros(N);
    for n = 1:(N^2)
        bootstrapMask(randi(N),randi(N)) = 1;        
    end
end

% prepare figure
mainFigure = figure('Position',[150 75 1600 900]);
colormap([1 1 1; jet(64)]);
brighten(0.4);

% create some variables
finalIterations = 1;
rmsError = 1e10;
iterationVector = [];
errorVector = [];
bestError = 1;
bestFROG = [];
bestPulse = [];
stopped = 0;

% run the main loop
while (stopped == 0)

    % intensity replacement step, with or without bootstrap mask
    if (useBootstrap == 1)
        temp = retrievedEFROG;
        retrievedEFROG = retrievedEFROG.*(sqrt(originalFROG./retrievedFROG));
        retrievedEFROG(bootstrapMask==1) = temp(bootstrapMask==1);
    else
        retrievedEFROG = retrievedEFROG.*(sqrt(originalFROG./retrievedFROG));
    end

	% extract pulse field from FROG complex amplitude
	retrievedPulse = makePulse(retrievedEFROG,retrievedPulse,whichMethod); 
    
    % set center of mass of a pulse to zero delay and zero phase at center of mass
    tau = -sum(delays' .* abs(retrievedPulse).^2)/sum(abs(retrievedPulse).^2);
    retrievedPulse = ifftshift(ifft(ifftshift(fftshift(fft(fftshift(retrievedPulse))).*exp(1i.*omegas'.*tau))));
    [~, centerPulse] = max(abs(retrievedPulse).^2);
    retrievedPulse = retrievedPulse./retrievedPulse(centerPulse);

    % phase and intensity flip if n2 would come out negative
    if ((trapz(gradient(gradient(angle(retrievedPulse((abs(retrievedPulse).^2/max(abs(retrievedPulse).^2)) > 0.3)))))>0) && (finalIterations > maxIterations/5) && (flipPhase  == 1))
        retrievedPulse = flipud(abs(retrievedPulse)).*exp(-1i*flipud(angle(retrievedPulse)));
    elseif ((trapz(gradient(gradient(angle(retrievedPulse((abs(retrievedPulse).^2/max(abs(retrievedPulse).^2)) > 0.3)))))<0) && (finalIterations > maxIterations/5) && (flipPhase  == -1))
        retrievedPulse = flipud(abs(retrievedPulse)).*exp(-1i*flipud(angle(retrievedPulse)));
    end

    % make a FROG trace from new fields
	[retrievedFROG, retrievedEFROG] = makeFROG(retrievedPulse);
    
	% normalize trace and calculate rms error
    retrievedFROG = retrievedFROG/max(max(retrievedFROG));
	rmsError = sqrt(mean(mean((originalFROG-retrievedFROG).^2)));

    % keeping track of error
    iterationVector = [iterationVector; finalIterations];
    errorVector = [errorVector; rmsError];
    
    % best pulse
    if ((rmsError < bestError) && (finalIterations > maxIterations/5))
        bestError = rmsError;
        bestFROG = retrievedFROG;
        bestPulse = retrievedPulse;
    end
    
    % show best fit after algorithm finishes
    if ~((rmsError > errorTolerance) && (finalIterations < maxIterations))
        stopped = 1;
        rmsError = bestError;
        retrievedFROG = bestFROG;
        retrievedPulse = bestPulse;
    end
           
    % text output when plots are turned off
    if (hidePlots == 1)
        if stopped == 1
            disp(['Best error: ' num2str(rmsError)]);
        end
        close all;
    end
    
    % draw the plots
    if (hidePlots == 0)

        figure(mainFigure)
        
        % original FROG trace plot
        subplot(2,3,1)
        h = pcolor(delays, 1000*omegas, sqrt(originalFROG));
        set(h, 'EdgeColor', 'none');
        caxis([0.01 1])
        title('Original FROG trace');
        xlabel('Delay [fs]');
        ylabel('Signal frequency [THz]');
        pbaspect([1 1 1])
        
        % retrieved FROG trace plot
		subplot(2,3,2) 
        h = pcolor(delays, 1000*omegas, sqrt(retrievedFROG));
        set(h, 'EdgeColor', 'none');
        caxis([0.01 1])
		title('Reconstructed FROG trace');
		xlabel('Delay [fs]');
		ylabel('Signal frequency [THz]');
        pbaspect([1 1 1])
        
        % difference between original and retrieved
        subplot(2,3,3)
        h = pcolor(delays, 1000*omegas, sqrt(originalFROG)-sqrt(retrievedFROG));
        set(h, 'EdgeColor', 'none');
        caxis([0.01 0.1])
		title('Difference between FROG traces');
		xlabel('Delay [fs]');
		ylabel('Signal frequency [THz]');
        pbaspect([1 1 1])
        
        % retrieved temporal profile
		subplot(2,3,4)
		plot(delays, 2*pi*abs(retrievedPulse).^2/max(abs(retrievedPulse))^2, delays, angle(retrievedPulse)+pi, 'LineWidth', 2);
		title('Reconstructed intensity and temporal phase');
		xlabel('Time [fs]');
		axis([[min(delays) max(delays)] 0 6.5]);
        xlim([-500 500]);

        % retrieved spectrum
		subplot(2,3,5)
        FFTPt=fftshift(fft(fftshift(retrievedPulse)));
		plot(1000*omegas, 2*pi*abs(FFTPt).^2/max(abs(FFTPt))^2, 1000*omegas, angle(FFTPt)+pi, 'LineWidth', 2);
		title('Reconstructed spectrum and spectral phase');
		xlabel('Frequency [THz]');
		axis([[min(1000*omegas) max(1000*omegas)] 0 6.5]);
        xlim([-40 40]);
        
        % frog error
        subplot(2,3,6)
        plot(iterationVector, errorVector, 'LineWidth', 2);
		title(['Iteration ' num2str(finalIterations) '; Error = ' num2str(rmsError)]);
		xlabel('Number of iterations');
		ylim auto;
        set(gca, 'YScale', 'log')
        
		drawnow;
        % print(gcf,'foo.png','-dpng','-r600');         
        % pressing q ends the program
		if (strcmp(get(mainFigure,'CurrentCharacter'),'q'))
            close all;
            break;
        end
        
    end
    
    finalIterations = finalIterations + 1;
    
end