function [retrievedPulse, retrievedFROG, finalGError, finalIterations] = mainFROG(originalFROG, errorTolerance, maxIterations, delays, omegas, whichMethod, hidePlots, useBootstrap)

% RMS difference in the entries of two real matrices/vectors
rmsdiff = @(F1, F2) sqrt(mean(mean((F1-F2).^2))); 

% calculates alpha, the positive number that minimizes rmsdiff(Fm,alpha*Fr). See DeLong1996
alpha = @(Fm,Fr) sum(sum(Fm.*Fr))/sum(sum(Fr.^2));

% get trace dimensions
N = size(originalFROG, 1);

% x axis labels and plot ranges
timeLabels = delays;
freqLabels = 1000*omegas;
timeRange = [min(timeLabels) max(timeLabels)];
freqRange = [min(freqLabels) max(freqLabels)];

% generate initial guess
% initialIntensity = awgn(exp(-2*log(2)*(((0:N-1)'-N/2)/(N/10)).^2),30);
initialIntensity = exp(-2*log(2)*(((0:N-1)'-N/2)/(N/10)).^2) + 2*rand(N,1);
initialPhase = 0.1*2*pi*rand(N,1);
retrievedPulse = initialIntensity.*exp(1i.*initialPhase);
[retrievedFROG, retrievedEFROG] = makeFROG(retrievedPulse);

% prepare bootstrap mask if bootstrap is to be used
if (useBootstrap == 1)
    bootstrapMask = zeros(N);
    for n = 1:(N^2)
        bootstrapMask(randi(N),randi(N)) = 1;        
    end
end

%   ------------------------------------------------------------
%   F R O G   I T E R A T I O N   A L G O R I T H M
%   ------------------------------------------------------------

mainFigure = figure('units','normalized','outerposition',[0 0 1 1]);
% mainFigure = figure('units','normalized','outerposition',[0.3 0.1 0.4 0.8]);
finalIterations = 1;
finalGError = 1e10;
testError = 1e10;
iterationVector = [];
errorVector = [];
bestError = 1;
bestFROG = [];
bestPulse = [];

stopped = 0;

while (stopped == 0)
    
    % for debugging
    %pause(1);
    
    % text output
    if hidePlots==1
		disp(['Iteration number: ' num2str(finalIterations) '  Error: ' num2str(finalGError)]);
    end
    
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
    
    % keep zero phase at zero (in power method its close to zero by itself)
    %if (whichMethod == 1)
    retrievedPulse = abs(retrievedPulse).*exp(1i*(angle(retrievedPulse) - angle(retrievedPulse(N/2))));
    %end
    
    % phase flip (and intensity flip) if n2 would come out negative
    flipRange = 3;
    if ((trapz(gradient(gradient(angle(retrievedPulse(round(N/2-flipRange):round(N/2+flipRange))))))>0) && (finalGError < 1e-2))
        retrievedPulse = flipud(abs(retrievedPulse)).*exp(-1i*flipud(angle(retrievedPulse)));
    end
    
    % add perturbation to the pulse if the error is stagnating
    if (mod(finalIterations,50) == 0)
        testError = finalGError;
    end
    if ((abs(testError - finalGError) < testError/10) && (mod(finalIterations,50) == 49))
        %change it so it doesnt use awgn but rand()
        retrievedPulse = retrievedPulse + (maxIterations - finalIterations)/10000;
    end

    	% use weighted average to keep peak centered at zero
% 	centerIndex = sum((1:N)'.*abs(retrievedPulse.^4))/sum(abs(retrievedPulse.^4));
% 	retrievedPulse = circshift(retrievedPulse, -round(centerIndex-N/2));
    
    % use maximum value to keep peak centered at zero
	[maxV, centerIndex] = max(abs(retrievedPulse));
	retrievedPulse = circshift(retrievedPulse, N/2-centerIndex);
    
    % make a FROG trace from new fields
	[retrievedFROG, retrievedEFROG] = makeFROG(retrievedPulse);
    
	% calculate FROG error G, scale Fr to best match Fm, see DeLong1996,
	% and femtosoft error - intensity weighted
	retrievedFROG = retrievedFROG*alpha(originalFROG,retrievedFROG);
	finalGError = rmsdiff(originalFROG,retrievedFROG);
    weightedError = sum(sum((retrievedFROG-originalFROG).^2 .* originalFROG))/sum(sum(originalFROG)); % check this

    % keeping track of error
    iterationVector = [iterationVector; finalIterations];
    errorVector = [errorVector; finalGError];
    
    % best pulse
    if (finalGError < bestError)
        bestError = finalGError;
        bestFROG = retrievedFROG;
        bestPulse = retrievedPulse;
    end
    
    % drawing
    if (hidePlots == 0)
        
        % show best fit after algorithm finishes
        if ~((finalGError > errorTolerance) && (finalIterations < maxIterations))
            stopped = 1;
            finalGError = bestError;
            retrievedFROG = bestFROG;
            retrievedPulse = bestPulse;
        end
            
        figure(mainFigure)
        colormap('hot');
        
        subplot(3,2,1) % original FROG trace plot
        imagesc([min(timeLabels) max(timeLabels)],[min(freqLabels) max(freqLabels)], sqrt(originalFROG));
        title('Original FROG trace');
        xlabel('Delay [fs]');
        ylabel('Signal frequency [THz]');
        pbaspect([1 1 1])
        
		subplot(3,2,3) % retrieved FROG trace plot
		imagesc([min(timeLabels) max(timeLabels)],[min(freqLabels) max(freqLabels)], sqrt(retrievedFROG));
		title('Reconstructed FROG trace');
		xlabel('Delay [fs]');
		ylabel('Signal frequency [THz]');
        pbaspect([1 1 1])
        
		subplot(3,2,2) % retrived temporal profile
		plot(timeLabels, 2*pi*abs(retrievedPulse).^2/max(abs(retrievedPulse))^2,timeLabels,angle(retrievedPulse)+pi);
		title('Reconstructed intensity and temporal phase');
		xlabel('Time [fs]');
		axis([timeRange 0 6.5]);
		
		subplot(3,2,4) % retrieved spectrum
        FFTPt=fftshift(fft(fftshift(retrievedPulse)));
		plot(freqLabels, 2*pi*abs(FFTPt).^2/max(abs(FFTPt))^2,freqLabels,angle(FFTPt)+pi);
		title('Reconstructed spectrum and spectral phase');
		xlabel('Frequency [THz]');
		axis([freqRange 0 6.5]);
        
        subplot(3,1,3) % frog error
        plot(iterationVector, errorVector );
		title(['Iterations=' num2str(finalIterations) ' Femtosoft error=' num2str(finalGError)]);
		xlabel('Number of iterations');
		ylim auto;
        set(gca, 'YScale', 'log')
        
		drawnow;

        % pressing q ends the program
		if (strcmp(get(mainFigure,'CurrentCharacter'),'q'))
            close all;
            break;
        end
        
        finalIterations = finalIterations + 1;
    end
end