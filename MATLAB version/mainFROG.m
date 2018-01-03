function [retrievedPulse, retrievedFROG, finalGError, finalIterations] = mainFROG(originalFROG, errorTolerance, maxIterations, deltaDelay, whichMethod, hidePlots, useBootstrap)

% RMS difference in the entries of two real matrices/vectors
rmsdiff = @(F1, F2) sqrt(mean(mean((F1-F2).^2))); 

% normalize a matrix or vector for its maximum to be 1. Must have real nonnegative entries
normalize = @(M) M/max(max(M));

% calculates alpha, the positive number that minimizes rmsdiff(Fm,alpha*Fr). See DeLong1996
alpha = @(Fm,Fr) sum(sum(Fm.*Fr))/sum(sum(Fr.^2));

% get trace dimensions
N = size(originalFROG, 1);

% normalize FROG trace to unity max intensity
originalFROG = normalize(originalFROG);

% frequency interval per pixel
deltaFreq = 1/(N*deltaDelay);

% x axis labels and plot ranges
timeLabels = (-deltaDelay*(N-1)/2:deltaDelay:deltaDelay*(N-1)/2)';
freqLabels = 1000*(-deltaFreq*(N-1)/2:deltaFreq:deltaFreq*(N-1)/2)';
timeRange = [min(timeLabels) max(timeLabels)];
freqRange = [min(freqLabels) max(freqLabels)];

% generate initial guess
initialIntensity = awgn(exp(-2*log(2)*(((0:N-1)'-N/2)/(N/10)).^2),30);
initialPhase = awgn(exp(0.1*2*pi*1i*rand(N,1)),10);
retrievedPulse = initialIntensity.*initialPhase;
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
finalIterations = 1;
finalGError = 1e10;
iterationVector = [];
errorVector = [];

while ((finalGError > errorTolerance) && (finalIterations < maxIterations))
    
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

	% use weighted average to keep peak centered at zero
	centerIndex = sum((1:N)'.*abs(retrievedPulse.^4))/sum(abs(retrievedPulse.^4));
	retrievedPulse = circshift(retrievedPulse,-round(centerIndex-N/2));
    
    % phase flip if n2 would come out negative
    if (trapz(gradient(gradient(angle(retrievedPulse(N/2-20:N/2+20)))))>0)
        retrievedPulse = abs(retrievedPulse).*exp(-1i*angle(retrievedPulse));
    end
    
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
    
    % drawing
    if (hidePlots == 0)
        figure(mainFigure)
        colormap('hot');
        
        subplot(3,2,1) % original FROG trace plot
        image([min(timeLabels) max(timeLabels)],[min(freqLabels) max(freqLabels)], sqrt(originalFROG)*64); %64 is colormap range, sqrt for electric field
        title('Original FROG trace');
        xlabel('Delay [fs]');
        ylabel('Signal frequency [THz]');
		ylim([-10 10]);
        
		subplot(3,2,3) % retrieved FROG trace plot
		image([min(timeLabels) max(timeLabels)],[min(freqLabels) max(freqLabels)], sqrt(retrievedFROG)*64);
		title(['Reconstructed FROG trace: iterations=' num2str(finalIterations) ' Femtosoft error=' num2str(finalGError)]);
		xlabel('Delay [fs]');
		ylabel('Signal frequency [THz]');
        ylim([-10 10]);
        
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
		axis([freqRange/5 0 6.5]);
        
        subplot(3,1,3) % frog error
        plot(iterationVector, errorVector );
		title('FROG error');
		xlabel('Number of iterations');
		ylim auto;
        set(gca, 'YScale', 'log')
        
		drawnow;
        finalIterations = finalIterations + 1;
        
        % pressing q ends the program
		if(strcmp(get(mainFigure,'CurrentCharacter'),'q'))
            close all;
            break;
		end
    end
end