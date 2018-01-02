function [retrievedPulse, retrievedTrace, finalGError, finalIterations] = mainFROG(originalTrace, errorTolerance, maxIterations, hidePlots, deltaDelay, bootStrap)

% define simple functions
rmsdiff = @(F1, F2) sqrt(mean(mean((F1-F2).^2))); %RMS difference in the entries of two real matrices/vectors.
normalize = @(M) M/max(max(M)); %normalize a matrix or vector for its maximum to be 1. Must have real nonnegative entries.
calcalpha = @(Fm,Fr) sum(sum(Fm.*Fr))/sum(sum(Fr.^2)); %calculates alpha, the positive number that minimizes rmsdiff(Fm,alpha*Fr). See DeLong1996

%   Get trace dimensions
N = size(originalTrace, 1);

% Default behavior for various inputs...
if (~exist('maxIterations', 'var')||isempty(maxIterations))
    maxIterations = inf;
end
if (~exist('errorTolerance', 'var')||isempty(errorTolerance))
    errorTolerance = 0;
end
if (~exist('deltaDelay','var')||isempty(deltaDelay))
	deltaDelay = 1;
end
if (~exist('bootStrap','var')||isempty(bootStrap))
	bootStrap = 0;
end

% frequency interval per pixel
deltaFreq = 1/(N*deltaDelay);

% x axis labels
timeLabels = (-deltaDelay*(N-1)/2:deltaDelay:deltaDelay*(N-1)/2)';
freqLabels = 1000*(-deltaFreq*(N-1)/2:deltaFreq:deltaFreq*(N-1)/2)';

% plot range to show
timeRange = [min(timeLabels) max(timeLabels)];
freqRange = [min(freqLabels) max(freqLabels)];

% generate initial guess
initialIntensity = awgn(exp(-2*log(2)*(((0:N-1)'-N/2)/(N/10)).^2),30);
initialPhase = awgn(exp(0.1*2*pi*1i*rand(N,1)),10);
retrievedPulse = initialIntensity.*initialPhase;

% normalize FROG trace to unity max intensity
originalTrace = normalize(originalTrace);

% EFr is reconstructed FROG trace complex amplitudes ( Fr=|EFr|.^2 )
[retrievedTrace, EFr] = makeFROG(retrievedPulse);

% prepare bootstrap mask if bootstrap is to be used
if (bootStrap == 1)
    bootstrapMask = zeros(N);
    for iter = 1:(N^2)
        bootstrapMask(randi(N),randi(N)) = 1;        
    end
end

%   ------------------------------------------------------------
%   F R O G   I T E R A T I O N   A L G O R I T H M
%   ------------------------------------------------------------
mainFigure = figure('units','normalized','outerposition',[0 0 1 1]);
finalIterations = 0;
finalGError = 1e10;
while ((finalGError > errorTolerance) && (finalIterations < maxIterations))
	finalIterations = finalIterations + 1;
	
    % text output
    if hidePlots==1
		disp(['Iteration number: ' num2str(finalIterations) '  Error: ' num2str(finalGError)]);
    end
    
    % intensity replacement step, with or without bootstrap mask
    if (bootStrap == 1)
        oldEFr = EFr;
        EFr = EFr.*(sqrt(originalTrace./retrievedTrace));
        EFr(bootstrapMask==1) = oldEFr(bootstrapMask==1);
    else
        EFr = EFr.*(sqrt(originalTrace./retrievedTrace));
    end
    
	% extract pulse field from FROG complex amplitude
	retrievedPulse = makePulse(EFr,retrievedPulse,0); 

	% use weighted average to keep peak centered at zero
	centerIndex = sum((1:N)'.*abs(retrievedPulse.^4))/sum(abs(retrievedPulse.^4));
	retrievedPulse = circshift(retrievedPulse,-round(centerIndex-N/2));
	
    % make a FROG trace from new fields
	[retrievedTrace, EFr] = makeFROG(retrievedPulse);
    
	% calculate FROG error G, scale Fr to best match Fm, see DeLong1996
	retrievedTrace = retrievedTrace*calcalpha(originalTrace,retrievedTrace);
	finalGError = rmsdiff(originalTrace,retrievedTrace);
    % intensity weighted error (femtosoft FROG - like?)
    weightedError = sum(sum((retrievedTrace-originalTrace).^2 .* originalTrace))/sum(sum(originalTrace));

    if hidePlots==0
        figure(mainFigure)
        colormap('hot');
        
        subplot(3,2,1) % original FROG trace plot
        image([min(timeLabels) max(timeLabels)],[min(freqLabels) max(freqLabels)], sqrt(originalTrace)*64); %64 is colormap range, sqrt for electric field
        title('Original FROG trace');
        xlabel('Delay [fs]');
        ylabel('Signal frequency [THz]');

		subplot(3,2,2) % retrieved FROG trace plot
		image([min(timeLabels) max(timeLabels)],[min(freqLabels) max(freqLabels)], sqrt(retrievedTrace)*64);
		title(['Reconstructed FROG trace: iterations=' num2str(finalIterations) ' Femtosoft error=' num2str(weightedError)]);
		xlabel('Delay [fs]');
		ylabel('Signal frequency [THz]');

		subplot(3,1,2) % retrived temporal profile
		plot(timeLabels, 2*pi*abs(retrievedPulse).^2/max(abs(retrievedPulse))^2,timeLabels,angle(retrievedPulse)+pi);
		title('Reconstructed intensity and temporal phase');
		xlabel('Time [fs]');
		axis([timeRange 0 6.5]);
		
		subplot(3,1,3) % retrieved spectrum
        FFTPt=fftshift(fft(fftshift(retrievedPulse)));
		plot(freqLabels, 2*pi*abs(FFTPt).^2/max(abs(FFTPt))^2,freqLabels,angle(FFTPt)+pi);
		title('Reconstructed spectrum and spectral phase');
		xlabel('Frequency [THz]');
		axis([freqRange 0 6.5]);
		
		drawnow;
    end
end
%   ------------------------------------------------------------
%   E N D   O F   A L G O R I T H M
%   ------------------------------------------------------------