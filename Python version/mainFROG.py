from numpy import size, mean, square, transpose, log, pi, zeros
from numpy as np
import matplotlib.pyplot as plt

def mainFROG(originalFROG, errorTolerance, maxIterations, deltaDelay, whichMethod, hidePlots, useBootstrap):

    # RMS difference in the entries of two real matrices/vectors
    rmsdiff = lambda f1, f2: sqrt(mean(mean(square(f1-f2))))

    # calculates alpha, the positive number that minimizes rmsdiff(Fm,alpha*Fr). See DeLong1996
    alpha = lambda Fm, Fr: np.sum(np.sum(multiply(Fm, Fr)))/np.sum(np.sum(square(Fr)))

    # get trace dimensions
    N = size(originalFROG, 1)

    # normalize FROG trace to unity max intensity
    originalFROG = originalFROG/amax(originalFROG)

    # frequency interval per pixel
    deltaFreq = 1/(N*deltaDelay)

    # x axis labels and plot ranges
    timeLabels = transpose(-deltaDelay*(N-1)/2:deltaDelay*(N-1)/2:deltaDelay)
    freqLabels = 1000*transpose(-deltaFreq*(N-1)/2:deltaFreq*(N-1)/2:deltaFreq)
    timeRange = [min(timeLabels) max(timeLabels)]
    freqRange = [min(freqLabels) max(freqLabels)]

    # generate initial guess
    #initialIntensity = awgn(exp(-2*log(2)*(((0:N-1)'-N/2)/(N/10)).^2),30); - version with intensity noise
    initialIntensity = exp(-2*log(2)*sqare((transpose(0:N-1)-N/2)/(N/10)))
    initialPhase = exp(0.1*2*pi*1j*np.random.uniform(0,1,N)) + np.random.normal(0, 10)
    retrievedPulse = multiply(initialIntensity, initialPhase)
    (retrievedFROG, retrievedEFROG) = makeFROG(retrievedPulse)

    # prepare bootstrap mask if bootstrap is to be used
    if (useBootstrap == 1):
        bootstrapMask = zeros(N)
        for n in range(1,N**2):
            bootstrapMask(np.random.randint(N),np.random.randint(N)) = 1

#   ------------------------------------------------------------
#   F R O G   I T E R A T I O N   A L G O R I T H M
#   ------------------------------------------------------------

    mainFigure = figure('units','normalized','outerposition',[0.05 0.1 0.9 0.8])
    finalIterations = 1
    finalGError = 1e10
    testError = 1e10
    iterationVector = []
    errorVector = []

    while ((finalGError > errorTolerance) and (finalIterations < maxIterations)):
        # for debugging
        #pause(1);

        # text output
        if (hidePlots==1):
            disp(['Iteration number: ' num2str(finalIterations) '  Error: ' num2str(finalGError)])

        # intensity replacement step, with or without bootstrap mask
        if (useBootstrap == 1):
            temp = retrievedEFROG
            retrievedEFROG = retrievedEFROG.*(sqrt(originalFROG./retrievedFROG))
            retrievedEFROG(bootstrapMask==1) = temp(bootstrapMask==1)
        else:
            retrievedEFROG = retrievedEFROG.*(sqrt(originalFROG./retrievedFROG))

	    # extract pulse field from FROG complex amplitude
	    retrievedPulse = makePulse(retrievedEFROG,retrievedPulse,whichMethod)

	    # use weighted average to keep peak centered at zero
	    centerIndex = sum((1:N)'.*abs(retrievedPulse.^4))/sum(abs(retrievedPulse.^4))
	    retrievedPulse = circshift(retrievedPulse,-round(centerIndex-N/2))

        # keep zero phase at zero (only needed for svd frog, in power it
        # somehow stays at zero by itself)
        if (whichMethod == 1):
            retrievedPulse = abs(retrievedPulse).*exp(1i*(angle(retrievedPulse) - angle(retrievedPulse(N/2))))

        # phase flip (and intensity flip) if n2 would come out negative
        if ((trapz(gradient(gradient(angle(retrievedPulse(N/2-25:N/2+25)))))>0) and (finalGError < 1e-3)):
            retrievedPulse = flipud(abs(retrievedPulse)).*exp(-1i*flipud(angle(retrievedPulse)))

        # add perturbation to the pulse if the error is stagnating
        if (mod(finalIterations,30) == 0):
            testError = finalGError

        if ((abs(testError - finalGError) < testError/10) and (mod(finalIterations,30) == 29)):
            retrievedPulse = awgn(retrievedPulse,(50+finalIterations/10))

        # make a FROG trace from new fields
	    [retrievedFROG, retrievedEFROG] = makeFROG(retrievedPulse)

	    # calculate FROG error G, scale Fr to best match Fm, see DeLong1996,
	    # and femtosoft error - intensity weighted
	    retrievedFROG = retrievedFROG*alpha(originalFROG,retrievedFROG)
	    finalGError = rmsdiff(originalFROG,retrievedFROG)
        weightedError = sum(sum((retrievedFROG-originalFROG).^2 .* originalFROG))/sum(sum(originalFROG)) # check this

        # keeping track of error
        iterationVector = [iterationVector; finalIterations]
        errorVector = [errorVector; finalGError]

        # drawing
        if (hidePlots == 0):
            figure(mainFigure)
            colormap('hot')

            subplot(3,2,1) # original FROG trace plot
            image([min(timeLabels) max(timeLabels)],[min(freqLabels) max(freqLabels)], sqrt(originalFROG)*64) #64 is colormap range, sqrt for electric field
            title('Original FROG trace')
            xlabel('Delay [fs]')
            ylabel('Signal frequency [THz]')
		    ylim([-10 10])

    		subplot(3,2,3) # retrieved FROG trace plot
		    image([min(timeLabels) max(timeLabels)],[min(freqLabels) max(freqLabels)], sqrt(retrievedFROG)*64)
		    title(['Reconstructed FROG trace: iterations=' num2str(finalIterations) ' Femtosoft error=' num2str(finalGError)])
		    xlabel('Delay [fs]')
		    ylabel('Signal frequency [THz]')
            ylim([-10 10])

		    subplot(3,2,2) # retrived temporal profile
		    plot(timeLabels, 2*pi*abs(retrievedPulse).^2/max(abs(retrievedPulse))^2,timeLabels,angle(retrievedPulse)+pi)
		    title('Reconstructed intensity and temporal phase')
		    xlabel('Time [fs]')
		    axis([timeRange 0 6.5])

		    subplot(3,2,4) # retrieved spectrum
            FFTPt=np.fft.fftshift(fft(np.fft.fftshift(retrievedPulse)))
		    plot(freqLabels, 2*pi*abs(FFTPt).^2/max(abs(FFTPt))^2,freqLabels,angle(FFTPt)+pi)
		    title('Reconstructed spectrum and spectral phase')
		    xlabel('Frequency [THz]')
		    axis([freqRange/5 0 6.5])

            subplot(3,1,3) # frog error
            plot(iterationVector, errorVector )
		    title('FROG error')
		    xlabel('Number of iterations')
		    ylim auto
            set(gca, 'YScale', 'log')

		    plt.draw()
            finalIterations = finalIterations + 1

            # pressing q ends the program
		    if (strcmp(get(mainFigure,'CurrentCharacter'),'q')):
                close all
                break

return (retrievedPulse, retrievedFROG, finalGError, finalIterations)
