from numpy import size, mean, log, pi, zeros, max, sum, round, diff, trapz, flipud, arange, exp, angle, roll, sqrt
from numpy.fft import fft, fftshift
from numpy.random import normal, uniform, randint
import numpy as np
from matplotlib.pyplot import close, figure, subplot, plot, yscale, xscale, title, grid, draw, xlabel, ylabel, imshow
from makeFROG import makeFROG
from makePulse import makePulse


def mainFROG(originalFROG, errorTolerance, maxIterations, deltaDelay, whichMethod, hidePlots, useBootstrap):

    # RMS difference in the entries of two real matrices/vectors
    rmsdiff = lambda f1, f2: (mean(mean((f1-f2)**2)))**2

    # calculates alpha, the positive number that minimizes rmsdiff(Fm,alpha*Fr). See DeLong1996
    alpha = lambda Fm, Fr: np.sum(np.sum(Fm*Fr))/np.sum(np.sum(Fr**2))

    # get trace dimensions
    N = size(originalFROG, 1)

    # normalize FROG trace to unity max intensity
    originalFROG = originalFROG/max(originalFROG)

    # frequency interval per pixel
    deltaFreq = 1/(N*deltaDelay)

    # x axis labels and plot ranges
    timeLabels = arange((-deltaDelay*(N-1)/2), (deltaDelay*((N-1)/2+1)), deltaDelay).T
    freqLabels = 1000*arange((-deltaFreq*(N-1)/2), (deltaFreq*((N-1)/2+1)), deltaFreq).T
    timeRange = [min(timeLabels), max(timeLabels)]
    freqRange = [min(freqLabels), max(freqLabels)]

    # generate initial guess
    initialIntensity = exp(-2*log(2)*((arange(0,N).T-N/2)/(N/10))**2)
    initialPhase = exp(0.1*2*pi*1j*uniform(0,1,N)) + normal(0, 10)
    retrievedPulse = initialIntensity*initialPhase
    (retrievedFROG, retrievedEFROG) = makeFROG(retrievedPulse)

    # prepare bootstrap mask if bootstrap is to be used
    if (useBootstrap == 1):
        bootstrapMask = zeros(N)
        for n in range(1,N**2):
            bootstrapMask[randint(N),randint(N)] = 1

#   ------------------------------------------------------------
#   F R O G   I T E R A T I O N   A L G O R I T H M
#   ------------------------------------------------------------

    close('all')
    mainFigure = figure(figsize=(1,1))
    finalIterations = 1
    finalGError = 1e10
    testError = 1e10
    iterationVector = []
    errorVector = []

    while ((finalGError > errorTolerance) and (finalIterations < maxIterations)):
        # for debugging
        #sleep(1)

        # text output
        if (hidePlots==1):
            print(f'Iteration number: {finalIterations} Error: {finalGError}')

        # intensity replacement step, with or without bootstrap mask
        if (useBootstrap == 1):
            temp = retrievedEFROG
            retrievedEFROG = retrievedEFROG*(sqrt(originalFROG/retrievedFROG))
            retrievedEFROG[bootstrapMask] = temp[bootstrapMask]
        else:
            retrievedEFROG = retrievedEFROG*(sqrt(originalFROG/retrievedFROG))

        # extract pulse field from FROG complex amplitude
        retrievedPulse = makePulse(retrievedEFROG, retrievedPulse, whichMethod)

        # use weighted average to keep peak centered at zero
        centerIndex = sum((arange(0, N)).T*abs(retrievedPulse**4))/sum(abs(retrievedPulse**4))
        retrievedPulse = roll(retrievedPulse, int(-round(centerIndex-N/2)))

        # keep zero phase at zero (only needed for svd frog, in power it
        # somehow stays at zero by itself)
        if (whichMethod == 1):
            retrievedPulse = abs(retrievedPulse)*exp(1j*(angle(retrievedPulse) - angle(retrievedPulse(N/2))))

        # phase flip (and intensity flip) if n2 would come out negative
        if ((trapz(diff(angle(retrievedPulse[arange(int(N/2)-25,int(N/2)+25)]),2))>0)-(finalGError < 1e-3)).all():
            retrievedPulse = flipud(abs(retrievedPulse))*exp(-1j*flipud(angle(retrievedPulse)))

        # add perturbation to the pulse if the error is stagnating
        if (finalIterations%30 == 0):
            testError = finalGError

        if ((abs(testError - finalGError) < testError/10) and (finalIterations%30 == 29)):
            retrievedPulse = retrievedPulse + normal(-1/finalIterations, 1/finalIterations)

        # make a FROG trace from new fields
        (retrievedFROG, retrievedEFROG) = makeFROG(retrievedPulse)

        # calculate FROG error G, scale Fr to best match Fm, see DeLong1996,
        # and femtosoft error - intensity weighted
        retrievedFROG = retrievedFROG*alpha(originalFROG,retrievedFROG)
        #TODO finalGError is wrong?
        finalGError = rmsdiff(originalFROG,retrievedFROG)
        weightedError = sum(sum((retrievedFROG-originalFROG)**2 * originalFROG))/sum(sum(originalFROG)) # check this

        # keeping track of error
        iterationVector.append(finalIterations)
        errorVector.append(finalGError)

        # drawing
        if (hidePlots == 0):
            #figure(mainFigure)
            #colormap('hot')
            print(f'Iteration number: {finalIterations} Error: {finalGError}')

            # original FROG trace plot
            subplot(321)
            imshow(sqrt(originalFROG)*64, extent=(min(timeLabels), max(timeLabels), min(freqLabels), max(freqLabels))) #64 is colormap range, sqrt for electric field in matlab
            title('Original FROG trace')
            xlabel('Delay [fs]')
            ylabel('Signal frequency [THz]')
            #ylim([-10 10])

            # retrieved FROG trace plot
            subplot(323)
            imshow(sqrt(retrievedFROG)*64, extent=(min(timeLabels), max(timeLabels), min(freqLabels), max(freqLabels)))
            title(f'Reconstructed FROG trace: iterations= {finalIterations}  Femtosoft error= {finalGError}')
            xlabel('Delay [fs]')
            ylabel('Signal frequency [THz]')
            #ylim([-10 10])

            # retrived temporal profile
            subplot(322)
            plot(timeLabels, 2*pi*abs(retrievedPulse)**2/max(abs(retrievedPulse))**2,timeLabels,angle(retrievedPulse)+pi)
            title('Reconstructed intensity and temporal phase')
            xlabel('Time [fs]')
            #axis([timeRange 0 6.5])

            # retrieved spectrum
            subplot(324)
            FFTPt = fftshift(fft(fftshift(retrievedPulse)))
            plot(freqLabels, 2*pi*abs(FFTPt)**2/max(abs(FFTPt))**2,freqLabels,angle(FFTPt)+pi)
            title('Reconstructed spectrum and spectral phase')
            xlabel('Frequency [THz]')
            #axis([freqRange/5 0 6.5])

            # subplot(313) # frog error
            # plot(iterationVector, errorVector )
            # title('FROG error')
            # xlabel('Number of iterations')
            # #ylim auto
            # set(gca, 'YScale', 'log')

            draw()
            finalIterations = finalIterations + 1

            return (retrievedPulse, retrievedFROG, finalGError, finalIterations)
