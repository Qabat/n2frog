from numpy import sqrt, exp, multiply, square amax, angle, savetxt
import matplotlib.pyplot as plt


# prepare FROG trace from pulse retrieved by Femtosoft FROG
Pulse = load('..\testfrog\result.txt')
Time = Pulse[:,0]
Intensity = Pulse[:,1]
Phase = Pulse[:,2]
(computedFROG, electricFROG) = makeFROG(multiply(sqrt(Intensity), exp(mutiply(1j, Phase))))

# input parameters for FROG algorithm
errorTolerance = 8e-6
maxIterations = 500
deltaDelay = 6.515
whichMethod = 0
hidePlots = 0
useBootstrap = 0

# main
howMany = 1
for n in range(1,howMany):

    (retrievedPulse, retrievedFROG, finalGError, finalIterations) = mainFROG(computedFROG, errorTolerance, maxIterations, deltaDelay, whichMethod, hidePlots, useBootstrap)

    retrievedIntensity = square(abs(retrievedPulse))
    retrievedIntensity = retrievedIntensity/amax(retrievedIntensity)
    retrievedPhase = angle(retrievedPulse)
    retrievedPhase(retrievedIntensity < 0.1) = 0 # phase blanking

    outputFile = (Time, retrievedIntensity, retrievedPhase)
    method = 'power'
    savetxt(('.\output_' method '\\' num2str(i) '.txt'),outputFile,'\t')

# compare retrieved pulses
# for n in range(1,howMany):
#     file = load(['.\output_' method '\' num2str(i) '.txt']);
#     figure()
#     plot(file(:,2))
#     hold on
