% which pulse to fit
reference = 'ref2 test';
sample = 'YAG 60 test';
blankingTreshold = 0.2;

% read in data from file
pulseToFit = dlmread(['../../fullruns/' sample '.txt']);
calibrationPulse = dlmread(['../../reference/' reference '.txt']);

calibrationIntensity = calibrationPulse(:,2);
calibrationPhase = calibrationPulse(:,3);
calibrationPhaseError = calibrationPulse(:,5);

time = pulseToFit(:,1);
intensity = pulseToFit(:,2);
phase = pulseToFit(:,3);
intensityError = pulseToFit(:,4);
phaseError = pulseToFit(:,5);

retrievedPulse = sqrt(intensity).*exp(-1i*phase);

hold on
plot(calibrationIntensity)
plot(calibrationPhase)
plot(intensity)
plot(phase)
hold off

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