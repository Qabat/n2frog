% interpolate calibration phase for use in n2 reconstruction
% calibrationPhase = dlmread('ref pulse.txt');
% calibrationPhase = interp1(calibrationPhase(:,1), calibrationPhase(:,3), denseDelays, 'spline', 0);
% something may be wrong with interpolating calibration phase?
% using it should change all phases in 100 runs in same way
% but it seems they are more noisy this way