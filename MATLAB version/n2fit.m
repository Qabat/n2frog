clc;
clear;
close all; 

% experimental parameters
sample = 'YAG';
P = 95;                             % average power mW
n = 1.8153;                         % linear index of refraction
initialGuess = 7;                   % starting value for fitting algorithm
reference = 'ref2';                 % which reference pulse to use
d = 2;                              % sample thickness mm
w = 254;                            % beam spot size in um
R = 1;                              % repetition rate in kHz
lambda = 1028;                      % wavelength in nm
errP = 0.03*P;                      % average power error in mW
errd = 0.02*d;                      % sample thickness error in mm
errw = 3;                           % beam spot size error in um
errlambda = 2;                      % lambda error in nm
Pf = P * (1 - ((1-n)/(1+n))^2 );    % power corrected for Fresnel reflection
blankingTreshold = 0.2;             % phase blanking treshold for fitting

% read reference pulse from the file
calibrationPulse = dlmread(['../../reference/' reference '.txt']);
calibrationIntensity = calibrationPulse(:,2);
calibrationIntensityError = calibrationPulse(:,4);
calibrationPhase = calibrationPulse(:,3);
calibrationPhaseError = calibrationPulse(:,5);

% read measured pulse from the file
pulseToFit = dlmread(['../../fullruns/' sample '/' sample ' ' num2str(P) '.txt']);
time = pulseToFit(:,1);
intensity = pulseToFit(:,2);
rawIntensity = intensity;
phase = pulseToFit(:,3);
intensityError = pulseToFit(:,4);
phaseError = pulseToFit(:,5);

% plotting inpute and reference pulse intensity
figure('Position',[450 90 1200 700]);
subplot(2,2,1)
errorbar(time, intensity, intensityError);
hold on
errorbar(time, calibrationIntensity, calibrationIntensityError);
xlim([-500 500]);
ylim([0 1.2]);
title('Intensity');
xlabel('time [fs]');
ylabel('normalized intensity');

% plotting input, reference and difference phases
subplot(2,2,2)
errorbar(time, phase, phaseError);
hold on
errorbar(time, calibrationPhase, calibrationPhaseError);

% subtract calibration phase
phase = phase - calibrationPhase;
phaseError = sqrt(phaseError.^2 + calibrationPhaseError.^2);

errorbar(time, phase, phaseError);
hold off
xlim([-500 500]);
ylim([-2 0.5]);
title('Phase');
xlabel('time [fs]');
ylabel('phase [rad]');

% calculating "almostPhase" from intensity profile, theoretical phase with n2 = 1
intensityError = intensityError/trapz(1e-15*time, intensity);
intensity = intensity/trapz(1e-15*time, intensity);
almostPhase = intensity * ((Pf*1e-3)/(R*1e3)) * (d*1e-3) * ((2/(lambda*1e-9))/(w*1e-6)^2);
almostPhaseError = intensityError * ((Pf*1e-3)/(R*1e3)) * (d*1e-3) * ((2/(lambda*1e-9))/(w*1e-6)^2);
fullAlmostPhase = almostPhase;
fullTime = time;

% blank phase at blanking treshold
time(rawIntensity < blankingTreshold) = [];
intensity(rawIntensity < blankingTreshold) = [];
intensityError(rawIntensity < blankingTreshold) = [];
almostPhase(rawIntensity < blankingTreshold) = [];
almostPhaseError(rawIntensity < blankingTreshold) = [];
phase(rawIntensity < blankingTreshold) = [];
phaseError(rawIntensity < blankingTreshold) = [];

% initially set both to max = 0
almostPhase = almostPhase - max(almostPhase);
phase = phase - max(phase);

%   ----------------------------------------------------------------------------------------------
%   FITTING
%   ----------------------------------------------------------------------------------------------

fitFun = @(b, x) b(1).*x(:,1) + b(2).*x(:,2) + b(3).*x(:,3);

almostPhase = almostPhase*1e-20;
almostPhaseError = almostPhaseError*1e-20;
fullAlmostPhase = fullAlmostPhase*1e-20;

weights = 1./phaseError.^2;
n2Model = fitnlm([almostPhase time*1e-2 ones(size(time))*1e-1], phase, fitFun, [initialGuess 1e-2 1e-2], 'Weights', weights);
n2 = n2Model.Coefficients.Estimate(1);

for ii = 1:3
    weights = 1./(phaseError.^2 + n2^2 * almostPhaseError.^2);
    n2Model = fitnlm([almostPhase time*1e-2 0.1*ones(size(time))], phase, fitFun, [initialGuess 1e-2 1e-2], 'Weights', weights);
    n2 = n2Model.Coefficients.Estimate(1);
end

n2Model = fitnlm([almostPhase time*1e-2 0.1*ones(size(time))], phase, fitFun, [initialGuess 1e-2 1e-2], 'Weights', weights);
n2 = n2Model.Coefficients.Estimate(1);
a = n2Model.Coefficients.Estimate(2);
b = n2Model.Coefficients.Estimate(3);

%   ----------------------------------------------------------------------------------------------

% experimental error
errExp = n2 * sqrt((errlambda/lambda)^2 + (errd/d)^2 + (errP/P)^2 + 4*(errw/w)^2);
% fitting error
errFit = n2Model.Coefficients.SE(1);
% total error
errn2 = sqrt(errExp^2 + errFit^2);

% adding linear and const for plotting
phase = phase - a*time*1e-2 - b/10;

% plotting fit
subplot(2,2,3);
errorbar(almostPhase, phase, phaseError, phaseError, almostPhaseError, almostPhaseError,'o')
hold on
plot(almostPhase, n2*almostPhase, 'LineWidth', 2);
title('Fitting phases to each other');
xlabel('almost phase [rad*W/cm^2]');
ylabel('phase [rad]');
subplot(2,2,4)
phase = phase + max(n2 * fullAlmostPhase);
fullAlmostPhase = n2 * fullAlmostPhase;
plot(fullTime, fullAlmostPhase, 'LineWidth', 2);
xlim([-500 500]);
hold on
plot(time, phase, 'LineWidth', 2);
xlim([-500 500]);
ylim([0 0.2 + max(fullAlmostPhase)]);
title('Fitted phases');
xlabel('time [fs]');
ylabel('phase [rad]');

% showing results on the console
disp(['sample: ' sample]);
disp(['n2: ' num2str(n2)]);
disp(['n2error: ' num2str(errn2)]);

% saving to file
time = [time; NaN([length(fullTime)-length(time), 1])];
phase = [phase; NaN([length(fullAlmostPhase)-length(phase), 1])];
dlmwrite(['../../fits/' sample ' ' num2str(P) ' ' num2str(n2) ' ' num2str(errn2) '.txt'], [fullTime, fullAlmostPhase, time, phase], '\t');