clc;
clear;
close all; 

set(0,'DefaultAxesFontSize',20)

% experimental parameters
sample = 'FS';
P = 160;                             % average power mW
n = 1.45;                         % linear index of refraction
initialGuess = 2.5;                 % starting value for fitting algorithm
reference = 'ref1';                 % which reference pulse to use
d = 3;                              % sample thickness mm
w = 254;                            % beam spot size in um
R = 1;                              % repetition rate in kHz
lambda = 1028;                      % wavelength in nm
errP = 0.03*P;                      % average power error in mW
errd = 0.02*d;                      % sample thickness error in mm
errw = 3;                           % beam spot size error in um
errlambda = 2;                      % lambda error in nm
Pf = P * (1 - ((1-n)/(1+n))^2 );    % power corrected for Fresnel reflection
blankingTreshold = 0.15;             % phase blanking treshold for fitting
% Pf = 0.7152 * Pf;              % correction for ND filter for YVO

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
errorbar(time, intensity, intensityError, 'LineWidth', 2);
hold on
errorbar(time, calibrationIntensity, calibrationIntensityError, 'LineWidth', 2);
xlim([-500 500]);
ylim([0 1.2]);
title('Intensity comparison');
xlabel('time [fs]');
ylabel('normalized intensity');
legend('measured','reference');

% chi2
% diffIntensity = intensity - calibrationIntensity;
% T = sum((diffIntensity./intensityError).^2)
% k = length(diffIntensity)
% display(T/k);

% chi p-value
% diffIntensity = intensity - calibrationIntensity;
% diffError = intensityError.^2 + calibrationIntensityError.^2);
% T = sum((diffIntensity./diffError).^2)
% k = length(diffIntensity)
% pValue = chi2cdf(T, k-1)

% geometric
% BDMIntensity = intensity / sum(intensity);
% BDMcalibrationIntensity = calibrationIntensity / sum(calibrationIntensity);
% TBDM = sqrt(sum((intensity .* calibrationIntensity)/(sum(intensity)*sum(calibrationIntensity))))
% k = length(BDMIntensity)
% pValue = chi2cdf(TBDM, k-1)

% plotting input, reference and difference phases
subplot(2,2,2)
errorbar(time, phase, phaseError, 'LineWidth', 2);
hold on
errorbar(time, calibrationPhase, calibrationPhaseError, 'LineWidth', 2);

% subtract calibration phase
phase = phase - calibrationPhase;
fullPhase = phase;
phaseError = sqrt(phaseError.^2 + calibrationPhaseError.^2);
errorbar(time, phase, phaseError, 'LineWidth', 2);
hold off
xlim([-500 500]);
ylim([-2 0.5]);
title('Reference phase subtraction');
xlabel('time [fs]');
ylabel('phase [rad]');
legend('measured', 'reference', 'SPM-induced', 'Location', 'northeast');

% calculating "almostPhase" from intensity profile, theoretical phase with n2 = 1
intensityError = intensityError/trapz(1e-15*time, intensity);
intensity = intensity/trapz(1e-15*time, intensity);

% error propagation on intensity
almostPhase = intensity * ((Pf*1e-3)/(R*1e3)) * (d*1e-3) * ((2/(lambda*1e-9))/(w*1e-6)^2);
almostPhaseError = almostPhase .* sqrt((intensityError./intensity).^2 + (errlambda/lambda)^2 + (errd/d)^2 + (errP/P)^2 + 4*(errw/w)^2); 
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
fullPhase = fullPhase - max(phase);
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
errn2 = n2Model.Coefficients.SE(1);

%   ----------------------------------------------------------------------------------------------

% adding linear and const for plotting
phase = phase - a*time*1e-2 - b/10;
fullPhase = fullPhase - a*fullTime*1e-2 - b/10;

% plotting fit
subplot(2,2,3);
errorbar(almostPhase*1e20, phase, phaseError, phaseError, almostPhaseError*1e20, almostPhaseError*1e20, 'o', 'LineWidth', 1)
hold on
plot(almostPhase*1e20, n2*almostPhase, 'LineWidth', 3);
title('Weighted least squares regression');
xlabel(' 2\pidI(t)/\lambda [radW/cm^{2}] ');
ylabel('phase [rad]');
legend('experiment', 'best fit', 'Location', 'southeast');
% box on;
% xlim([-5 0.5]*1e19);
ax = get(gca);
ax.XAxis.Exponent = 20;

subplot(2,2,4)
fullAlmostPhase = n2 * fullAlmostPhase;
phase = phase + max(fullAlmostPhase);
fullPhase = fullPhase + max(fullAlmostPhase);
plot(fullTime, fullAlmostPhase, 'LineWidth', 3);
xlim([-500 500]);
hold on
plot(time, phase, 'LineWidth', 3);
legend('intensity', 'phase', 'Location', 'northeast');
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
dlmwrite(['../../fits/' sample ' ' num2str(0.7152*P) ' ' num2str(n2) ' ' num2str(errn2) '.txt'], [fullTime, fullAlmostPhase, fullPhase], '\t');
% print(gcf,'-dpng','-r600',['../../fits/' sample ' ' num2str(P) ' ' num2str(n2) ' ' num2str(errn2) '.png'])