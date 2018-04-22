clc;
clear;
close all; 

% which pulse to fit
reference = 'ref2 test';
sample = 'YAG 60 test';
initialGuess = 7.5e-20;
blankingTreshold = 0.2;

% experimental parameters
n = 1.8153;     % linear index of refraction
P = 60;         % average power mW
d = 2;          % sample thickness mm
w = 254;        % beam spot size in um
R = 1;          % repetition rate in kHz
lambda = 1028;  % wavelength in nm
errP = 1;       % average power error in mW
errd = 0.1;     % sample thickness error in mm
errw = 2;       % beam spot size error in um
errlambda = 2;  % lambda error in nm
P = P * (1 - ((1-n)/(1+n))^2 ); % power corrected for Fresnel reflection

% read reference pulse from the file
calibrationPulse = dlmread(['../../reference/' reference '.txt']);
calibrationPhase = calibrationPulse(:,3);
calibrationPhaseError = calibrationPulse(:,5);

% read measured pulse from the file
pulseToFit = dlmread(['../../fullruns/' sample '.txt']);
time = pulseToFit(:,1);
intensity = pulseToFit(:,2);
rawIntensity = intensity;
phase = pulseToFit(:,3);
intensityError = pulseToFit(:,4);
phaseError = pulseToFit(:,5);

% subtract calibration phase
phase = phase - calibrationPhase;
phaseError = sqrt(phaseError.^2 + calibrationPhaseError.^2);

% plotting input pulse with SPM phase only
figure('Position',[150 75 1600 900]);
subplot(2,2,1)
errorbar(time, intensity, intensityError);
xlim([-500 500]);
ylim([0 1.2]);
title('Intensity');
xlabel('time [fs]');
ylabel('normalized intensity');
subplot(2,2,2)
errorbar(time, phase, phaseError);
xlim([-500 500]);
ylim([-1.5 0.2]);
title('Phase');
xlabel('time [fs]');
ylabel('phase [rad]');

% calculating theoretical phase change from intensity profile with n2 = 1
intensityError = intensityError/trapz(1e-15*time, intensity);
intensity = intensity/trapz(1e-15*time, intensity);
almostPhase = intensity * ((P*1e-3)/(R*1e3)) * (d*1e-3) * ((2/(lambda*1e-9))/(w*1e-6)^2);
almostPhaseError = intensityError * ((P*1e-3)/(R*1e3)) * (d*1e-3) * ((2/(lambda*1e-9))/(w*1e-6)^2);
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

%   ----------------------------------------------------------------------------------------------
%   ITERATIVE NONLINEAR REGRESSION
%   ----------------------------------------------------------------------------------------------

weights = 1./phaseError.^2;
n2FitFun = @(b, x) b*x;
linearFitFun = @(b, x) x(:,1) + b*x(:,2);
constFitFun = @(b, x) x + b;
const = 1;

% initially set both to max = 0
almostPhase = almostPhase - max(almostPhase);
phase = phase - max(phase);

n2Model = fitnlm(almostPhase, phase, n2FitFun, initialGuess, 'Weights', weights);
n2 = n2Model.Coefficients.Estimate(1);

while (const > 0.0001)

for ii = 1:3
    weights = 1./(phaseError.^2 + n2^2 * almostPhaseError.^2);
    n2Model = fitnlm(almostPhase, phase, n2FitFun, initialGuess, 'Weights', weights);
    n2 = n2Model.Coefficients.Estimate(1);
end

constModel = fitnlm(n2*almostPhase, phase, constFitFun, 0.001, 'Weights', weights);
const = constModel.Coefficients.Estimate(1);
phase = phase - const;

linearModel = fitnlm([n2*almostPhase time], phase, linearFitFun, 0.001, 'Weights', weights);
linear = linearModel.Coefficients.Estimate(1);
phase = phase - linear * time;

constModel = fitnlm(n2*almostPhase, phase, constFitFun, 0.001, 'Weights', weights);
const = constModel.Coefficients.Estimate(1);
phase = phase - const;

end

n2Model = fitnlm(almostPhase, phase, n2FitFun, initialGuess, 'Weights', weights);
n2 = n2Model.Coefficients.Estimate(1);

%   ----------------------------------------------------------------------------------------------

% experimental error
errExp = n2 * sqrt((errlambda/lambda)^2 + (errd/d)^2 + (errP/P)^2 + 4*(errw/w)^2);
% fitting error
errFit = n2Model.Coefficients.SE(1);
% total error
errn2 = sqrt(errExp^2 + errFit^2);

% plotting fit
subplot(2,2,3);
errorbar(almostPhase, phase, phaseError, phaseError, almostPhaseError, almostPhaseError,'o')
hold on
plot(almostPhase, n2*almostPhase, 'LineWidth', 3);
title('Fitting phases to each other');
xlabel('almost phase [rad*W/cm^2]');
ylabel('phase [rad]');
subplot(2,2,4)
phase = phase + max(n2 * fullAlmostPhase);
fullAlmostPhase = n2 * fullAlmostPhase;
plot(fullTime, fullAlmostPhase);
xlim([-500 500]);
hold on
plot(time, phase);
xlim([-500 500]);
ylim([0 1]);
title('Fitted phases');
xlabel('time [fs]');
ylabel('phase [rad]');

% showing results on the console
disp(['sample: ' sample]);
disp(['n2: ' num2str(n2*1e20)]);
disp(['n2error: ' num2str(errn2*1e20)]);

% saving to file
time = [time; NaN([length(fullTime)-length(time), 1])];
phase = [phase; NaN([length(fullAlmostPhase)-length(phase), 1])];
dlmwrite(['../../fits/' sample ' ' num2str(n2*1e20) ' ' num2str(errn2*1e20) '.txt'], [fullTime, fullAlmostPhase, time, phase+max(fullAlmostPhase)], '\t');