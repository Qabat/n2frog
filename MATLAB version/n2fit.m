clc;
clear;
close all; 

% which pulse to fit
reference = 'ref2';
sample = 'YAG 60';

% experimental parameters
n = 1.8153;     % linear index of refraction
P = 60;        % average power mW
d = 2;          % sample thickness mm
w = 254;        % beam spot size in um
R = 1;          % repetition rate in kHz
lambda = 1028;  % wavelength in nm
errP = 1;       % average power error in mW
errd = 0.1;     % sample thickness error in mm
errw = 2;       % beam spot size error in um
errlambda = 2;  % lambda error in nm
P = P * (1 - ((1-n)/(1+n))^2 ); % power corrected for Fresnel reflection

% read in data from file
pulseToFit = dlmread(['../../fullruns/' sample '.txt']);
calibrationPulse = dlmread(['../../reference/' reference '.txt']);
calibrationPhase = calibrationPulse(:,3);
calibrationPhaseError = calibrationPulse(:,5);
time = pulseToFit(:,1);
intensity = pulseToFit(:,2);
phase = pulseToFit(:,3);
intensityError = pulseToFit(:,4);
phaseError = pulseToFit(:,5);

% subtract calibration phase
phase = phase - calibrationPhase;
phase = phase - max(phase);
phaseError = sqrt(phaseError.^2 + calibrationPhaseError.^2);

% plotting input pulse
figure('Position',[150 75 1600 900]);
subplot(2,2,1)
errorbar(time, intensity, intensityError);
xlim([-1000 1000]);
ylim([0 1.2]);
title('Intensity');
xlabel('time [fs]');
ylabel('normalized intensity');
subplot(2,2,3)
errorbar(time, phase, phaseError);
title('Phase');
xlabel('time [fs]');
ylabel('phase [rad]');

% calculating theoretical phase change from intensity profile with n2 = 1
intensity = intensity/trapz(1e-15*time, intensity);
tPhase = intensity * ((P*1e-3)/(R*1e3)) * (d*1e-3) * ((2/(lambda*1e-9))/(w*1e-6)^2);
fulltPhase = tPhase;
fullTime = time;

% cut data to the size of phase vector
time(isnan(phase)) = [];
tPhase(isnan(phase)) = [];
tPhase = tPhase - max(tPhase);
intensityError(isnan(phase)) = [];
phaseError(isnan(phase)) = [];
phase(isnan(phase)) = [];

% fitting curves without errors for now
n2 = lsqcurvefit(@(n2, tPhase) n2 * tPhase, 3e-20, tPhase, phase);
% b³¹d eskperymentalny
errExp = n2 * sqrt((errlambda/lambda)^2 + (errd/d)^2 + (errP/P)^2 + 4*(errw/w)^2);
% b³¹d fitowania
errFit = 0;
% b³¹d ca³kowity
errn2 = sqrt(errExp^2 + errFit^2);

% b³¹d w labview liczony by³ jako iloczyn kowariancji i residuum
% dodatkowo by³ liczony b³¹d eksperymentalny
% tu jako kolejne Ÿród³o b³êdu trzeba uwzglêdniæ b³edy z bootstrapa

% plotting fit
subplot(1,2,2)
fulltPhase = n2 * fulltPhase;
plot(fullTime, fulltPhase);
hold on
plot(time, phase + max(fulltPhase));
xlim([-1000 1000]);
ylim([0 max(fulltPhase) + 0.1]);
title('Phase fitting');
xlabel('time [fs]');
ylabel('phase [rad]');

display(n2*1e20)
display(errn2*1e20)

temp = (length(fulltPhase)-length(phase))/2;
nany = NaN(temp);
phase = [nany(1,:)  phase'  nany(1,:)]';
dlmwrite(['../../fits/' sample ' ' num2str(n2*1e20) ' ' num2str(errn2*1e20) '.txt'], [fullTime, fulltPhase, phase+max(fulltPhase)], '\t');