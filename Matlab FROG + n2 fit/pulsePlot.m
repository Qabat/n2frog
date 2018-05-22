close all;
warning('off','all');

pulse = dlmread('YVO 90.txt');
ref = dlmread('ref2.txt');

time = pulse(:,1);
intensity = pulse(:,2);
phase = unwrap(pulse(:,3))-18.85;

intensityRef = ref(:,2);
phaseRef = unwrap(ref(:,3)-12.665);

% interpolate
denseTime = linspace(time(1), time(end), length(time)*1)';
intensity = interp1(time, intensity, denseTime, 'Spline');
phase = interp1(time, phase, denseTime, 'Spline');
intensityRef = interp1(time, intensityRef, denseTime, 'Spline');
phaseRef = interp1(time, phaseRef, denseTime, 'Spline');


% plot(denseTime, sqrt(intensity));
% hold on
% plot(denseTime, abs(fftshift(fft(fftshift(sqrt(intensity).*exp(1i*phase)))))/max(abs(fftshift(fft(fftshift(sqrt(intensity).*exp(1i*phase)))))));


% t = -100;
sigma = 200;

%     gate = exp(-((t-denseTime).^2)/(2*sigma.^2));
%     gate = gate/max(gate);
%     windowedE = gate.*sqrt(intensity).*exp(1i*phase);
%     
%     plot(denseTime, sqrt(intensity));
%     hold on
%     plot(denseTime, abs(windowedE)/max(abs(windowedE)));
% 



gabor = [];
tau = linspace(-3000, 3000, 2000);
for t = tau
    gate = exp(-((t-denseTime).^2)/(2*sigma.^2));
    gate = gate/max(gate);
    windowedE = gate.*sqrt(intensity).*exp(1i*phase);
    spectrum = abs(fftshift(fft(fftshift(windowedE))));
%     spectrum = spectrum/max(spectrum);
    gabor = [gabor spectrum];
end

colormap(jet(1024))
h = pcolor(tau, denseTime, gabor);
set(h, 'EdgeColor', 'none');
shading interp
% % caxis([0.01 1])
% ylim([-60 30]);
xlim([-1500 1000]);

% plot(denseTime, sum(gabor,2)/max(sum(gabor,2)));
% hold on
% plot(tau, sum(gabor,1)/max(sum(gabor,1)));


%     gate = exp(-((t-denseTime).^2)/(2*sigma.^2));
%     gate = gate/max(gate);
%     windowedE = gate.*sqrt(intensity);
%     windowedE = windowedE/max(windowedE);
% plot(denseTime, gate)
% hold on
% plot(denseTime, windowedE)
% hold on
% plot(denseTime, sqrt(intensity));


% fig1 = figure();
% 
% subplot(1,3,1)
% plot(denseTime, 6*sqrt(intensity)*pi-8);
% hold on
% plot(denseTime, phase+pi/2);
% plot(denseTime, zeros(size(denseTime))-8);
% xlim([-1500 1000]);
% ylim([-15 15]);
% 
% subplot(1,3,2)
% plot(denseTime, 6*intensityRef*pi-8);
% hold on
% plot(denseTime, phaseRef+pi/2);
% plot(denseTime, zeros(size(denseTime))-8);
% xlim([-1500 1000]);
% ylim([-15 15]);
% 
% subplot(1,3,3)
% plot(denseTime, 6*intensity*pi-8);
% hold on
% plot(denseTime, (phase-phaseRef)+pi/2);
% plot(denseTime, zeros(size(denseTime))-8);
% xlim([-1500 1000]);
% ylim([-15 15]);
% 
% fig2 = figure();
% 
% P1 = sqrt(intensity).*exp(-1i*(2*pi*291*1e12.*denseTime*1e-15 - phase));
% P2 = sqrt(intensity).*exp(1i*(2*pi*291*1e12.*denseTime*1e-15 - phase));
% 
% P1NP = sqrt(intensity).*exp(-1i*(2*pi*291*1e12.*denseTime*1e-15));
% P2NP = sqrt(intensity).*exp(1i*(2*pi*291*1e12.*denseTime*1e-15));
% 
% pulse = (P1 + P2)/2;
% pulseNP = (P1NP + P2NP)/2;
% 
% subplot(1,3,1)
% plot(denseTime, pulse);
% hold on
% plot(denseTime, pulseNP);
% xlim([-1500 1000]);
% 
% P1REF = sqrt(intensityRef).*exp(-1i*(2*pi*291*1e12.*denseTime*1e-15 - phaseRef));
% P2REF = sqrt(intensityRef).*exp(1i*(2*pi*291*1e12.*denseTime*1e-15 - phaseRef));
% 
% P1REFNP = sqrt(intensityRef).*exp(-1i*(2*pi*291*1e12.*denseTime*1e-15));
% P2REFNP = sqrt(intensityRef).*exp(1i*(2*pi*291*1e12.*denseTime*1e-15));
% 
% pulseREF = (P1REF + P2REF)/2;
% pulseREFNP = (P1REFNP + P2REFNP)/2;
% 
% subplot(1,3,2)
% plot(denseTime, pulseREF);
% hold on
% plot(denseTime, pulseREFNP);
% xlim([-1500 1000]);
% 
% P1SPM = sqrt(intensity).*exp(-1i*(2*pi*291*1e12.*denseTime*1e-15 - (phase-phaseRef)));
% P2SPM = sqrt(intensity).*exp(1i*(2*pi*291*1e12.*denseTime*1e-15 - (phase-phaseRef)));
% 
% P1SPMNP = sqrt(intensity).*exp(-1i*(2*pi*291*1e12.*denseTime*1e-15));
% P2SPMNP = sqrt(intensity).*exp(1i*(2*pi*291*1e12.*denseTime*1e-15));
% 
% pulseSPM = (P1SPM + P2SPM)/2;
% pulseSPMNP = (P1SPMNP + P2SPMNP)/2;
% 
% subplot(1,3,3)
% plot(denseTime, pulseSPM);
% hold on
% plot(denseTime, pulseSPMNP);
% xlim([-1500 1000]);
