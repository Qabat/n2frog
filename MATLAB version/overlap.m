pulse1 = dlmread('../../overlap/1.txt');
pulse2 = dlmread('../../overlap/3.txt');

time1 = pulse1(:,1);
time2 = pulse2(:,1);

intensity1 = pulse1(:,2);
intensity2 = pulse2(:,2);

phase1 = pulse1(:,3);
phase2 = pulse1(:,3);

omegas1 = pulse1(:,4)/1000;

bestToverlap = 100;
bestTau = 0;
for tau = -200:200
    Pulse = sqrt(intensity1).*exp(-1i*phase1);
    Pulse = ifftshift(ifft(ifftshift(fftshift(fft(fftshift(Pulse))).*exp(1i.*omegas1.*tau))));
    [~, centerPulse] = max(abs(Pulse).^2);
    Pulse = Pulse./Pulse(centerPulse);
    Tintensity1 = abs(Pulse).^2;
    
    Toverlap = trapz(abs(Tintensity1-intensity2));
    if (Toverlap < bestToverlap)
        bestToverlap = Toverlap;
        bestTau = tau;
    end
end
disp(bestTau);

    Pulse = sqrt(intensity1).*exp(-1i*phase1);
    Pulse = ifftshift(ifft(ifftshift(fftshift(fft(fftshift(Pulse))).*exp(1i.*omegas1.*bestTau))));
    [~, centerPulse] = max(abs(Pulse).^2);
    Pulse = Pulse./Pulse(centerPulse);
intensity1 = abs(Pulse).^2;
phase1 = angle(Pulse);

figure('Position',[150 75 1600 900]);
plot(time1, intensity1)
hold on 
plot(time2, intensity2)
hold on
plot(time1, abs(intensity1-intensity2))
% xlim([-700 700]);

disp(trapz(abs(intensity1-intensity2)));