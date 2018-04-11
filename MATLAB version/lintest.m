%   flat test 
%
% sig = [-25 -16 -9 -4 -1 0 -1 -4 -9 -16 -25];
% trend = [0 2 4 6 8 10 12 14 16 18 20]-10;
% sum = sig + trend;
% 
% 
% plot(trend)
% hold on
% plot(sig)
% hold on
% plot(sum)
% hold on
% plot(detrend(sum)-max(detrend(sum)), '*')

pulse0 = dlmread('../../porownanie faz do odjecia liniowego skladnika/0.txt');
pulse1 = dlmread('../../porownanie faz do odjecia liniowego skladnika/1.txt');
pulse11 = dlmread('../../porownanie faz do odjecia liniowego skladnika/-1.txt');

phase0 = pulse0(:,3);
phase1 = pulse1(:,3);
phase11 = pulse11(:,3);

intensity0 = pulse0(:,2);
intensity1 = pulse1(:,2);
intensity11 = pulse11(:,2);

phase0(intensity0 < 0.1) = [];
phase1(intensity1 < 0.1) = [];
phase11(intensity11 < 0.1) = [];

phase0 = phase0 - max(phase0);
phase1 = phase1 - max(phase1);
phase11 = phase11 - max(phase11);

figure()
hold on
plot(phase0+1);
plot(phase1+1);
plot(phase11+1);

Dphase0 = detrend(phase0);
Dphase1 = detrend(phase1);
Dphase11 = detrend(phase11);

Dphase0 = Dphase0 - max(Dphase0);
Dphase1 = Dphase1 - max(Dphase1);
Dphase11 = Dphase11 - max(Dphase11);

plot(Dphase0+1,'*');
plot(Dphase1+1,'*');
plot(Dphase11+1,'*');
% hold off

intensity0(intensity0 < 0.1) = [];
intensity1(intensity1 < 0.1) = [];
intensity11(intensity11 < 0.1) = [];

% figure()
% hold on

intensity0 = intensity0 - intensity0(1);
intensity1 = intensity1 - intensity1(1);
intensity11 = intensity11 - intensity11(1);

intensity0 = intensity0 / max(intensity0);
intensity1 = intensity1 / max(intensity1);
intensity11 = intensity11 / max(intensity11);

plot(intensity0,'--');
plot(intensity1,'--');
plot(intensity11,'--');

Dintensity0 = detrend(intensity0);
Dintensity1 = detrend(intensity1);
Dintensity11 = detrend(intensity11);

Dintensity0 = Dintensity0 - Dintensity0(1);
Dintensity1 = Dintensity1 - Dintensity1(1);
Dintensity11 = Dintensity11 - Dintensity11(1);

Dintensity0 = Dintensity0 / max(Dintensity0);
Dintensity1 = Dintensity1 / max(Dintensity1);
Dintensity11 = Dintensity11 / max(Dintensity11);

plot(Dintensity0,'o');
plot(Dintensity1,'o');
plot(Dintensity11,'o');
hold off