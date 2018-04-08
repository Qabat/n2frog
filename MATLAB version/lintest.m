
sig = [-25 -16 -9 -4 -1 0 -1 -4 -9 -16 -25];
trend = [0 2 4 6 8 10 12 14 16 18 20]-10;
sum = sig + trend;


plot(trend)
hold on
plot(sig)
hold on
plot(sum)
hold on
plot(detrend(sum)-max(detrend(sum)), '*')