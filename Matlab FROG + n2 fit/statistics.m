clear;
n2results = dlmread('../../fits/fs/results/fs.csv');

power = n2results(:,1);
n2 = n2results(:,2);
error = n2results(:,3);

errorbar(power, n2, error,'o')
hold on
xlim([10 180]);
ylim([0 20]);
xlabel('Energy [\muJ]');
ylabel('n_2 [10^{-16} cm^{2}/W]');
title('n_2 values for all measured samples');

weights = 1./error.^2;

finalResult = sum(weights.*n2)/sum(weights);
sigmaInternal = 1/sqrt(sum(weights));
sigmaExternal = sigmaInternal * sqrt(sum(weights.*(n2-finalResult).^2)./(length(power)-1));

if (sigmaInternal > sigmaExternal)
    finalError = sigmaInternal;
else
    finalError = sigmaExternal;
end

plot(power, ones(size(power)).*finalResult)

disp(['Final fs n2 value: ' num2str(finalResult) ' +/- ' num2str(finalError)]);

clear;
n2results = dlmread('../../fits/yag/results/yag.csv');

power = n2results(:,1);
n2 = n2results(:,2);
error = n2results(:,3);

errorbar(power, n2, error,'o')
hold on

weights = 1./error.^2;

finalResult = sum(weights.*n2)/sum(weights);
sigmaInternal = 1/sqrt(sum(weights));
sigmaExternal = sigmaInternal * sqrt(sum(weights.*(n2-finalResult).^2)./(length(power)-1));

if (sigmaInternal > sigmaExternal)
    finalError = sigmaInternal;
else
    finalError = sigmaExternal;
end

plot(power, ones(size(power)).*finalResult)

disp(['Final yag n2 value: ' num2str(finalResult) ' +/- ' num2str(finalError)]);

clear;
n2results = dlmread('../../fits/caf2/results/caf2.csv');

power = n2results(:,1);
n2 = n2results(:,2);
error = n2results(:,3);

errorbar(power, n2, error,'o')
hold on

weights = 1./error.^2;

finalResult = sum(weights.*n2)/sum(weights);
sigmaInternal = 1/sqrt(sum(weights));
sigmaExternal = sigmaInternal * sqrt(sum(weights.*(n2-finalResult).^2)./(length(power)-1));

if (sigmaInternal > sigmaExternal)
    finalError = sigmaInternal;
else
    finalError = sigmaExternal;
end

plot(power, ones(size(power)).*finalResult)

disp(['Final caf2 n2 value: ' num2str(finalResult) ' +/- ' num2str(finalError)]);

clear;
n2results = dlmread('../../fits/bibo/results/bibo.csv');

power = n2results(:,1);
n2 = n2results(:,2);
error = n2results(:,3);

errorbar(power, n2, error,'o')
hold on

weights = 1./error.^2;

finalResult = sum(weights.*n2)/sum(weights);
sigmaInternal = 1/sqrt(sum(weights));
sigmaExternal = sigmaInternal * sqrt(sum(weights.*(n2-finalResult).^2)./(length(power)-1));

if (sigmaInternal > sigmaExternal)
    finalError = sigmaInternal;
else
    finalError = sigmaExternal;
end

plot(power, ones(size(power)).*finalResult)

disp(['Final bibo n2 value: ' num2str(finalResult) ' +/- ' num2str(finalError)]);

clear;
n2results = dlmread('../../fits/CALC no/results/calc no.csv');

power = n2results(:,1);
n2 = n2results(:,2);
error = n2results(:,3);

errorbar(power, n2, error,'o')
hold on

weights = 1./error.^2;

finalResult = sum(weights.*n2)/sum(weights);
sigmaInternal = 1/sqrt(sum(weights));
sigmaExternal = sigmaInternal * sqrt(sum(weights.*(n2-finalResult).^2)./(length(power)-1));

if (sigmaInternal > sigmaExternal)
    finalError = sigmaInternal;
else
    finalError = sigmaExternal;
end

plot(power, ones(size(power)).*finalResult)

disp(['Final calc no n2 value: ' num2str(finalResult) ' +/- ' num2str(finalError)]);

clear;
n2results = dlmread('../../fits/CALC ne/results/calc ne.csv');

power = n2results(:,1);
n2 = n2results(:,2);
error = n2results(:,3);

errorbar(power, n2, error,'o')
hold on

weights = 1./error.^2;

finalResult = sum(weights.*n2)/sum(weights);
sigmaInternal = 1/sqrt(sum(weights));
sigmaExternal = sigmaInternal * sqrt(sum(weights.*(n2-finalResult).^2)./(length(power)-1));

if (sigmaInternal > sigmaExternal)
    finalError = sigmaInternal;
else
    finalError = sigmaExternal;
end

plot(power, ones(size(power)).*finalResult)

disp(['Final calc ne value: ' num2str(finalResult) ' +/- ' num2str(finalError)]);


clear;
n2results = dlmread('../../fits/YVO no/results/yvo no.csv');

power = n2results(:,1);
n2 = n2results(:,2);
error = n2results(:,3);

errorbar(power, n2, error,'o')
hold on

weights = 1./error.^2;

finalResult = sum(weights.*n2)/sum(weights);
sigmaInternal = 1/sqrt(sum(weights));
sigmaExternal = sigmaInternal * sqrt(sum(weights.*(n2-finalResult).^2)./(length(power)-1));

if (sigmaInternal > sigmaExternal)
    finalError = sigmaInternal;
else
    finalError = sigmaExternal;
end

plot(power, ones(size(power)).*finalResult)

disp(['Final yvo no value: ' num2str(finalResult) ' +/- ' num2str(finalError)]);


clear;
n2results = dlmread('../../fits/YVO ne/results/YVO ne.csv');

power = n2results(:,1);
n2 = n2results(:,2);
error = n2results(:,3);

errorbar(power, n2, error,'o')
hold on

weights = 1./error.^2;

finalResult = sum(weights.*n2)/sum(weights);
sigmaInternal = 1/sqrt(sum(weights));
sigmaExternal = sigmaInternal * sqrt(sum(weights.*(n2-finalResult).^2)./(length(power)-1));

if (sigmaInternal > sigmaExternal)
    finalError = sigmaInternal;
else
    finalError = sigmaExternal;
end

plot(power, ones(size(power)).*finalResult)

disp(['Final yvo ne value: ' num2str(finalResult) ' +/- ' num2str(finalError)]);