n2results = dlmread('../../fits/yag/yag.txt');

power = n2results(:,1);
n2 = n2results(:,2);
error = n2results(:,3);

errorbar(power, n2, error,'o')
hold on
xlim([0 120]);

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
ylim([0 10]);

disp(['Final n2 value: ' num2str(finalResult) ' +/- ' num2str(finalError)]);