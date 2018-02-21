function [omegaFROG, lenDelay, lenOmega, delDelay, delOmega, centOmega] = interpFROG(lambdaFROG)
    normalize = @(M) M/max(max(M));
    c = 299.792458; % speed of light in nm/fs    
    
    lenDelay = lambdaFROG(1);
    lenLambda = lambdaFROG(2);
    delDelay = lambdaFROG(3);
    delLambda = lambdaFROG(4);
    centLambda = lambdaFROG(5);
    
    lambdaFROG(1:5,:) = [];
    omegaFROG = lambdaFROG;

    vectLambda = (centLambda-delLambda*lenLambda/2):delLambda:(centLambda+delLambda*(lenLambda/2-1));
    delOmegaVector = 2*pi*c*delLambda./vectLambda(1:end-1).^2;

	iter = 1;
    vectOmega = 2*pi*c / vectLambda(1); % initialized with first value of omega
    for delOmega = delOmegaVector
        omega = vectOmega(iter) - delOmega;
        vectOmega = [vectOmega omega];
        iter = iter + 1;
    end

    iter = 1;
    for lambda = vectLambda
        omegaFROG(iter,:) = lambdaFROG(iter,:) * lambda^2 / (2*pi*c);
        iter = iter + 1;
    end
        
%     check if integrals in both domain are equal
%     A = sum(lambdaFROG(:,100).*vectLambda')
%     B = sum(omegaFROG(:,100).*vectOmega')
    
    delOmega = (vectOmega(end) - vectOmega(1)) / 255;
    omegaFROG = interp1(vectOmega, omegaFROG, vectOmega(1):delOmega:vectOmega(end),'spline');
    
    omegaFROG = normalize(omegaFROG);    
    lenOmega = length(vectOmega);
    centOmega = 2*pi*c / centLambda;
%     firstOmega = vectOmega(1);
end

