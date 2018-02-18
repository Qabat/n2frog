function [ omegaFROG ] = lambdatoomega( lambdaFROG )
%LAMBDATOOMEGA Summary of this function goes here
%   function converting input FROG matrix (femtosoft compatible)
%   from wavelength domain to frequency domain
%   femtosoft: TReadIn.GridSpectrum
%   zrobic sprawdzenie ze calka z widma w lambda i w omega musza byc równe

    c = 299.792458; % speed of light in nm/fs

    lenDelay = lambdaFROG(1);
    lenLambda = lambdaFROG(2);
    delDelay = lambdaFROG(3);
    delLambda = lambdaFROG(4);
    centLambda = lambdaFROG(5);
    
    lambdaFROG(1:5,:) = [];

    centOmega = c/centLambda;
    
    testSpectrum = lambdaFROG(128,:);
    plot(testSpectrum);
    
    for i=1:lenLambda
        delOmega = (c/lambda^2) * delLambda;
        omega = centOmega + (i-N/2)*delOmega;
        lambda = c / omega;
        
    end
    
    omegaFROG = lambdaFROG;
end

