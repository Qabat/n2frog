clc;
clear;
close all;

lambdaFROG = dlmread('../testfrog/60.txt');
omegaFROG = interpFROG(lambdaFROG);
    lenDelay = lambdaFROG(1);
    lenLambda = lambdaFROG(2);
    delDelay = lambdaFROG(3);
    delLambda = lambdaFROG(4);
    centLambda = lambdaFROG(5);
vectLambda = (centLambda-delLambda*lenLambda/2):delLambda:(centLambda+delLambda*(lenLambda/2-1));
% lambdaFROG(1:5,:) = [];
%  subplot(1,2,1)
%  image( lambdaFROG*150);
%  
%  subplot(1,2,2)
%  image( omegaFROG*150);
 
     
dlmwrite('../testfrog/60 omega.txt', omegaFROG, '\t');