clc;
clear;
close all;

Pulse = dlmread('..\frog trace for test\result.txt');

Time = Pulse(:,1);
Intensity = Pulse(:,2);
Phase = Pulse(:,3);
complexPulse = sqrt(Intensity).*exp(1i.*Phase);

[computedFROG, EF] = makeFROG(complexPulse);

% rekonstrukcja frogowa
for i=1:1
    % sprawdzic tu czy ja nie daje jako input do kolejknych algorytmow tego
    % co poprzedni wyrzucil, to by wyjasnialo czemu nie widac wplywu
    % bootstrapa
    
% startPulse can be switched for [] to use default pulse
[Pt, Fr, G, iter] = mainFROG(computedFROG, 0.0000001, 100, 0, 6.515, 0);

recIntensity = abs(Pt.^2);
recIntensity = recIntensity/max(recIntensity);
recPhase = angle(Pt);
%recPhase(recIntensity<0.1) = 0;

outputFile = [Time, recIntensity, recPhase];
method = 'svd';
dlmwrite(['.\output_' method '\'  num2str(i)  '.txt'],outputFile,'\t');

end