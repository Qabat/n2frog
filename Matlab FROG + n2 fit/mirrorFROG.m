function symmetrizedFROG = mirrorFROG(rawFROG, delays, omegas, mirror)

N = size(rawFROG, 1);
temporalShift = (delays(2)-delays(1))/2;

if (~strcmp(mirror, 'none'))
    
    rawFROG = abs(ifftshift(ifft(ifftshift(fftshift(fft(fftshift(sqrt(rawFROG), 2), [], 2), 2).*exp(1i.*omegas.*temporalShift), 2), [], 2), 2)).^2;
    
	if (strcmp(mirror, 'left'))
        rawFROG(:,(N/2+1):end) = [];
        symmetrizedFROG = [rawFROG fliplr(rawFROG)];
    elseif (strcmp(mirror, 'right'))
        rawFROG(:,1:N/2) = [];
        symmetrizedFROG = [fliplr(rawFROG) rawFROG];
    elseif (strcmp(mirror, 'both'))
        right = rawFROG(:,(N/2+1):end);
        left = rawFROG(:,1:N/2);
        both = (left + fliplr(right))/2;
        symmetrizedFROG = [both fliplr(both)];
    end
    
    symmetrizedFROG = abs(ifftshift(ifft(ifftshift(fftshift(fft(fftshift(sqrt(symmetrizedFROG), 2), [], 2), 2).*exp(-1i.*omegas.*temporalShift), 2), [], 2), 2)).^2;
    
else
	symmetrizedFROG = rawFROG;
end

symmetrizedFROG = symmetrizedFROG/max(max(symmetrizedFROG));

end

