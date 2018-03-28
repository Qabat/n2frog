function symmetrizedFROG = mirrorFROG(rawFROG, mirror)

N = size(rawFROG, 1);

if (~strcmp(mirror, 'none'))
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
else
    symmetrizedFROG = rawFROG;
end

end

