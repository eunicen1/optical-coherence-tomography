% tiltCorr.m
function OCTTcorr = tiltCorr(OCTMcorr, nLines, usfac)
for I = 1:nLines
 %Every for loop, reference frame will be the middle frame
 [output, ~] = dftregistration(fft2(squeeze(OCTMcorr(:, round(nLines./2),:))),...
 fft2(squeeze(OCTMcorr(:,I,:))), usfac);
 %Assign and save the shifting value along axial direction
 axialShift(I) = round(output(3));
 %Motion correction processing as shifting value along axial direction
 OCTTcorr(:,I,:) = circshift(OCTMcorr(:,I,:), [output(3),0]);
end
end
