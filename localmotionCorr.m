%% Motion correction
% localmotionCorr.m
function [OCTMcorr, x_axialShift, y_axialShift] = localmotionCorr(nFrames, bm, procdOCT_BM, usfac)
%Set the shifting variable to save and analyze
x_axialShift=zeros([nFrames 1]);
y_axialShift=zeros([nFrames 1]);
cplxOCT_ROI = procdOCT_BM; %(depthROI(1):depthROI(2), :, :);
OCTMcorr = cplxOCT_ROI; %cplxOCT_ROI;
for I = 1:bm:nFrames
    %Every for loop, reference frame will be the middle frame

    [output, ~] = dftregistration(fft2(20.*log10(abs(cplxOCT_ROI(:,:,I)))), ...
        fft2(20.*log10(abs(cplxOCT_ROI(:,:,I+1)))), usfac);
    %Assign and save the shifting value along axial direction
    x_axialShift(I) = round(output(4));
    y_axialShift(I) = round(output(3));
    %Motion correction processing as shifting value along axial direction
    OCTMcorr(:, :, I+1) = circshift(cplxOCT_ROI(:,:,I+1), [output(3), output(4)]);
end
end