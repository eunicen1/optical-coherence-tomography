fn = 'RawOCT_BM';
loadloc = 'Data';
load(fullfile(loadloc,fn));
addpath('Functions');

%%% Preset parameter %%%
dispMaxOrder = 5;
coeffRange = 10;
depthROI = [45, 300];
calSigWindow = [29, 40];
ref_Frame = 499;
calSigOffIdx = 34;

%%% Reference frame process %%%
ref_RawData = rawOCT_BM(:, :, ref_Frame);
ref_FFTData = fft(hilbert(ref_RawData));
%%% Extract calibration signal with binary window function %%%
%
winFunc = zeros(size(ref_FFTData));
winFunc(calSigWindow(1):calSigWindow(2), :) = 1;
cal_FFTData = ref_FFTData.*winFunc;
cal_RawData = ifft(cal_FFTData);
%hilbert transform - finds analytical signal using real term => generation
%of complex signal 
%plot(abs(ref_FFTData)); %include plot in report
%plot(abs(cal_FFTData)); %include plot
%plot(real(cal_RawData)); %include plot
ref_RawData_Rescaled = resampCal(ref_FFTData, cal_RawData);
ref_FFTData_Rescaled = fft(ref_RawData_Rescaled);

%%% Phase shift estimation & compensation %%%
ref_Ascan = ref_FFTData_Rescaled(:, end/2);
ref_RawData_comp = compPhaseShift(ref_Ascan, ref_FFTData_Rescaled, calSigOffIdx);
ref_FFTData_comp = fft(ref_RawData_comp);

%remove spectral noise
ref_RawData_FPNSub = fpnSubWin(ref_RawData_comp); 
ref_FFTData_FPNSub = fft(ref_RawData_FPNSub);

%%% Dispersion estimation %%%
dispCoeffs = setDispCoeffs(ref_RawData_FPNSub, depthROI, dispMaxOrder, coeffRange);
ref_RawData_DispComp = compDisPhase(ref_RawData_FPNSub, dispMaxOrder, dispCoeffs);
ref_FFTData_DispComp = fft(ref_RawData_DispComp);

plot(ref_RawData); %include plot in report
xlabel("Sampling index [pixel]");
ylabel("Amplitude [a.u.]");
title("Retina OCT fringe data - reference B-scan");
plot(winFunc); %include plot in report
xlabel("Depth index [pixel]");
ylabel("Amplitude [a.u.]");
title("Binary windowing for calibration signal");

plot(real(ref_RawData_Rescaled));%include plot in report
xlabel("Sampling index [pixel]");
ylabel("Amplitude [a.u.]");
title("Resampled retina OCT fringe data - reference B-scan");

plot(real(ref_RawData_FPNSub));%include plot in report
xlabel("Sampling index [pixel]");
ylabel("Amplitude [a.u.]");
title("FPN Removed OCT fringe data");

OCTImg1 = 20.*log10(abs(ref_FFTData_Rescaled));
imagesc(imadjust(mat2gray(OCTImg1(1:end/2,:)))); colormap("gray")
axis off;

OCTImg2 = 20.*log10(abs(ref_FFTData_FPNSub));
imagesc(imadjust(mat2gray(OCTImg2(1:end/2,:)))); colormap("gray")
axis off;

OCTImg3 = 20.*log10(abs(ref_FFTData_comp));
imagesc(imadjust(mat2gray(OCTImg3(1:end/2,:)))); colormap("gray")
axis off;

OCTImg4 = 20.*log10(abs(ref_FFTData_DispComp));
imagesc(imadjust(mat2gray(OCTImg4(1:end/2,:)))); colormap("gray")
axis off;

% Splitting spectrum window 
Fraction = 2;
numPoints = size(ref_RawData_DispComp,1);
numBands = 2*Fraction-1;
winWidth = round(numPoints/Fraction);
scale = sqrt(sum(hanning(numPoints).^2)/sum(hanning(winWidth).^2));
for I = 1:numBands
    wins(:, I) = circshift(cat(1, hanning(winWidth)...
        , zeros(numPoints-winWidth, 1))...
        , round((I-1)*numPoints/(numBands+1)))*scale;
end

splitWins = permute(shiftdim(wins, -1), [2,1,3]);
sub_split = bsxfun(@times, ref_RawData_DispComp, splitWins); %*** check this

plot(splitWins(:,:,1)); hold on; plot(splitWins(:,:,2)); plot(splitWins(:,:,3));
title("Split spectrum windows");
xlabel("Sampling index [pixel]");
ylabel("Amplitude [a.u.]");

plot(real(sub_split(:,:,1)));
title("Split spectrum OCT fringe data - 1");
xlabel("Sampling index [pixel]");
ylabel("Amplitude [a.u.]");

plot(real(sub_split(:,:,2)));
title("Split spectrum OCT fringe data - 2");
xlabel("Sampling index [pixel]");
ylabel("Amplitude [a.u.]");

plot(real(sub_split(:,:,3)));
title("Split spectrum OCT fringe data - 3");
xlabel("Sampling index [pixel]");
ylabel("Amplitude [a.u.]");

%need array to store full dimension, split 1, split 2, and split 3
imagesc(squeeze(20.*log10(abs(fft(ref_RawData_DispComp)))));colormap(gray); axis off;
imagesc(squeeze(20.*log10(abs(fft(sub_split(:,:,1))))));colormap(gray); axis off;
imagesc(squeeze(20.*log10(abs(fft(sub_split(:,:,2))))));colormap(gray); axis off;
imagesc(squeeze(20.*log10(abs(fft(sub_split(:,:,3))))));colormap(gray); axis off;

ref_rawData_Full = rawHann(ref_RawData_DispComp);
ref_FFTData_Full = fft(ref_RawData_DispComp);
ref_rawData_Split = sub_split;

OCTImg5 = 20.*log10(abs(ref_FFTData_Full));
imagesc(imadjust(mat2gray(OCTImg5(1:end/2,:)))); colormap("gray")
axis off;

ref_FFTData_Full = fft(ref_rawData_Full);
ref_FFTData_Split = fft(ref_rawData_Split);

ref_OCT_Full = 20.*log10(abs(ref_FFTData_Full(1:end/2,:)));
ref_OCT_Split_1 = 20.*log10(abs(ref_FFTData_Split(1:end/2,:,1)));
ref_OCT_Split_2 = 20.*log10(abs(ref_FFTData_Split(1:end/2,:,2)));
ref_OCT_Split_3 = 20.*log10(abs(ref_FFTData_Split(1:end/2,:,3)));

imagesc(imadjust(mat2gray(squeeze(ref_OCT_Split_1)))); colormap(gray); axis off;
imagesc(imadjust(mat2gray(squeeze(ref_OCT_Split_2)))); colormap(gray); axis off;
imagesc(imadjust(mat2gray(squeeze(ref_OCT_Split_3)))); colormap(gray); axis off;

for FrameNum = 1:size(rawOCT_BM, 3)
    rawData = rawOCT_BM(:,:,FrameNum);
    ref_FFTData_n = fft(hilbert(rawData));
    
    %resampling
    ref_RawData_Rescaled_n = resampCal(ref_FFTData_n, cal_RawData);
    ref_FFTData_Rescaled_n = fft(ref_RawData_Rescaled_n);
    
    %comp Phase Shift    
    ref_RawData_comp = compPhaseShift(ref_Ascan, ref_FFTData_Rescaled_n, calSigOffIdx);
    
    %remove spectral noise
    ref_RawData_FPNSub_n = fpnSubWin(ref_RawData_comp);
    %proccess_n = rawHann(ref_RawData_FPNSub_n);
    
    %spectral shift correction 
    ref_RawData_comp_n = compDisPhase(ref_RawData_FPNSub_n, dispMaxOrder, dispCoeffs);
    raw_Full = rawHann(ref_RawData_comp_n);
    
    % Splitting spectrum window
    raw_split_n = bsxfun(@times, ref_RawData_comp_n, splitWins); %*** check this

    fftD = fft(raw_Full); %full dimension
    fftD_s1 = squeeze(fft(raw_split_n(:,:,1))); % split 1
    fftD_s2 = squeeze(fft(raw_split_n(:,:,2))); % split 2
    fftD_s3 = squeeze(fft(raw_split_n(:,:,3))); % split 3
    
    ProcdD(:,:,FrameNum) = fftD(depthROI(1):depthROI(2),:); %full dimension
    ProcdD_s1(:,:,FrameNum) = fftD_s1(depthROI(1):depthROI(2),:); % split 1
    ProcdD_s2(:,:,FrameNum) = fftD_s2(depthROI(1):depthROI(2),:); % split 2
    ProcdD_s3(:,:,FrameNum) = fftD_s3(depthROI(1):depthROI(2),:); % split 3
end

imagesc(abs(20.*log10(ProcdD(:,:,ref_Frame)))); colormap(gray); axis off;
imagesc(abs(20.*log10(ProcdD_s1(:,:,ref_Frame)))); colormap(gray); axis off;
imagesc(abs(20.*log10(ProcdD_s2(:,:,ref_Frame)))); colormap(gray); axis off;
imagesc(abs(20.*log10(ProcdD_s3(:,:,ref_Frame)))); colormap(gray); axis off;
usfac = 1;
bm = 2; 
% apply motion correction (full)
nFrames_full = size(ProcdD, 3);
[raw_motionCorr_full, x_axialShift, y_axialShift] = localmotionCorr(nFrames_full, bm, ProcdD, usfac);
ref_OCT = 20.*log10(abs(squeeze(raw_motionCorr_full(:,250,:))));
imagesc(ref_OCT); axis off; colormap(gray);
ref_OCT = 20.*log10(abs(squeeze(raw_motionCorr_full(:,:,250))));
imagesc(ref_OCT); axis off; colormap(gray);

nl_full = size(ProcdD, 2);
tcorr_full = tiltCorr(raw_motionCorr_full, nl_full, usfac);
t_full = squeeze(max(abs(raw_motionCorr_full(26:180,:,:)),[],1));
imagesc(t_full); axis off; colormap(gray);

%local motion correction plot
plot(x_axialShift);
xlabel("Frame Index");
ylabel("X Axial Shifts [pixel]");
title("X Local Motion Correction Plot");

plot(y_axialShift);
xlabel("Frame Index");
ylabel("Y Axial Shifts [pixel]");
title("Y Local Motion Correction Plot");

% apply motion correction (s1)
nFrames_s1 = size(ProcdD_s1, 3);
[raw_motionCorr_s1, x_axialShift1, y_axialShift1] = localmotionCorr(nFrames_s1, bm, ProcdD_s1, usfac);
ref_OCT_s1 = 20.*log10(abs(squeeze(raw_motionCorr_s1(:,250,:))));
imagesc(ref_OCT_s1); axis off; colormap(gray);
ref_OCT_s1 = 20.*log10(abs(squeeze(raw_motionCorr_s1(:,:,250))));
imagesc(ref_OCT_s1); axis off; colormap(gray);

nl_s1 = size(ProcdD_s1, 2);
tcorr_s1 = tiltCorr(raw_motionCorr_s1, nl_s1, usfac);
t_s1 = squeeze(max(abs(raw_motionCorr_s1(26:180,:,:)),[],1));
imagesc(t_s1); axis off; colormap(gray);


% apply motion correction (s2)
nFrames_s2 = size(ProcdD_s2, 3);
[raw_motionCorr_s2, x_axialShift2, y_axialShift2] = localmotionCorr(nFrames_s2, bm, ProcdD_s2, usfac);
ref_OCT_s2 = 20.*log10(abs(squeeze(raw_motionCorr_s2(:,250,:))));
imagesc(ref_OCT_s2); axis off; colormap(gray);
ref_OCT_s2 = 20.*log10(abs(squeeze(raw_motionCorr_s2(:,:,250))));
imagesc(ref_OCT_s2); axis off; colormap(gray);

nl_s2 = size(ProcdD_s2, 2);
tcorr_s2 = tiltCorr(raw_motionCorr_s2, nl_s2, usfac);
t_s2 = squeeze(max(abs(raw_motionCorr_s2(26:180,:,:)),[],1));
imagesc(t_s2); axis off; colormap(gray);

% apply motion correction (s3)
nFrames_s3 = size(ProcdD_s3, 3);
[raw_motionCorr_s3, x_axialShift3, y_axialShift3] = localmotionCorr(nFrames_s3, bm, ProcdD_s3, usfac);
ref_OCT_s3 = 20.*log10(abs(squeeze(raw_motionCorr_s3(:,250,:))));
imagesc(ref_OCT_s3); axis off; colormap(gray);
ref_OCT_s3 = 20.*log10(abs(squeeze(raw_motionCorr_s3(:,:,250))));
imagesc(ref_OCT_s3); axis off; colormap(gray);

nl_s3 = size(ProcdD_s3, 2);
tcorr_s3 = tiltCorr(raw_motionCorr_s3, nl_s3, usfac);
t_s3 = squeeze(max(abs(raw_motionCorr_s3(26:180,:,:)),[],1));
imagesc(t_s3); axis off; colormap(gray);

% OCT Angiography processing
for I = 1:bm:nFrames_full
   K=((I-1)/bm)+1;
   xconj = ProcdD(:, :, I+1).*conj(ProcdD(:,:,I));
   bulkoff = repmat(angle(sum(xconj,1)), [size(xconj,1) 1]);
   bscan1 = ProcdD(:, :, I);
   bscan2 = ProcdD(:, :, I+1).*exp(-1j*bulkoff);
   
   % average OCT
   avgOCT(:, :, K) = (bscan1 + bscan2)./2;
   % variance OCT
   varOCT(:, :, K) = abs(var(cat(3, bscan1, bscan2), 0, 3));
   % subtraction 
   sub(:, :, K) = abs(bscan1 - bscan2);
   % decorrelation 
   dec(:, :, K) = 1 - ((abs(bscan1).*abs(bscan2))...
       ./((abs(bscan1).^2 + abs(bscan2).^2)./2));
end
ref_OCT_avg = 20.*log10(abs(squeeze(avgOCT(:,250,:))));
imagesc(ref_OCT_avg); axis off; colormap(gray);
ref_OCT_avg = 20.*log10(abs(squeeze(avgOCT(:,:,250))));
imagesc(ref_OCT_avg); axis off; colormap(gray);

ref_OCT_var = 20.*log10(abs(squeeze(varOCT(:,250,:))));
imagesc(ref_OCT_var); axis off; colormap(gray);
ref_OCT_var = 20.*log10(abs(squeeze(varOCT(:,:,250))));
imagesc(ref_OCT_var); axis off; colormap(gray);

ref_OCT_sub = 20.*log10(abs(squeeze(sub(:,250,:))));
imagesc(ref_OCT_sub); axis off; colormap(gray);
ref_OCT_sub = 20.*log10(abs(squeeze(sub(:,:,250))));
imagesc(ref_OCT_sub); axis off; colormap(gray);

ref_OCT_dec = 20.*log10(abs(squeeze(dec(:,250,:))));
imagesc(ref_OCT_dec); axis off; colormap(gray);
ref_OCT_dec = 20.*log10(abs(squeeze(dec(:,:,250))));
imagesc(ref_OCT_dec); axis off; colormap(gray);


% avgOCT log scale
avgOCT_log = 20.*log10(abs(flipud(avgOCT)));
% OCTA with intensity thresholding
OCTA_Var = flipud(varOCT); OCTA_Var(avgOCT_log<90) = min(varOCT(:));
OCTA_Sub = flipud(sub); OCTA_Sub(avgOCT_log<90) = min(sub(:));
OCTA_Dec = flipud(dec); OCTA_Dec(avgOCT_log<90) = min(dec(:));
% smoothening and global motion correction
avgOCTA_var = mov2DAvg(OCTA_Var, [2, 2]);
avgOCTA_sub = mov2DAvg(OCTA_Sub, [2, 2]);
avgOCTA_dec = mov2DAvg(OCTA_Dec, [2, 2]);

for I = 1:nFrames_full/2
   [output, ~] = dftregistration(fft2(avgOCT_log(:,:,nFrames_full/2)), fft2(avgOCT_log(:,:,I)), usfac);
   yshift_global(I) = round(output(3));
   avgOCTA_mcorr(:, :, I) = circshift(avgOCT_log(:,:,I), [output(3), 0]);
   varOCTA_mcorr(:, :, I) = circshift(avgOCTA_var(:,:,I), [output(3), 0]);
   subOCTA_mcorr(:, :, I) = circshift(avgOCTA_sub(:,:,I), [output(3), 0]);
   decOCTA_mcorr(:, :, I) = circshift(avgOCTA_dec(:,:,I), [output(3), 0]);  
end

nLines = size(ProcdD, 2);
avgOCTA_tcorr = tiltCorr(avgOCTA_mcorr, nLines, usfac);
varOCTA_tcorr = tiltCorr(varOCTA_mcorr, nLines, usfac);
subOCTA_tcorr = tiltCorr(subOCTA_mcorr, nLines, usfac);
decOCTA_tcorr = tiltCorr(decOCTA_mcorr, nLines, usfac);

plot(yshift_global);
xlabel("Frame Index");
ylabel("Y Global Shifts [pixel]");
title("Y Global Motion Correction Plot");

ref_OCT_avg_corr = 20.*log10(abs(squeeze(avgOCTA_mcorr(:,250,:))));
imagesc(ref_OCT_avg_corr); axis off; colormap(gray);
ref_OCT_avg_corr = 20.*log10(abs(squeeze(avgOCTA_mcorr(:,:,250))));
imagesc(ref_OCT_avg_corr); axis off; colormap(gray);
average = abs(avgOCTA_mcorr);

ref_OCT_var_corr = 20.*log10(abs(squeeze(varOCTA_mcorr(:,:,250))));
imagesc(ref_OCT_var_corr); axis off; colormap(gray);
variance = abs(varOCTA_mcorr);

ref_OCT_sub_corr = 20.*log10(abs(squeeze(subOCTA_mcorr(:,:,250))));
imagesc(ref_OCT_sub_corr); axis off; colormap(gray);
subtraction = abs(subOCTA_mcorr);

ref_OCT_dec_corr = 20.*log10(abs(squeeze(decOCTA_mcorr(:,:,250))));
imagesc(ref_OCT_dec_corr); axis off; colormap(gray);
decor = abs(decOCTA_mcorr);

%tilt correct, for octa no need for log scales, enface for split data

for I = 1:bm:nFrames_full
   K=((I-1)/bm)+1;
   xconj_s1 = ProcdD_s1(:, :, I+1).*conj(ProcdD_s1(:,:,I));
   bscan1_s1 = ProcdD_s1(:, :, I);
   bscan2_s1 = ProcdD_s1(:, :, I+1);
   
   xconj_s2 = ProcdD_s2(:, :, I+1).*conj(ProcdD_s2(:,:,I));
   bscan1_s2 = ProcdD_s2(:, :, I);
   bscan2_s2 = ProcdD_s2(:, :, I+1);
 
   xconj_s3 = ProcdD_s3(:, :, I+1).*conj(ProcdD_s3(:,:,I));
   bscan1_s3 = ProcdD_s3(:, :, I);
   bscan2_s3 = ProcdD_s3(:, :, I+1);
 
   % decorrelation 
   dec_s1(:, :, K) = 1 - ((abs(bscan1_s1).*abs(bscan2_s1))...
       ./((abs(bscan1_s1).^2 + abs(bscan2_s1).^2)./2));
   dec_s2(:, :, K) = 1 - ((abs(bscan1_s2).*abs(bscan2_s2))...
       ./((abs(bscan1_s2).^2 + abs(bscan2_s2).^2)./2));
   dec_s3(:, :, K) = 1 - ((abs(bscan1_s3).*abs(bscan2_s3))...
       ./((abs(bscan1_s3).^2 + abs(bscan2_s3).^2)./2));
end

average_dec = (dec_s1 + dec_s2 + dec_s3)/3;

decor_avg = abs(decOCTA_mcorr);
nLines_davg = size(average_dec, 2);
decor_tilt_avg = tiltCorr(decor_avg, nLines_davg, usfac);

enface_avg = squeeze(max(abs(avgOCTA_tcorr(26:180,:,:)),[],1));
imagesc(enface_avg); axis off; colormap(gray);

enface_var = squeeze(max(abs(varOCTA_tcorr(26:180,:,:)),[],1));
imagesc(enface_var); axis off; colormap(gray);

enface_sub = squeeze(max(abs(subOCTA_tcorr(26:180,:,:)),[],1));
imagesc(enface_sub); axis off; colormap(gray);

enface_dec = squeeze(max(abs(decOCTA_tcorr(26:180,:,:)),[],1));
imagesc(enface_dec); axis off; colormap(gray);

xy_avg = squeeze(max(abs(avgOCTA_tcorr(:,:,26:180)),[],3));
imagesc(xy_avg); axis off; colormap(gray);

xy_var = squeeze(max(abs(varOCTA_tcorr(:,:,26:180)),[],3));
imagesc(xy_var); axis off; colormap(gray);

xy_sub = squeeze(max(abs(subOCTA_tcorr(:,:,26:180)),[],3));
imagesc(xy_sub); axis off; colormap(gray);

xy_dec = squeeze(max(abs(decOCTA_tcorr(:,:,26:180)),[],3));
imagesc(xy_dec); axis off; colormap(gray);

enface_ssada_davg = squeeze(max(abs(decor_tilt_avg(26:180,:,:)),[],1));
imagesc(enface_ssada_davg); axis off; colormap(gray);

xy_ssada_davg = squeeze(max(abs(decor_tilt_avg(:,:,26:180)),[],3));
imagesc(xy_ssada_davg); axis off; colormap(gray);
