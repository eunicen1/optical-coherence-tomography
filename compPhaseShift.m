%% Spectral shift compensation
% compPhaseShift.m
function rawPhaseComp = compPhaseShift(ref_fft_1d, fft_2d, depthi)
%Reference A-scan in OCT volume data
refAscan = ref_fft_1d;

%Complex conjugate multiplication
cxConjX = fft_2d...
    .* repmat(conj(refAscan), [1 size(fft_2d,2)]);

%Calibration signal depth index
calsdepth = depthi;

%Phase slope
phasem = (angle(cxConjX(calsdepth,:))/calsdepth)...
    .*linspace(1, length(refAscan), length(refAscan))';

%Phase compensation with phase slope
for i = 1:size(fft_2d, 2)
    fftPhaseComp(:,i) = fft_2d(:,i)...
        .*exp(-1j*phasem(:,i));
end
rawPhaseComp = ifft(fftPhaseComp);
end