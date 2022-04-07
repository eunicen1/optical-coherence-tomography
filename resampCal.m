%% Resampling function
% resamplecal.m
function rawResc = resampCal(fftd, cald) % fft data and calibration
raw_r = real(ifft(fftd));
raw_i = imag(ifft(fftd));
rawResc = zeros(size(fftd));

idx = (1:length(cald))';

for i = 1:size(fftd, 2)
   %Phase extraction
   phase = unwrap(angle(cald(:,i)));
   phaseNorm = phase - min(phase(:));
   phaseNorm = phaseNorm ./ max(phaseNorm(:));
   
   %Polynomial fitting 
   p = fit(phaseNorm, idx, 'cubicinterp');
   
   %Linear phase
   phaseLin = linspace(0,1,length(cald))';
   
   %Linear sampling index along wave numbers
   idxLinK = p(phaseLin);
   
   rawResc(:,i) = (interp1([1:size(raw_r, 1)]', raw_r(:,i), idxLinK, 'spline'))...
       +1j .* (interp1([1:size(raw_i,1)]', raw_i(:,i), idxLinK, 'spline'));
end

end
