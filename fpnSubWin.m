%% FPN substraction & windowing 
% fpnSubWin.m
function fpnRem = fpnSubWin(raw)
    fpnRem = raw...
        - (repmat(median(real(raw), 2), [1, size(raw,2)])...
        +1j.*repmat(median(imag(raw),2), [1, size(raw,2)]));
    FPN = (repmat(median(real(raw), 2), [1, size(raw,2)])...
        +1j.*repmat(median(imag(raw),2), [1, size(raw,2)]));
    
    %plot(real(FPN)); % to plot

end