%% Hanning windowing processing 
% rawHann.m
function rawHanWin = rawHann(raw)
    rawHanWin = raw ... 
    .* repmat(hann(size(raw, 1)), ...
    [1 size(raw,2)]);
end