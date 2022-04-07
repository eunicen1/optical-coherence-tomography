%% Dispersion compensation process
% compDisPhase.m
function rawDispComp = compDisPhase(raw, maxDispOrders, dispCoeffs)
scanPts = size(raw,1);
linePerFrame = size(raw,2);
kLinear = linspace(-1,1,scanPts);
kaxis = repmat(kLinear', 1, linePerFrame);
rawDispComp = raw;
for i = 1:maxDispOrders-1
    rawDispComp = rawDispComp.*exp(1j.*(dispCoeffs(i)*(kaxis.^(i+1))));
end

end