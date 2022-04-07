%% Entropy based cost value calculation
% calCostFun.m
function cost = calCostFun(raw, depthROI, maxDispOrders, dispCoeffs)
rawDispComp = compDisPhase(raw, maxDispOrders, dispCoeffs);
OCT = (abs(fft(rawDispComp))).^2;
OCTROI = OCT(depthROI(1):depthROI(2),:);
OCTNorm = OCTROI./sum(OCTROI(:));
entropy = -1*((OCTNorm).*log(OCTNorm));
cost = sum(entropy(:));
end