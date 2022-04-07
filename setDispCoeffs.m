%% Dispersion coefficient estimation
% setDispCoeffs.m
function dispCoeffs = setDispCoeffs(raw, depthROI, maxDispOrders, coeffRange)
dispCoeffs = zeros(1, maxDispOrders-1);
for i = 1:length(dispCoeffs)
    dispCoeffRng = [coeffRange, -1 * coeffRange];
    costs = zeros(size(dispCoeffRng));
    for j = 1:length(dispCoeffRng)
        dispCoeffs(i) = dispCoeffRng(j);
        costs(j) = calCostFun(raw, depthROI, maxDispOrders, dispCoeffs);
    end
    for k=1:20
        [~, idx] = sort(costs);
        dispcMin1 = dispCoeffRng(idx(1)); % index for the first minimum cost value
        dispcMin2 = dispCoeffRng(idx(2)); % index for the second minimum cost value
        dispcNew = (dispcMin1 + dispcMin2)/2;
        dispCoeffRng = [dispCoeffRng, dispcNew];
        dispCoeffs(i) = dispcNew;
        costNew = calCostFun(raw, depthROI, maxDispOrders, dispCoeffs);
        costs = [costs, costNew];
    end
    [~, argmin] = min(costs);
    dispCoeffs(i) = dispCoeffRng(argmin);
end
end