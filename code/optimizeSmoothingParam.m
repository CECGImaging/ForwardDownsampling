function lambda = optimizeSmoothingParam(data, lapTranLap, interpMat, subInd, minLambda, maxLambda)

options = [];
% options = optimset('PlotFcns',@optimplotfval);
expon = fminbnd(@objFct, log10(minLambda), log10(maxLambda), options);
lambda = 10^expon;

function RMSE = objFct(expon)
    lambda = 10^expon;
    S = speye(size(lapTranLap)) + lambda*lapTranLap;
    P = ichol(S, struct('diagcomp',1e-1));
    [data_smoothed,~] = pcg(S, data, 1e-6, 1e4, P, P', data);
    RMSE = rms(data - interpMat*data_smoothed(subInd));
end

end