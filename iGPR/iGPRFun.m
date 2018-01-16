function [GP dataOut] = iGPRFun(data, options)
GP = []; 
dataOut = []; 
PtsNumThr = options.PtsNumThr;
epsilon = options.epsilon;
MinPts = options.MinPts;
tData = options.tData; 
tDist = options.tDist; 
if size(data, 2) <= PtsNumThr
    return; 
end
[~, ValidIdx] = Cluster(data, MinPts, epsilon, PtsNumThr);
if length(ValidIdx) <= options.PtsNumThr % point number is smaller than PtsNumThr
    return; 
end
GP = GPR_CurbNew_IV(data, ValidIdx, options, []);
dataOut =  GP.data(GP.ValidIdx, : )';   % data0 is valid road edge points.
end

