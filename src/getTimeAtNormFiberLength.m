function measurementTime = getTimeAtNormFiberLength(simRecord, normFiberLength, count)

idxBestFit = 0;
errBestFit = Inf;
for z=1:1:size(simRecord.standardResults.normFiberLength,1)
  errLength = abs(simRecord.standardResults.normFiberLength(z,count) ...
                  - normFiberLength);
  if( errLength < errBestFit)
    idxBestFit=z;
    errBestFit=errLength;
  end                       
end

idxFitRange = [(idxBestFit-3):1:(idxBestFit+3)];

measurementTime = interp1( ...
  simRecord.standardResults.normFiberLength(idxFitRange,count),...
  simRecord.standardResults.time(idxFitRange,count),...
  normFiberLength);