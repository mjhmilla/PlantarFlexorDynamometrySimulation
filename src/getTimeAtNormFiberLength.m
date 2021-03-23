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



idxFitRange = [max((idxBestFit-3),1):1:min((idxBestFit+3), size(simRecord.standardResults.normFiberLength,1))];



measurementTime = interp1( ...
  simRecord.standardResults.normFiberLength(idxFitRange,count),...
  simRecord.standardResults.time(idxFitRange,count),...
  normFiberLength);