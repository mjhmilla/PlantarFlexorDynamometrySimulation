function [fiberVelocityError] = calcFiberTendonVelocityDifference(...
  vel, dataPathVelocity, dataFiberVelocityAlongTendon,dataTendonVelocity)

tendonVelocity = interp1(dataPathVelocity,...
                         dataTendonVelocity,vel);
fiberVelocityAlongTendon = interp1(dataPathVelocity,...
                                   dataFiberVelocityAlongTendon,vel);
fiberVelocityError = fiberVelocityAlongTendon-tendonVelocity;