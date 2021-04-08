
function pathState = calcPathFunctionState(t,initialPauseTime,...
                                            timeSeries,...
                                            lengthSeries,...
                                            velocitySeries)


pathState = zeros(2,1);

dydt = 0;
y    = lengthSeries(1,1);

if(t > initialPauseTime && t <= (initialPauseTime + timeSeries(end)))
  dydt = interp1(timeSeries,velocitySeries, t-initialPauseTime,'linear','extrap');
  y    = interp1(timeSeries,  lengthSeries, t-initialPauseTime,'linear','extrap');
end
if(t>(initialPauseTime+timeSeries(end)))
  dydt = 0;
  y    = lengthSeries(end,1);
end

pathState(1) = dydt;
pathState(2) = y;
