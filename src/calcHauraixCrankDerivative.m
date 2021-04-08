function der = calcHauraixCrankDerivative(t, x, dataHauraix)

der= interp1(dataHauraix.x,dataHauraix.y,x,'linear','extrap');
