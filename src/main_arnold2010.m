%%
% This function reads in Table 1 from Arnold 2010 into a muscle
% architecture structure that is compatible with the benchmark
%
% Arnold, E. M., Ward, S. R., Lieber, R. L., & Delp, S. L. (2010). 
% A model of the lower limb for analysis of human movement. 
% Annals of biomedical engineering, 38(2), 269-279.
%
% 
%%

clc;
close all;
clear all;

unitsMKSN = 1;
arnold2010LegArch = getArnold2010LegMuscleArchitecture(unitsMKSN);

index.soleus = getMuscleIndex(muscleAbbr,...
    arnold2010LegArch.abbrevation);

soleusArchitecture = [];
    soleusArchitecture.name                = arnold2010LegArch.names{index.soleus};
    soleusArchitecture.abbr                = arnold2010LegArch.abbrevation{index.soleus};
    soleusArchitecture.fiso                = arnold2010LegArch.peakForce(index.soleus);
    soleusArchitecture.optimalFiberLength  = arnold2010LegArch.optimalFiberLength(index.soleus);
    soleusArchitecture.maximumNormalizedFiberVelocity = maximumNormalizedFiberVelocity;
    soleusArchitecture.pennationAngle      = arnold2010LegArch.pennationAngle(index.soleus);
    soleusArchitecture.tendonSlackLength   = arnold2010LegArch.tendonSlackLength(index.soleus);
        
    minimumActiveFiberNormalizedLength = ...
        normMuscleCurves.activeForceLengthCurve.xEnd(1);
    
    
    minFiberKinematics = calcFixedWidthPennatedFiberMinimumLength(...
                minimumActiveFiberNormalizedLength,...
                maximumPennationAngle,...
                soleusArchitecture.optimalFiberLength,...
                soleusArchitecture.pennationAngle);
    
    soleusArchitecture.minimumFiberLength = ...
                                        minFiberKinematics.fiberLength;
                                    
    soleusArchitecture.minimumFiberLengthAlongTendon =...
                             minFiberKinematics.fiberLengthAlongTendon;
                         
    soleusArchitecture.pennationAngleAtMinumumFiberLength = ...
                                        minFiberKinematics.pennationAngle;