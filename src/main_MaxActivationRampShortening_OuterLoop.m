clc;
close all;
clear all;

flag_useFlatActiveForceLengthCurve = 0;
flag_useConstantTendonStiffness    = 0;
flag_useLinearForceVelocityCurve   = 0;

flag_rampType = 1;
flag_useTendonDampingDampedEq = 1;   

% 0: matches the -15-to-15 degrees of ankle movement of Holzer et al.
% 1: matches the -20-to-35 degrees of ankle movement of Hauraix et al.
%    this also changes the angular velocites simulated to match Hauraix:
%    30,90,150,210,270,330

preloadHauraixReplication = 0.15; 
normalizedTendonDamping          = 0.05*(2/3);
normalizedTendonDampingConstant  = 0.05*(1/3);

flag_useHauraixVmax = 0;
maximumNormalizedFiberVelocity = 10;

ankleAngleMaxPlantarFlexion   = -17; %Holzer

scaleLceOpt = 1;

% if(flag_rampType==1)
%   scaleLceOpt                 = 1;
%   ankleAngleMaxPlantarFlexion = -17;
%   maximumNormalizedFiberVelocity = 6;
% end


standardMomentArm = 0.054;
smallMomentArm    = standardMomentArm/1.18;

standardTendonElasticity = 0.049; %Magnusson et al. 2001
highTendonElasticity     = 0.092; %Waugh et al. 2012 (adults results)

fractionOfFastTwitchFibers          = 0.5; 
% Fast: norm force of 0.22 at half vmax
% Slow: norm force of 0.15 at hal vmax
% From Ranatunga 1984 - experiments done on slow/fast twitch muscle in rats

flag_measurementSetting = 0;
%0: ankle angle
%1: fiber length

%Configure this run
 flag_rigidTendon         = 0;
 flag_standardTendon      = 1;
 flag_highlyElasticTendon = 1;

flag_standardMomentArm   = 1;
flag_smallMomentArm      = 0;

%Do not touch: used by main_MaxActivationRampShortening
flag_runRigidBench               = 0;
flag_runClassicElasticBench      = 0;
flag_runDampedFiberElasticBench  = 0;


%Rigid tendon, standard moment arm
if(flag_rigidTendon==1 && flag_standardMomentArm == 1)
  ankleAchillesTendonMomentArm        = standardMomentArm;
  tendonStrainAtOneNormForceOverride  =-1;
  flag_runRigidBench                  = 1;
  flag_runClassicElasticBench         = 0;
  flag_runDampedFiberElasticBench     = 0;
  flag_useTendonDamping =0;
  

  flag_simForceVelocityExpWithPreload = 0;
  preloadFraction = 0;
  main_MaxActivationRampShortening;
  
  flag_simForceVelocityExpWithPreload = 1;
  preloadFraction = 1;
  main_MaxActivationRampShortening;

  flag_simForceVelocityExpWithPreload = 1;
  preloadFraction = 0.5;
  main_MaxActivationRampShortening;

  flag_simForceVelocityExpWithPreload = 1;
  preloadFraction = 0.1;
  main_MaxActivationRampShortening;  
  
end

%Standard tendon, standard moment arm
if(flag_standardTendon==1 && flag_standardMomentArm == 1)
  ankleAchillesTendonMomentArm        = standardMomentArm;
  tendonStrainAtOneNormForceOverride  = standardTendonElasticity;
  flag_runRigidBench                  = 0;
  flag_runClassicElasticBench         = 0;
  flag_runDampedFiberElasticBench     = 1;
  
  flag_useTendonDamping = flag_useTendonDampingDampedEq;

  if(flag_rampType==1)
    flag_simForceVelocityExpWithPreload = 0;
    preloadFraction = preloadHauraixReplication;
    main_MaxActivationRampShortening;
    
  else
    flag_simForceVelocityExpWithPreload = 0;
    preloadFraction = 0;
    main_MaxActivationRampShortening;

    flag_simForceVelocityExpWithPreload = 1;
    preloadFraction = 1;
    main_MaxActivationRampShortening;

    flag_simForceVelocityExpWithPreload = 1;
    preloadFraction = 0.5;
    main_MaxActivationRampShortening;

    flag_simForceVelocityExpWithPreload = 1;
    preloadFraction = 0.1;
    main_MaxActivationRampShortening;    
  end
end

%More elastic tendon, standard moment arm
if(flag_highlyElasticTendon==1 && flag_standardMomentArm == 1)
  ankleAchillesTendonMomentArm        = standardMomentArm;
  tendonStrainAtOneNormForceOverride  = highTendonElasticity;
  flag_runRigidBench                  = 0;
  flag_runClassicElasticBench         = 0;
  flag_runDampedFiberElasticBench     = 1;

  flag_useTendonDamping = flag_useTendonDampingDampedEq;  
  

  if(flag_rampType==1)
    flag_simForceVelocityExpWithPreload = 0;
    preloadFraction = preloadHauraixReplication;
    main_MaxActivationRampShortening;
    
  else
    flag_simForceVelocityExpWithPreload = 0;
    preloadFraction = 0;
    main_MaxActivationRampShortening;

    flag_simForceVelocityExpWithPreload = 1;
    preloadFraction = 1;
    main_MaxActivationRampShortening;

    flag_simForceVelocityExpWithPreload = 1;
    preloadFraction = 0.5;
    main_MaxActivationRampShortening;

    flag_simForceVelocityExpWithPreload = 1;
    preloadFraction = 0.1;
    main_MaxActivationRampShortening;    
    
  end
  
end  
  

%Rigid tendon, shorter moment arm
if(flag_rigidTendon == 1 && flag_smallMomentArm == 1)
 ankleAchillesTendonMomentArm        = smallMomentArm;
  tendonStrainAtOneNormForceOverride  = -1;
  flag_runRigidBench                  = 1;
  flag_runClassicElasticBench         = 0;
  flag_runDampedFiberElasticBench     = 0;

  flag_useTendonDamping = 0;
  
  flag_simForceVelocityExpWithPreload = 0;
  preloadFraction = 0.;
  main_MaxActivationRampShortening;

  flag_simForceVelocityExpWithPreload = 1;
  preloadFraction = 1;
  main_MaxActivationRampShortening;
  
  flag_simForceVelocityExpWithPreload = 1;
  preloadFraction = 0.5;
  main_MaxActivationRampShortening;
  
  flag_simForceVelocityExpWithPreload = 1;
  preloadFraction = 0.1;
  main_MaxActivationRampShortening;  
end

%Standard tendon, shorter moment arm
if(flag_standardTendon == 1 && flag_smallMomentArm == 1)
  ankleAchillesTendonMomentArm        = smallMomentArm;
  tendonStrainAtOneNormForceOverride  = standardTendonElasticity;
  flag_runRigidBench                  = 0;
  flag_runClassicElasticBench         = 0;
  flag_runDampedFiberElasticBench     = 1;

  flag_useTendonDamping = flag_useTendonDampingDampedEq;  
  
  if(flag_rampType==1)
    flag_simForceVelocityExpWithPreload = 0;
    preloadFraction = preloadHauraixReplication;
    main_MaxActivationRampShortening;
    
  else
    flag_simForceVelocityExpWithPreload = 0;
    preloadFraction = 0;
    main_MaxActivationRampShortening;

    flag_simForceVelocityExpWithPreload = 1;
    preloadFraction = 1;
    main_MaxActivationRampShortening;

    flag_simForceVelocityExpWithPreload = 1;
    preloadFraction = 0.5;
    main_MaxActivationRampShortening;
    
    flag_simForceVelocityExpWithPreload = 1;
    preloadFraction = 0.1;
    main_MaxActivationRampShortening;
  end

end

%Elastic tendon, shorter moment arm
if(flag_highlyElasticTendon == 1 && flag_smallMomentArm == 1)
  ankleAchillesTendonMomentArm        = smallMomentArm;
  tendonStrainAtOneNormForceOverride  = highTendonElasticity;
  flag_runRigidBench                  = 0;
  flag_runClassicElasticBench         = 0;
  flag_runDampedFiberElasticBench     = 1;

  flag_useTendonDamping = flag_useTendonDampingDampedEq;  
  
  if(flag_rampType==1)
    flag_simForceVelocityExpWithPreload = 0;
    preloadFraction = preloadHauraixReplication;
    main_MaxActivationRampShortening;
    
  else
    flag_simForceVelocityExpWithPreload = 0;
    preloadFraction = 0;
    main_MaxActivationRampShortening;

    flag_simForceVelocityExpWithPreload = 1;
    preloadFraction = 1;
    main_MaxActivationRampShortening;

    flag_simForceVelocityExpWithPreload = 1;
    preloadFraction = 0.5;
    main_MaxActivationRampShortening;
    
    flag_simForceVelocityExpWithPreload = 1;
    preloadFraction = 0.1;
    main_MaxActivationRampShortening;    
  end

end