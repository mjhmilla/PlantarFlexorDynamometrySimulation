clc;
close all;
clear all;


standardMomentArm = 0.054;
smallMomentArm    = standardMomentArm/1.18;

standardTendonElasticity = 0.049;
highTendonElasticity     = 0.14;

fractionOfFastTwitchFibers          = 0.5; 
% Fast: norm force of 0.22 at half vmax
% Slow: norm force of 0.15 at hal vmax
% From Ranatunga 1984 - experiments done on slow/fast twitch muscle in rats

flag_measurementSetting = 0;
%0: ankle angle
%1: fiber length

%Standard tendon, standard moment arm

  ankleAchillesTendonMomentArm        = standardMomentArm;
  tendonStrainAtOneNormForceOverride  = standardTendonElasticity;

  flag_simForceVelocityExpWithPreload = 0;
  preloadFraction = 0;
  main_MaxActivationRampShortening;
  
  flag_simForceVelocityExpWithPreload = 1;
  preloadFraction = 1;
  main_MaxActivationRampShortening;

  flag_simForceVelocityExpWithPreload = 1;
  preloadFraction = 0.5;
  main_MaxActivationRampShortening;

%More elastic tendon, standard moment arm

  ankleAchillesTendonMomentArm        = standardMomentArm;
  tendonStrainAtOneNormForceOverride  = highTendonElasticity;

  flag_simForceVelocityExpWithPreload = 0;
  preloadFraction = 0;
  main_MaxActivationRampShortening;

  flag_simForceVelocityExpWithPreload = 1;
  preloadFraction = 1;
  main_MaxActivationRampShortening;

  flag_simForceVelocityExpWithPreload = 1;
  preloadFraction = 0.5;
  main_MaxActivationRampShortening;
  
  
%Standard tendon, shorter moment arm
  ankleAchillesTendonMomentArm        = smallMomentArm;
  tendonStrainAtOneNormForceOverride  = standardTendonElasticity;

  flag_simForceVelocityExpWithPreload = 0;
  preloadFraction = 0.;
  main_MaxActivationRampShortening;

  flag_simForceVelocityExpWithPreload = 1;
  preloadFraction = 1;
  main_MaxActivationRampShortening;
  
  flag_simForceVelocityExpWithPreload = 1;
  preloadFraction = 0.5;
  main_MaxActivationRampShortening;

%Elastic tendon, shorter moment arm
  ankleAchillesTendonMomentArm        = smallMomentArm;
  tendonStrainAtOneNormForceOverride  = highTendonElasticity;

  flag_simForceVelocityExpWithPreload = 0;
  preloadFraction = 0.;
  main_MaxActivationRampShortening;

  flag_simForceVelocityExpWithPreload = 1;
  preloadFraction = 1;
  main_MaxActivationRampShortening;
  
  flag_simForceVelocityExpWithPreload = 1;
  preloadFraction = 0.5;
  main_MaxActivationRampShortening;
