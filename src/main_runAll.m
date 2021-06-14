clc;
close all;
clear all;

flag_rigidTendon         = 1;
flag_standardTendon      = 1;
flag_highlyElasticTendon = 1;


flag_rampType                 = 0;
flag_useTendonDampingDampedEq = 1;   
main_MaxActivationRampShortening_OuterLoop;

flag_rampType                 = 0;
flag_useTendonDampingDampedEq = 0;   
main_MaxActivationRampShortening_OuterLoop;

flag_rampType                 = 1;
flag_useTendonDampingDampedEq = 1;   
main_MaxActivationRampShortening_OuterLoop;

flag_rampType                 = 1;
flag_useTendonDampingDampedEq = 0;   
main_MaxActivationRampShortening_OuterLoop;


