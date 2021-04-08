clc;
close all;
clear all;

flag_rampType                 = 1;
flag_useTendonDampingDampedEq = 1;   
main_MaxActivationRampShortening_OuterLoop;


flag_rampType                 = 1;
flag_useTendonDampingDampedEq = 0;   
main_MaxActivationRampShortening_OuterLoop;

flag_rampType                 = 0;
flag_useTendonDampingDampedEq = 1;   
main_MaxActivationRampShortening_OuterLoop;


flag_rampType                 = 0;
flag_useTendonDampingDampedEq = 0;   
main_MaxActivationRampShortening_OuterLoop;
