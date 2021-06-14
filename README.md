# PlantarFlexorDynamometrySimulation
A simulation of in-vivo dynamometry experiments on the plantar flexors carried out to estimate the force-velocity curve. 

This code is an extension of https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort

# Installation

1. Clone/download this repository to your machine.
2. Make sure these files are in the 'data' and 'dataKGR' directories:

- Isos_for_Matt.mat
- Isos_for_Matt.xlsx
- Own_Study.xlsx

If this data is missing, look at the supplementary data for the paper that accompanies this simulation

3. Check that everything works:
  - Run main_runAll.m, then run main_plotAll. This should save 3 pdf files to the "plots" folder:

    fig_ForceVelocity_Simulation_Vs_Experiment_lceOptScale_100p_vmax_10_ramp_Holzer_FastTwitch_50_Tdn_0p0_TendonDamping_0p_5p4cm_fixedAnkleAngle.pdf  

    fig_ForceVelocity_Simulation_Vs_Experiment_lceOptScale_100p_vmax_10_ramp_Holzer_FastTwitch_50_Tdn_9p2_TendonDamping_0p_5p4cm_fixedAnkleAngle

    fig_ForceVelocity_Simulation_Vs_Experiment_lceOptScale_100p_vmax_10_ramp_Holzer_FastTwitch_50_Tdn_4p9_TendonDamping_0p_5p4cm_fixedAnkleAngle

  - In addition a large number of MAT files will be written to the "data" folder:

# Overview

A number of functions have been added to so replicate plantar flexor dynamometry experiments using a gastrocnemius muscle and Achilles tendon:

1. main_MaxActivationRampShortening_OuterLoop.m and main_MaxActivationRampShortening.m
    This function simulates two different experimental protocols that are used to estimate the force-velocity curve of the plantar flexors in-vivo: one protocol has the participant fully activate their plantar flexors prior to shortening (preload), the other begins shortening as soon as the plantar flexors develop force (no preload). After simulation the results are written to mat files that appear in the data/ folder. As this function may need to be edited please see the Section "Guided Code-Tour: main_MaxActivationRampShortening.m" below for details.

2. main_runAll.m
    Runs simulations of the GM contracting during a dynamometry experiment under a variety of different preloads and tendon elasticities.

3. main_plotAll.m
    This function generates a detailed set of plots for each combination of model and preload.


# Code Tour: main_MaxActivationRampShortening.m

1. This file will loop through a grid of dynamometer angular velocities (for loop beginning on line 263). In each loop it will construct a muscle that has the architectural properties of the soleus muscle as documented by Arnold et al. 2010 and is consistent with the settings at the top of the file, activate it, and run it through the specified path function. The experiment will be repeated using angular velocities that start at 1 degree per second up to 211 degrees per second. The velocities tested are set using countMax near line 182 and mapping that to n the angular velocity in degrees per second near line 266.

2. If you look at the first if statement of the file "if(flag_outerLoopMode==0)" you will see all of the different input parameters that you can use to alter both the musculotendon model and the test. When the script is called with 'flag_outerLoopMode=1' then a parent script is used to set these parameters.

3. Read the comments below that begins with the title "Musculotendon Ramp-Shortening Function: Two different protocols" near line 470. Here you will find details on how I constructed a musclotendon path-length function that should be equivalent to the experiment you described to me. To alter the ramp simply change the parameters that appear below near the block of comments titled "Ramp-Shortening Function Parameters" near line 286.

 Note: I have given you a more advanced function to create the force velocity curve than is present in OpenSim: it fits the Bezier curve to Hill's hyperbola. For details read the comment in 'createFiberForceVelocityCurve2018.m' between lines 27-39. Now the muscle generates a bit more active force, though not much more. As far as I can tell the elasticity of the Achilles tendon has the stongest influnce on the active force created by the muscle during this ramp experiment.

4. If you want to see the tendon-force-length curve, the fiber-force-length curve, the active-fiber-force length curve, and the force velocity curve just set this flag to one (near line 274)

 flag_plotNormMuscleCurves            = 1;

 and re-run the script. This will give you the curves that are used to simulate the damped-fiber muscle model.

 Arnold, E. M., Ward, S. R., Lieber, R. L., & Delp, S. L. (2010). A model of the lower limb for analysis of human movement. Annals of biomedical engineering, 38(2), 269-279.

 Arnold, E. M., Hamner, S. R., Seth, A., Millard, M., & Delp, S. L. (2013). How muscle fiber lengths and velocities affect muscle force generation as humans walk and run at different speeds. Journal of Experimental Biology, jeb-075697.

<b>License</b>
```
Licensed under the MIT License (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License in this repository 
```
