# PlantarFlexorDynamometrySimulation
A simulation of in-vivo dynamometry experiments on the plantar flexors carried out to estimate the force-velocity curve. 

This code is an extension of https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort

# Installation

1. Clone/download this repository to your machine.
2. Put these files into the 'data' folder (which comes empty):

- Isos_for_Matt.mat
- Isos_for_Matt.xlsx
- Own_Study.xlsx

Contact Denis Holzer for this data

3. Check that everything works:
  - Run main_ActivationDynamicsTuning.m: this should save a pdf file to the "plots" folder:

    fig_IsometricNormMusculotendonForce_ExpAndSim.pdf  

  - Run main_MaxActivationRampShortening.m: this will save 4 files to the "data" folder:

    simFv_gaslat_preload_0_dampedFiberElasticTendon.mat
    simFv_gaslat_preload_0_rigidTendon.mat
    simFv_gaslat_preload_1_dampedFiberElasticTendon.mat
    simFv_gaslat_preload_1_rigidTendon.mat

  - Run main_PlotMaxActivationRampShortening.m: this should save a pdf file to the "plots" folder:

    fig_ForceVelocity_Simulation_Vs_Experiment.pdf

# Overview

There are three functions that have been added to simulate the in-vivo 
force-velocity dynamometry experiments that are being analyzed by Denis Holzer 
and Wolfgang Sieberl:

1. main_ActivationDynamicsTuning.m
    Compares the normalized force recorded during isometric trials to the isometric force developed by a simulated muscle with the activation dynamics equations and parameters that appear in Thelen. The time-course of the active force developed by the gaslat is dominated entirely by the force-velocity curve and the long elastic Achilles tendon. The experimental data and the simulated time course agree well enough that it is quite reasonable to use the standard model of activation dynamics to simulate these experiments.

2. main_MaxActivationRampShortening.m
    This function simulates two different experimental protocols that are used to estimate the force-velocity curve of the plantar flexors in-vivo: one protocol has the participant fully activate their plantar flexors prior to shortening (preload), the other begins shortening as soon as the plantar flexors develop force (no preload). After simulation the results are written to mat files that appear in the data/ folder. As this function may need to be edited please see the Section "Guided Code-Tour: main_MaxActivationRampShortening.m" below for details.

3. main_PlotMaxActivationRampShortening.m
    This function plots the simulation results.

Thelen DG. Adjustment of muscle mechanics model parameters to 
simulate dynamic contractions in older adults. J. Biomech. Eng. 
2003 Feb 1;125(1):70-7.


# Code Tour: main_MaxActivationRampShortening.m

1. Read the comments below that begins with the title "Musculotendon Ramp-Shortening Function: Two different protocols" near line 259. Here you will find details on how I constructed a musclotendon path-length function that should be equivalent to the experiment you described to me. To alter the ramp simply change the parameters that appear below near the block of comments titled "Ramp-Shortening Function Parameters" near line 117.

2. If the ramp function I have constructed looks fine to you, run this file. It will construct a muscle that has the architectural properties of the soleus muscle as documented by Arnold et al. 2010, maximally activate it, and run it through the specified path function. I've set the initial angular velocity of the ramp to be 300 degrees per second. At this velocity the force generated by the musculotendon does not go to zero. From this simulation 3 figures will be produced. Note the 'DFE' acronym stands for 'damped fiber equilibrium' which is type of the muscle model I recommend you simulate.

  - Fig 1: Basic properties. The normalized fiber force along the tendon is probably most useful to you in subplot 1,2.
  - Fig 2: Potential energy and work. Here you can see how the tendon/fiber are storing potential energy and how much work is being done by the fiber, the light damping in the fiber, and by the muscle in total. Note: Subplot 1,3 should be close to zero, and it is: the title is likely obscuring the order of magnitude which is around 10^-4.
  - Fig 3: Power of the passive and active components.

3. Try setting the shortening velocity to be very high - higher than is possible with a Biodex machine, say 1000 deg/sec, by setting rampAngularVelocity = 1000; Now when you re-run the script you will see the force developed by the muscle does go to zero. It may not in your experimental subjects because of the parameters used for the elasticity of the Achilles tendon and the force velocity curve. As noted by Arnold et al. 2013 an EMG driven model of walking and running produced ankle torques that best matched experimental data when the tendon stretched by 10% at one isometric force (eIso = 0.10). By default eIso is set to 4.9%. To change this

 - Open 'createDefaultNormalizedMuscleCurves.m'
 - Go to line 131 where the tendon curve is created
 - Change the eIso variable from 0.048 to 0.10.
 - Re-run the script.

Now the musculotendon force goes to zero but only at 15 degrees of plantar flexion. At an ankle angle of 0 (at a time of 0.215 s) the muscle is developing a normalized force of 0.2 - which is huge! This confirms the idea, to me, that your force-velocity curve is not going to zero due to series elasticity.

4. There are at least 3 other parameters that could differ between your subjects and the model which will affect the amount of force the muscle develops during your experiment:

 a. maximumNormalizedFiberVelocity : set to 10 by default. But this could be slower or faster
 b. Force-Velocity Curve: you can change the curve to be consistent with a slow twitch muscle or a fast twitch muscle by

 - Open 'createDefaultNormalizedMuscleCurves.m'
 - Go to line 91 & read the comment below fvAtHalfVMax
 - Set fvAtHalfVMax to 0.22 - now it has the normalized force velocity
   curve of a fast twitch muscle
 - Re-run the script

 Note: I have given you a more advanced function to create the force velocity curve than is present in OpenSim. For details read the comment in 'createFiberForceVelocityCurve2018.m' between lines 27-39. Now the muscle generates a bit more active force, though not much more. As far as I can tell the elasticity of the Achilles tendon has the stongest influnce on the active force created by the muscle during this ramp experiment.

5. If you want to see the tendon-force-length curve, the fiber-force-length curve, the active-fiber-force length curve, and the force velocity curve just set this flag to one (line 116)

 flag_plotNormMuscleCurves            = 1;

 and re-run the script. This will give you the curves that are used to simulate the damped-fiber muscle model.


 Arnold, E. M., Ward, S. R., Lieber, R. L., & Delp, S. L. (2010). A model of the lower limb for analysis of human movement. Annals of biomedical engineering, 38(2), 269-279.

 Arnold, E. M., Hamner, S. R., Seth, A., Millard, M., & Delp, S. L. (2013). How muscle fiber lengths and velocities affect muscle force generation as humans walk and run at different speeds. Journal of Experimental Biology, jeb-075697.

<b>License</b>
```
Licensed under the MIT License (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License in this repository 
```
