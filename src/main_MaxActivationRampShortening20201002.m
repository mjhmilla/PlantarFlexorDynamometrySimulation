
%%
% Denis: To generate a force vs. crank velocity with different tendon elasticities
%
% 1. Set the initial 'tendonStrainAtOneNormForceOverride' - e.g. 0.049
% 2. Set flag_updateExistingPlots = 0;
% 3. Run the script.
% 4. Update  'tendonStrainAtOneNormForceOverride' to be the next value -
%    e.g. 0.12
% 5. Set flag_updateExistingPlots = 1;
% 6. Run the script.
% 7. Repeat steps 4 and 6 as needed.
%%
flag_updateExistingPlots = 0;

for i = 1
%legendString = {''};
if(flag_updateExistingPlots == 0)
  clc;
  close all;
  clear all;
  flag_updateExistingPlots = 0;
  
  path = ['../data/Own_Study.xlsx'];
    Norm_Muscle_Force_Individual = xlsread(path,'Muscle_Force_Individual');
    Vicon_angular_Velocity = xlsread(path,'Vicon_angular_Velocity');
    FSZ_Velocity_Individual = xlsread(path,'Velocity_Individual');
    Isomed_angular_Velocity = xlsread(path,'Isomed_angular_Velocity');
    
    Hauraix_velocities = xlsread(path,'Hauraix_velocities');
    Hauraix_angular_velo = Hauraix_velocities(:,2);
    Hauraix_FSZ_velo = Hauraix_velocities(:,1);
    
    Chino_Joint_Velo = [0 51.642 98.657 134.254 188.657];
    Chino_Joint_Velo_STD = [0 5.038 9.739 20.821 16.119];

    Chino_force = [1 .868 .768 .697 .603];
    Chino_force_STD = [0 .143 .203 .242 .178];

    Chino_Fascicle_Velo = [80.168 54.246 42.107 22.128 0];
    Chino_Fascicle_Velo = flip(Chino_Fascicle_Velo);
    Chino_Fascicle_Velo_STD = [23.52 19.221 12.266 5.564 0];
    Chino_Fascicle_Velo_STD = flip(Chino_Fascicle_Velo_STD);

    Hauraix_Joint_Velo = [0 30 90 150 210 290 330];
    Hauraix_Force = [431.76 380.33 262.52 198.54 159.07 139.33 129.77];
    Hauraix_Force_norm = Hauraix_Force./Hauraix_Force(1);
    Hauraix_Fascicle_Velo = [0 10.61 46.06 72.85 98.53 120.58 131.19];

    figure(101)
    hold on
    m = errorbar(mean(FSZ_Velocity_Individual'),mean(-Vicon_angular_Velocity'),...
        std(-Vicon_angular_Velocity'),std(-Vicon_angular_Velocity'),...
        std(FSZ_Velocity_Individual'),std(FSZ_Velocity_Individual'),'ko');  
    m.Marker = 'o';
    m.Color = 'k';
%    m.LineStyle = '';    
    
     [p] = polyfit(mean(FSZ_Velocity_Individual'),mean(-Vicon_angular_Velocity'),1);
     xxx = linspace(min(mean(FSZ_Velocity_Individual')),max(mean(FSZ_Velocity_Individual')));
     yyy = polyval(p,xxx);
     
    ft = fittype( 'poly1' );
    [fit_A,a] = fit(mean(FSZ_Velocity_Individual')',mean(-Vicon_angular_Velocity')',ft);
    txt = num2str(a.rsquare);
    
    
    plot(xxx,yyy,'k:')
    txt = {'R^2 of linear fit:',txt};
    text(10,160,txt,'FontSize',14)
    
    xlabel('Fascicle shortening speed [mm/s]','FontSize',14)
    ylabel('Ankle joint angular velocity [�/s]','FontSize',14)   

    

    figure(1)
    hold on
    k = errorbar(Chino_Fascicle_Velo,Chino_force,...
        Chino_force_STD,Chino_force_STD,...
        Chino_Fascicle_Velo_STD,Chino_Fascicle_Velo_STD)

    k.Marker = '^';
    k.Color = 'r';
    k.LineStyle = ':';

     plot(Hauraix_Fascicle_Velo,Hauraix_Force_norm,':sb')

    hold on
    e = errorbar(mean(FSZ_Velocity_Individual'),mean(Norm_Muscle_Force_Individual'),...
        std(Norm_Muscle_Force_Individual'),std(Norm_Muscle_Force_Individual'),...
        std(FSZ_Velocity_Individual'),std(FSZ_Velocity_Individual'))

    e.Marker = 'o';
    e.Color = 'k';
    e.LineStyle = ':';
    
    ylim([0 1.1])
    legend('Chino et al., 2008','Hauraix et al., 2015','current study')
    %legend('Hauraix et al., 2015','current study')
    xlabel('Fascicle velocity [mm/s]')
    ylabel('normalized fascicle force along tendon [%]')
  
    figure(2)
    hold on
    k = errorbar(Chino_Joint_Velo,Chino_force,...
        Chino_force_STD,Chino_force_STD,...
        Chino_Joint_Velo_STD,Chino_Joint_Velo_STD)

    k.Marker = '^';
    k.Color = 'r';
    k.LineStyle = ':';
% 
     plot(Hauraix_Joint_Velo,Hauraix_Force_norm,':sb')

    hold on
    e = errorbar(mean(-Vicon_angular_Velocity'),mean(Norm_Muscle_Force_Individual'),...
        std(Norm_Muscle_Force_Individual'),std(Norm_Muscle_Force_Individual'),...
        std(-Vicon_angular_Velocity'),std(-Vicon_angular_Velocity'))

    e.Marker = 'o';
    e.Color = 'k';
    e.LineStyle = ':';    
    
    ylim([0 1.1])
    legend('Chino et al., 2008','Hauraix et al., 2015','current study')
    %legend('Hauraix et al., 2015','current study')
    xlabel('Angular velocity [�/s]')
    ylabel('normalized fascicle force along tendon [%]')
    
    Hauraix13_Fascicle_Velo = [56 85 113 131 153 168];
    Hauraix13_Joint_Velo = [30 90 150 210 270 330];
    
    figure(12354323)
    plot(Hauraix_Fascicle_Velo,Hauraix_Joint_Velo)
    hold on
    plot(Chino_Fascicle_Velo,Chino_Joint_Velo)
    plot(Hauraix13_Fascicle_Velo,Hauraix13_Joint_Velo)
end
end
tendonStrainAtOneNormForceOverride = 0.1;
muscle = 'gaslat';
initialHoldTime = 1;
shiftFiberForceLengthCurve = 0.;
%%
% Denis: You can change the tendon strain when one norm force is applied by 
%        changing the above parameter. Here are some values to try
%
% 0.049     : The default value used in OpenSim 
% 0.10-0.15 : A range of values that is typical for the Achilles tendon 
%             for EMG-driven models of running 
%            
%%


flag_runRigidBench               = 1;
flag_runClassicElasticBench      = 0;
flag_runDampedFiberElasticBench  = 1;

if(flag_runRigidBench==1 && flag_updateExistingPlots==0)
  figRigidTendonBasic  = figure;
  figRigidTendonEnergy = figure;
  figRigidTendonPower  = figure;
end

if(flag_runClassicElasticBench==1 && flag_updateExistingPlots==0)
  figClassicElasticBasic  = figure;
  figClassicElasticEnergy = figure;
  figClassicElasticPower  = figure;
end

if(flag_runDampedFiberElasticBench==1 && flag_updateExistingPlots==0)
  figDampedFiberElasticBasic  = figure;
  figDampedFiberElasticEnergy = figure;
  figDampedFiberElasticPower  = figure;
end
%%
% Dear Wolfgang,
%
% To start, take a guided tour through the code:
%
%
% 1. Read the comments below that begins with the title "Musculotendon
% Ramp-Shortening Function" near line 262. Here you will find details on
% how I constructed a musclotendon path-length function that should be
% equivalent to the experiment you described to me.
%
% To alter the ramp simply change the parameters that appear below near
% the block of comments titled "Ramp-Shortening Function Parameters" near
% line 124.
%
% 2. If the ramp function I have constructed looks fine to you, run this
% file. It will construct a muscle that has the architectural properties
% of the soleus muscle as documented by Arnold et al. 2010, maximally
% activate it, and run it through the specified path function. I've set
% the initial angular velocity of the ramp to be 300 degrees per second.
% At this velocity the force generated by the musculotendon does not go
% to zero.
%
% From this simulation 3 figures will be produced. Note the 'DFE' acronym
% stands for 'damped fiber equilibrium' which is type of the muscle model I
% recommend you simulate.
%
% Fig 1: Basic properties. The normalized fiber force along the tendon is
%        probably most useful to you in subplot 1,2.
%
% Fig 2: Potential energy and work. Here you can see how the tendon/fiber
%        are storing potential energy and how much work is being done
%        by the fiber, the light damping in the fiber, and by the
%        muscle in total.
%
%        Note: Subplot 1,3 should be close to zero, and it is: the title
%              is likely obscuring the order of magnitude which is around
%              10^-4.
% Fig 3: Power of the passive and active components.
%
% 3. Try setting the shortening velocity to be very high - higher than
% is possible with a Biodex machine, say 1000 deg/sec, by setting
%
% rampAngularVelocity = 1000;
%
% Now when you re-run the script you will see the force developed by the
% muscle does go to zero. It may not in your experimental subjects because
% of the parameters used for the elasticity of the Achilles tendon and the
% force velocity curve. As noted by Arnold et al. 2013 an EMG driven model
% of walking and running produced ankle torques that best matched
% experimental data when the tendon stretched by 10% at one isometric
% force (eIso = 0.10). By default eIso is set to 4.9%. To change this
%
% - Open 'createDefaultNormalizedMuscleCurves.m'
% - Go to line 131 where the tendon curve is created
% - Change the eIso variable from 0.048 to 0.10.
% - Re-run the script.
%
% Now the musculotendon force goes to zero but only at 15 degrees of
% plantar flexion. At an ankle angle of 0 (at a time of 0.215 s) the muscle
% is developing a normalized force of 0.2 - which is huge! This confirms
% the idea, to me, that your force-velocity curve is not going to zero due
% to series elasticity.
%
% 4. There are at least 3 other parameters that could differ between your
% subjects and the model which will affect the amount of force the muscle
% develops during your experiment:
%
% 4a. maximumNormalizedFiberVelocity : set to 10 by default. But this could
%       be slower or faster
% 4b. Force-Velocity Curve: you can change the curve to be consistent with
%     a slow twitch muscle or a fast twitch muscle by
%
% - Open 'createDefaultNormalizedMuscleCurves.m'
% - Go to line 91 & read the comment below fvAtHalfVMax
% - Set fvAtHalfVMax to 0.22 - now it has the normalized force velocity
%   curve of a fast twitch muscle
% - Re-run the script
%
% *Note: I have given you a more advanced function to create the force
% velocity curve than is present in OpenSim. For details read the comment
% in 'createFiberForceVelocityCurve2018.m' between lines 27-39.
%
% Now the muscle generates a bit more active force, though not much more.
% As far as I can tell the elasticity of the Achilles tendon has the
% stongest influnce on the active force created by the muscle during this
% ramp experiment.
%
% 4. If you want to see the tendon-force-length curve, the
% fiber-force-length curve, the active-fiber-force length curve, and the
% force velocity curve just set this flag to one (line 116)
%
% flag_plotNormMuscleCurves            = 1;
%
% and re-run the script. This will give you the curves that are used to
% simulate the damped-fiber muscle model.
%
%
% Arnold, E. M., Ward, S. R., Lieber, R. L., & Delp, S. L. (2010). A model
% of the lower limb for analysis of human movement. Annals of biomedical
% engineering, 38(2), 269-279.
%
% Arnold, E. M., Hamner, S. R., Seth, A., Millard, M., & Delp, S. L.
% (2013). How muscle fiber lengths and velocities affect muscle force
% generation as humans walk and run at different speeds. Journal of
% Experimental Biology, jeb-075697.
%%
% Loop with changing velocities
count=0;
countMax=21;

colorStart = [0,0,1];
colorEnd   = [1,0,0];
colorRecord = zeros(countMax,3);

ankleAngularVelocity = zeros(countMax,1);
fiberVelocity        = zeros(countMax,1);
normFiberForceAlongTendonIsometric = 0; %Isometric force generated by the muscle
                                        %mid-ramp
normFiberForceAlongTendonIsometric_rt = 0;

npts = 100;

%Elastic tendon data
ankle_angularVelocity       = zeros(npts,countMax);
fiber_force                 = zeros(npts,countMax);
fiber_length                = zeros(npts,countMax);
fiber_velocity              = zeros(npts,countMax);
fiber_velocityAlongTendon   = zeros(npts,countMax);
path                        = zeros(npts,countMax);
pennation_angle             = zeros(npts,countMax);
sim_time                    = zeros(npts,countMax);
measurement_time            = zeros(npts,1);

%Rigid tendon data
ankle_angularVelocity_rt       = zeros(npts,countMax);
fiber_force_rt                 = zeros(npts,countMax);
fiber_length_rt                = zeros(npts,countMax);
fiber_velocity_rt              = zeros(npts,countMax);
fiber_velocityAlongTendon_rt   = zeros(npts,countMax);
path_rt                        = zeros(npts,countMax);
pennation_angle_rt             = zeros(npts,countMax);
sim_time_rt                    = zeros(npts,countMax);
measurement_time_rt            = zeros(npts,1);


for count=1:1:countMax
    n = (count-1)*10+1;

    %%
    flag_useArnold2010SoleusArchitecture = 1;
    muscleAbbrArnold2010                 = muscle;
    
    flag_plotNormMuscleCurves            = 0;
    flag_updateNormMuscleCurves          = 0;
    if(count == 1)
      flag_updateNormMuscleCurves=1;
    end
    
    
    trialColor =   colorEnd.*((count-1)/(countMax-1)) ...
               + colorStart.*(1-((count-1)/(countMax-1)));
    colorRecord(count,:)=trialColor;
             
    %%
    % Ramp-Shortening Function Parameters
    %%
    
    %Path function inputs
    initialHoldTime              = initialHoldTime;   %Initial time [s] at the starting angle
    rampStartAngle               = -15;   %degrees, -ve: dorsiflexion
    rampEndAngle                 =  15;    %degrees
    rampAngularVelocity          =  n;  %degrees per second.
    ankleAngularVelocity(count)  = n;
    %Parameters from the literature
    ankleAchillesTendonMomentArm = 0.054; %moment arm in m
    %from Rugg et al.
    ankleAngleMaxPlantarFlexion  = -23;   %in degrees
    %From Anderson et al.
    
    %%
    %Parameters that apply to all muscles
    %%
    maximumNormalizedFiberVelocity = 10; %in units of norm fiber lengths/second
    % Note:
    %  A max. normalized fiber velocity of 10 is a starting point. Mammalian
    %  slow twitch fibers can have substantially slower maximum shortening
    %  velocities and fast twitch fibers can have (slightly) faster shortening
    %  velocities.
    %
    maximumPennationAngle          = 89*(pi/180); %if we go to 90 the
    %classic formulation goes
    %singular.
    
    
    
    
    %%
    %Get a copy of the default muscle curves
    %%
    muscleAbbr  = [];
    if(flag_useArnold2010SoleusArchitecture == 1)
        muscleAbbr = muscleAbbrArnold2010;
    else
        muscleAbbr = 'compBench';
    end
    
    
    normMuscleCurves = ...
        createDefaultNormalizedMuscleCurves(muscleAbbr,...
        tendonStrainAtOneNormForceOverride,...
        shiftFiberForceLengthCurve,...
        flag_updateNormMuscleCurves,...
        flag_plotNormMuscleCurves);
    
    
    %%
    %Get a muscle and extract out its architecture information
    %%
    
    muscleName  = [];
    fiso        = [];
    lceOpt      = [];
    alphaOpt    = [];
    ltSlk       = [];
    
    
    if(flag_useArnold2010SoleusArchitecture ==1)
        unitsMKSN = 1;
        arnold2010LegArch = getArnold2010LegMuscleArchitecture(unitsMKSN);
        
        idx =  getArnold2010MuscleIndex(muscleAbbrArnold2010,...
            arnold2010LegArch.abbrevation);
        
        muscleName  = arnold2010LegArch.names{idx};
        fiso        = arnold2010LegArch.peakForce(idx);
        lceOpt      = arnold2010LegArch.optimalFiberLength(idx);
        alphaOpt    = arnold2010LegArch.pennationAngle(idx);
        ltSlk       = arnold2010LegArch.tendonSlackLength(idx);        
    else
        muscleName                     = 'compBenchMillard2010';
        fiso    = 1;
        lceOpt  = 0.02;
        alphaOpt= 30*(pi/180);
        ltSlk   = 0.20;
    end
    
    %%
    % Wolfgang: If you want to alter the achitectural properties of the
    % musculotendon you can do so here. Simply set the parameters
    %
    % fiso      : max. isometric active force (N)
    % lceOpt    : optimal fiber length (m)
    % alphaOpt  : pennation angle at the optimal fiber length (rad)
    % ltSlk     : tendon slack length
    %
    % to what you want.
    %
    %%
    
    
    muscleArch = [];
    muscleArch.name                = muscleName;
    muscleArch.abbr                = muscleAbbr;
    muscleArch.fiso                = fiso;
    muscleArch.optimalFiberLength  = lceOpt;
    muscleArch.maximumNormalizedFiberVelocity = ...
        maximumNormalizedFiberVelocity;
    muscleArch.pennationAngle      = alphaOpt;
    muscleArch.tendonSlackLength   = ltSlk;
    
    minimumActiveFiberNormalizedLength = ...
        normMuscleCurves.activeForceLengthCurve.xEnd(1);
    
    minFiberKinematics = calcFixedWidthPennatedFiberMinimumLength(...
        minimumActiveFiberNormalizedLength,...
        maximumPennationAngle,...
        muscleArch.optimalFiberLength,...
        muscleArch.pennationAngle);
    
    muscleArch.minimumFiberLength = ...
        minFiberKinematics.minimumFiberLength;
    
    muscleArch.minimumFiberLengthAlongTendon =...
        minFiberKinematics.minimumFiberLengthAlongTendon;
    
    muscleArch.pennationAngleAtMinumumFiberLength = ...
        minFiberKinematics.pennationAngleAtMinimumFiberLength;
    
    
    
    %%
    %Run the computational benchmark described in:
    %  Millard, M., Uchida, T., Seth, A., & Delp, S. L. (2013).
    %    Flexing computational muscle: modeling and simulation of
    %    musculotendon dynamics. Journal of biomechanical engineering,
    %    135(2), 021005.
    %%
    benchConfig.npts   = npts;
    benchConfig.nSim   = 1;
    benchConfig.relTol = 1e-5;
    benchConfig.absTol = 1e-5;
    benchConfig.initialActivationState = 0;%[0:0.1:1]';
    benchConfig.color = trialColor;
    
    
    %%
    % Musculotendon Ramp-Shortening Function
    %
    %  1. Start at a static hold of 15 degrees of dorsi flexion
    %  2. Proceed at a constant shortening velocity
    %
    % We must map this description from being at the kinematic level of the
    % ankle to being in terms of the path length of the muscle. To do this we
    % must make use of the following facts:
    %
    % According to the measurements of Anderson et al. the plantar flexors
    % deliver their peak ankle torque at a dorsi flexion angle of -23 degrees
    % (from the young male category of Anderson et al.'s measurements.)
    %
    % According to Rugg et al. the moment arm of the Achilles tendon about the
    % ankle joint (from an MRI study) between plantar and dorsiflexion is
    % 5.4 cm +/- 0.3 cm. Note in Rugg et al. the border between plantar and
    % dorsiflexion occurrs at an ankle ankle of 1.92 radians using the the
    % anatomical axis they defined.
    %
    % From here on in we are going to assume that the moment arm length
    % reflects the subjects used in the the present study, and that the moment
    % arm length stays constant over the angles of interest. The first
    % assumption is probably fine if the subjects are young males and are not
    % elite sprinters/jumpers. The second assumption is reasonable according
    % to Rugg et al.'s data as the moment arm of the Achilles tendon varies
    % from by only +/- 3mm for a change in +/- 15 degrees about an anatomical
    % flat foot position.
    %
    % I. Starting Musculotendon Path Length
    %
    % Thus at the starting angle of 15 degrees of dorsi flexion the
    % musculotendon path length will be 0.754 cm shorter than it is at the
    % optimal length:  
    %
    % delta Ls = (delta theta)*r
    %         = (15-23)*(pi/180)*5.4 cm
    %         = -0.754 cm
    %
    % This gives us the starting musculotendon path length: 0.754 cm shorter
    % than optimal. Note 'delta Ls' stands for change in start length
    %
    % II. Change Musculotendon Path Length
    %
    % As the ankle rotates from 15 degrees of dorsi flexion to 15 degrees of
    % plantar flexion the path changes in length by:
    %
    % delta Lc = (delta theta)*r
    %          = (30)*(pi/180)*5.4cm
    %          = 2.8274 cm
    %
    % where delta Lc means the change in path length during the ramp.
    %
    % III. Time required to complete the ramp
    %
    % To compute the time required to complete the ramp we can again make use
    % of the moment arm:
    %
    % dt = 30 deg/omega
    %
    % To map between angular velocity and path shortening velocity we can
    % simply use the moment arm
    %
    % dl/dt = omega*r
    %
    % For example:
    %
    % @100 deg/sec
    % dl/dt = 100*(pi/180)*5.4
    %       = 9.42 cm/s
    % dt    = 30/100
    %       = 0.3s
    %
    % @500 deg/sec
    % dl/dt = 500*(pi/180)*5.4
    %       = 47.12 cm/s
    % dt    = 30/500
    %       = 0.06s
    %
    %  D. Anderson, M. Madigan, M. Nussbaum, Maximum voluntary
    %  joint torque as a function of joint angle and angular velocity:
    %  model development and application to the lower limb, Journal
    %  of Biomechanics 40 (14) (2007) 3105–3113.
    %
    %  Rugg, S. G., Gregor, R. J., Mandelbaum, B. R., & Chiu, L. (1990).
    %  In vivo moment arm calculations at the ankle using magnetic resonance
    %  imaging (MRI). Journal of biomechanics, 23(5), 495-497.
    %
    %%
    lceOpt   = muscleArch.optimalFiberLength;
    alphaOpt = muscleArch.pennationAngle;
    ltSlk    = muscleArch.tendonSlackLength;
    fiso     = muscleArch.fiso;
    
    
    
    %%
    % I. Starting Musculotendon Path Length
    %%
    
    %Solve for the length of the tendon when the fiber is at its optimal
    %length, optimal pennation angle, and is fully activated
    ftNOpt     = 1*cos(alphaOpt);
    ltNAtFiso  = normMuscleCurves.tendonForceLengthCurve.xEnd(2);
    ftNerr     = 1;
    tol        = 1e-6;
    iter       = 1;
    iterMax    =100;
    
    
    while(abs(ftNerr) > tol && iter < iterMax)
        ftNerr = calcBezierYFcnXDerivative(ltNAtFiso,...
            normMuscleCurves.tendonForceLengthCurve,0) ...
            - ftNOpt;
        DftNerr_DltN = calcBezierYFcnXDerivative(ltNAtFiso,...
            normMuscleCurves.tendonForceLengthCurve,1);
        dltN      = -ftNerr/DftNerr_DltN;
        ltNAtFiso = ltNAtFiso+dltN;
        iter=iter+1;
    end
    
    lpOpt = lceOpt*cos(alphaOpt) + ltSlk*ltNAtFiso;
    
    lpDeltaStart = -(rampStartAngle-ankleAngleMaxPlantarFlexion ...
        )*(pi/180)*ankleAchillesTendonMomentArm;
    
    lpRampStart = lpOpt + lpDeltaStart;
    
    %%
    % II. Change in Path Length
    %%
    
    lpDelta   = -(rampEndAngle-rampStartAngle) ...
        *(pi/180)*ankleAchillesTendonMomentArm;
    lpRampEnd = lpRampStart + lpDelta;
    
    lpRampMid = lpRampStart + 0.5*lpDelta;
    
    %%
    % III. Time required to complete the ramp
    %%
    
    rampTime = (rampEndAngle-rampStartAngle)/rampAngularVelocity;
    
    
    
    totalSimTime = 1;
    minHoldTime  = 0.25;
    if( (totalSimTime-(initialHoldTime+rampTime)) < minHoldTime)
      totalSimTime = (initialHoldTime+rampTime)+minHoldTime;
    end
    
    measurementTime = initialHoldTime + 0.5*rampTime;
    %%
    % Construct the path function
    %%
    
    pathFcn = @(t)calcRampFunctionState(t,initialHoldTime,...
        lpRampStart,lpRampEnd,rampTime);
    
    benchConfig.pathFcn = pathFcn;
    benchConfig.tspan   = [0, (totalSimTime)];
    
    %
    % Thelen DG. Adjustment of muscle mechanics model parameters to 
    % simulate dynamic contractions in older adults. J. Biomech. Eng. 
    % 2003 Feb 1;125(1):70-7.
    %
    tact   = 0.015;
    tdeact = 0.050;
    tmin   = 0.;
    
    ton = 0.;
    toff = totalSimTime;
    
    excitationFcn = @(argTime)calcStepFunction(argTime, ton, toff);    
    activationFcn = @(argEx, argAct)calcFirstOrderActivationDerivative(...
                                    argEx, argAct, tact,tdeact,tmin);
    
    benchConfig.excitationFcn = excitationFcn;
    benchConfig.activationFcn = activationFcn;
    
    %%=========================================================================
    %Rigid Tendon Model Benchmark
    %%=========================================================================
    
    if(flag_runRigidBench == 1)
        disp('Rigid-tendon model: ramp-shortening, constant activation');
        %figRigidTendonBasic  = figure;
        %figRigidTendonEnergy = figure;
        %figRigidTendonPower  = figure;
        
        rigidConfig.useFiberDamping  = 0;
        rigidConfig.useElasticTendon = 0;
        rigidConfig.damping          = 0;
        rigidConfig.iterMax          = 100;
        rigidConfig.tol              = 1e-12;
        rigidConfig.minActivation    = 0.0;
        
        calcRigidTendonMuscleInfoFcn =...
            @(actState1,pathState2,mclState3)...
            calcMillard2012DampedEquilibriumMuscleInfo(  ...
            actState1,...
            pathState2, ...
            mclState3,...
            muscleArch,...
            normMuscleCurves,...
            rigidConfig);
        
        
        
        benchConfig.numberOfMuscleStates  = 0;
        benchConfig.minimumActivation     = rigidConfig.minActivation;
        benchConfig.name                  = 'RT';
        
        calcInitalRigidMuscleState = [];
        
        %Calculate the normalization factor: normFiberForceAlongTendonIsometric
        if(count == 1)
          activation=1;
          pathState=[0;lpRampMid];
          muscleState = [];
          muscleInfo = calcRigidTendonMuscleInfoFcn(activation,pathState,muscleState);
          cosPennationAngle = muscleInfo.muscleLengthInfo.cosPennationAngle;
          normFiberForce = muscleInfo.muscleDynamicsInfo.normFiberForce;
          normFiberForceAlongTendonIsometric_rt = normFiberForce*cosPennationAngle;
        end
        
        benchRecordRigid = ...
            runMillard2012ComputationalBenchmark(calcRigidTendonMuscleInfoFcn,...
            calcInitalRigidMuscleState ,...
            benchConfig,...
            figRigidTendonBasic,...
            figRigidTendonEnergy,...
            figRigidTendonPower);
        
        save('benchRecordRigid.mat','benchRecordRigid');
        saveas(figRigidTendonBasic,  'figRigidTendonBasic.fig','fig');
        saveas(figRigidTendonEnergy, 'figRigidTendonEnergy.fig','fig');
        saveas(figRigidTendonPower,  'figRigidTendonPower.fig','fig');
        
        
        fiber_force_rt(:,count) = benchRecordRigid.normFiberForceAlongTendon;
        fiber_length_rt(:,count)= benchRecordRigid.normFiberLength;
        path_rt(:,count)        = benchRecordRigid.pathLength;

        fiber_velocity_rt(:,count)            = benchRecordRigid.fiberVelocity;
        fiber_velocityAlongTendon_rt(:,count) = benchRecordRigid.fiberVelocityAlongTendon;

        measurement_time_rt(count,1) = measurementTime;
        sim_time_rt(:,count)         = ([0:(1/(npts-1)):1]').*(totalSimTime);

        ankle_angularVelocity_rt(:,count)  = benchRecordRigid.pathVelocity./ankleAchillesTendonMomentArm;
        pennation_angle_rt(:,count)        =  benchRecordRigid.pennationAngle;          
    end
    
    
    
    
    %%=========================================================================
    %Classic Elastic Tendon Model Benchmark
    %%=========================================================================
    
    if flag_runClassicElasticBench == 1
        disp('Classic elastic-tendon model: ramp-shortening, constant activation');

        
        classicElasticTendonConfig.useFiberDamping  = 0;
        classicElasticTendonConfig.useElasticTendon = 1;
        classicElasticTendonConfig.damping          = 0;
        classicElasticTendonConfig.iterMax          = 100;
        classicElasticTendonConfig.tol              = 1e-6;
        classicElasticTendonConfig.minActivation    = 0.05;
        
        calcClassicElasticTendonMuscleInfoFcn =...
            @(actState1,pathState2,mclState3)...
            calcMillard2012DampedEquilibriumMuscleInfo(  ...
            actState1,...
            pathState2, ...
            mclState3,...
            muscleArch,...
            normMuscleCurves,...
            classicElasticTendonConfig);
        
        calcClassicElasticTendonInitialMuscleStateFcn = ...
            @(actState1,pathState2,calcMuscleInfo3, initConfig4) ...
            calcInitialMuscleState(actState1,...
            pathState2,...
            muscleArch,...
            calcMuscleInfo3,...
            initConfig4);
        
        %Calculate the normalization factor: normFiberForceAlongTendonIsometric
        if(count == 1)
          activation  = 1;
          pathState   = [0;lpRampMid];
                    
          initConfig.iterMax = 100;
          initConfig.tol     = 1e-8;
          initConfig.useStaticFiberSolution = 0;
          initSoln = calcClassicElasticTendonInitialMuscleStateFcn(...
                        activation,...
                        pathState,...                                          
                        calcClassicElasticTendonMuscleInfoFcn,...
                        initConfig);   
                      
          assert(initSoln.converged == 1 || initSoln.isClamped == 1,...
                 'Failed to bring the muscle to a valid initial solution');
               
          muscleState = initSoln.muscleState(:);
          
          muscleInfo = calcClassicElasticTendonMuscleInfoFcn(...
                          activation,pathState,muscleState);
                        
          cosPennationAngle = muscleInfo.muscleLengthInfo.cosPennationAngle;
          normFiberForce    = muscleInfo.muscleDynamicsInfo.normFiberForce;
          normFiberForceAlongTendonIsometric = normFiberForce*cosPennationAngle;
        end

          
        benchConfig.numberOfMuscleStates = 1;
        benchConfig.minimumActivation    = ...
            classicElasticTendonConfig.minActivation;
        benchConfig.name = 'CE';
        
        benchRecordClassicElastic = ...
            runMillard2012ComputationalBenchmark(...
            calcClassicElasticTendonMuscleInfoFcn,...
            calcClassicElasticTendonInitialMuscleStateFcn,...
            benchConfig,...
            figClassicElasticBasic,...
            figClassicElasticEnergy,...
            figClassicElasticPower);
        
        
        save('benchRecordClassicElastic.mat','benchRecordClassicElastic');
        saveas(figClassicElasticBasic,'figClassicElasticBasic.fig','fig');
        saveas(figClassicElasticEnergy,'figClassicElasticEnergy.fig','fig');
        saveas(figClassicElasticPower,'figClassicElasticPower.fig','fig');
        
    end
    %%=========================================================================
    %Damped Equilibrum Elastic Model Benchmark
    %%=========================================================================
    
    if flag_runDampedFiberElasticBench == 1
        disp('Damped-fiber elastic-tendon model: ramp shortening, constant activation');
        
        dampedFiberElasticTendonConfig.useFiberDamping  = 1;
        dampedFiberElasticTendonConfig.useElasticTendon = 1;
        dampedFiberElasticTendonConfig.damping          = 0.1;
        dampedFiberElasticTendonConfig.iterMax          = 100;
        dampedFiberElasticTendonConfig.tol              = 1e-6;
        dampedFiberElasticTendonConfig.minActivation    = 0.0;
        
        calcDampedFiberElasticTendonMuscleInfoFcn =...
            @(actState1,pathState2,mclState3)...
            calcMillard2012DampedEquilibriumMuscleInfo(  ...
            actState1,...
            pathState2, ...
            mclState3,...
            muscleArch,...
            normMuscleCurves,...
            dampedFiberElasticTendonConfig);
        
        calcDampedFiberElasticTendonInitialMuscleStateFcn = ...
            @(actState1,pathState2,calcMuscleInfo3, initConfig4) ...
            calcInitialMuscleState(actState1,...
            pathState2,...
            muscleArch,...
            calcMuscleInfo3,...
            initConfig4);
        
        %Calculate the normalization factor muscleForceNormalization
        if(count == 1)
          activation  = 1;
          pathState   = [0;lpRampMid];
                    
          initConfig.iterMax = 100;
          initConfig.tol     = 1e-8;
          initConfig.useStaticFiberSolution = 0;
          initSoln = calcDampedFiberElasticTendonInitialMuscleStateFcn(...
                        activation,...
                        pathState,...                                          
                        calcDampedFiberElasticTendonMuscleInfoFcn,...
                        initConfig);   
                      
          assert(initSoln.converged == 1 || initSoln.isClamped == 1,...
                 'Failed to bring the muscle to a valid initial solution');
               
          muscleState = initSoln.muscleState(:);
          
          muscleInfo = calcDampedFiberElasticTendonMuscleInfoFcn(...
                          activation,pathState,muscleState);
                        
          cosPennationAngle = muscleInfo.muscleLengthInfo.cosPennationAngle;
          normFiberForce    = muscleInfo.muscleDynamicsInfo.normFiberForce;
          normFiberForceAlongTendonIsometric = normFiberForce*cosPennationAngle;          
        end        
          
        benchConfig.numberOfMuscleStates = 1;
        benchConfig.minimumActivation    = ...
            dampedFiberElasticTendonConfig.minActivation;
        benchConfig.name = 'DFE';
        
        benchRecordDampedFiberElasticTendon = ...
            runMillard2012ComputationalBenchmark(...
            calcDampedFiberElasticTendonMuscleInfoFcn,...
            calcDampedFiberElasticTendonInitialMuscleStateFcn,...
            benchConfig,...
            figDampedFiberElasticBasic,...
            figDampedFiberElasticEnergy,...
            figDampedFiberElasticPower);
        
        save('benchRecordDampedFiberElasticTendon.mat',...
            'benchRecordDampedFiberElasticTendon');
        saveas(figDampedFiberElasticBasic,...
            'figDampedFiberElasticBasic.fig','fig');
        saveas(figDampedFiberElasticEnergy,...
            'figDampedFiberElasticEnergy.fig','fig');
        saveas(figDampedFiberElasticPower,...
            'figDampedFiberElasticPower.fig','fig');
        
        fiber_force(:,count)=...
          benchRecordDampedFiberElasticTendon.normFiberForceAlongTendon;
        fiber_length(:,count)= ...
          benchRecordDampedFiberElasticTendon.normFiberLength;
        path(:,count) = ...
          benchRecordDampedFiberElasticTendon.pathLength;

        fiber_velocity(:,count) =...
          benchRecordDampedFiberElasticTendon.fiberVelocity;
        fiber_velocityAlongTendon(:,count) = ...
          benchRecordDampedFiberElasticTendon.fiberVelocityAlongTendon;

        measurement_time(count,1)=measurementTime;
        sim_time(:,count) = ([0:(1/(npts-1)):1]').*(totalSimTime);

        ankle_angularVelocity(:,count)= ...
          benchRecordDampedFiberElasticTendon.pathVelocity./ankleAchillesTendonMomentArm;

        pennation_angle(:,count) = ...
          benchRecordDampedFiberElasticTendon.pennationAngle;          
          
    end
    

    
end %Velocity loop

%% figures
if(flag_updateExistingPlots==0)
  fig_force = figure;
  fig_lengthA= figure;
  fig_lengthB=figure;
end

for count=1:1:countMax
  figure(fig_force);
  plot( path(:,count),...
        fiber_force(:,count),...
        'Color',colorRecord(count,:));
  hold on;
  figure(fig_lengthA);
  plot( fiber_length(:,count),...
        fiber_force(:,count),...
        'Color',colorRecord(count,:));
  hold on;
  figure(fig_lengthB);
  plot( fiber_length(:,count),...
        path(:,count),'Color',colorRecord(count,:));
  hold on;
end

figure(fig_force);
xlabel('pathLength')
ylabel('normFiberForceAlongTendon')

figure(fig_lengthA);
xlabel('normFiberLength')
ylabel('normFiberForceAlongTendon')

figure(fig_lengthB);
xlabel('normFiberLength')
ylabel('pathLength')


%%
%
% Denis - the plot you are looking for is below
%
%%
if(countMax > 2)
  measurementLength                   = mean(path(:,1));
  measuredForce                       = zeros(size(path,2),1);
  measuredFiberVelocity               = zeros(size(path,2),1);
  measuredFiberVelocityAlongTendon    = zeros(size(path,2),1);
  measuredAnkleAngularVelocity        = zeros(size(path,2),1);
  pennationAngle                      = zeros(size(path,2),1);

  measurementLength_rt                = mean(path_rt(:,1));
  measuredForce_rt                    = zeros(size(path_rt,2),1);
  measuredFiberVelocity_rt            = zeros(size(path_rt,2),1);
  measuredFiberVelocityAlongTendon_rt = zeros(size(path_rt,2),1);
  measuredAnkleAngularVelocity_rt     = zeros(size(path_rt,2),1);
  pennationAngle_rt                   = zeros(size(path_rt,2),1);
  
  
  for n=1:1:size(path,2)

      measuredForce(n,1) = interp1(sim_time(:,n), ...
                                   fiber_force(:,n), ...
                                   measurement_time(n));

      measuredFiberVelocity(n,1) = interp1(sim_time(:,n), ...
                                           fiber_velocity(:,n), ...
                                           measurement_time(n));
                                         
      measuredFiberVelocityAlongTendon(n,1) = interp1(sim_time(:,n), ...
                                                fiber_velocityAlongTendon(:,n), ...
                                                measurement_time(n));                                         
 
      measuredAnkleAngularVelocity(n,1) = interp1(sim_time(:,n), ...
                                                ankle_angularVelocity(:,n), ...
                                                measurement_time(n));   
                                              
      pennationAngle(n,1) = interp1(sim_time(:,n), ...
                                    pennation_angle(:,n), ...
                                    measurement_time(n));   
                                              
  end
  
  for n=1:1:size(path_rt,2)

      measuredForce_rt(n,1) = interp1(sim_time_rt(:,n), ...
                                      fiber_force_rt(:,n), ...
                                      measurement_time_rt(n));

      measuredFiberVelocity_rt(n,1) = interp1(sim_time_rt(:,n), ...
                                              fiber_velocity_rt(:,n), ...
                                              measurement_time_rt(n));
                                         
      measuredFiberVelocityAlongTendon_rt(n,1) = ...
        interp1(sim_time_rt(:,n), ...
                fiber_velocityAlongTendon_rt(:,n), ...
                measurement_time_rt(n));                                         
 
      measuredAnkleAngularVelocity_rt(n,1) = ...
        interp1(sim_time_rt(:,n), ...
                ankle_angularVelocity_rt(:,n), ...
                measurement_time_rt(n));   
                                              
      pennationAngle_rt(n,1) = ...
        interp1(sim_time_rt(:,n), ...
                pennation_angle_rt(:,n), ...
                measurement_time_rt(n));   
                                              
  end  
  
  if(flag_updateExistingPlots==0)
    fig_forceAngularVelocity = figure;
    fig_forceFiberVelocity=figure;
    fig_ankleOmegaVsFiberVelocity=figure;   
    fig_ankleOmegaVsPennationAngle=figure;
  end
  figure(fig_forceAngularVelocity)
  legendEntry = sprintf('Elastic tendon (%1.3f)',tendonStrainAtOneNormForceOverride);
  legendEntry_rt = 'Rigid tendon';
  
  plot(ankleAngularVelocity,...
       measuredForce./normFiberForceAlongTendonIsometric,'o-',...
       'DisplayName',legendEntry);
  hold on;
  
  plot(ankleAngularVelocity,...
       measuredForce_rt./normFiberForceAlongTendonIsometric_rt,'o-',...
       'DisplayName',legendEntry_rt);
  hold on;

  
  xlabel('Ankle Angular Velocity (deg/s)');
  ylabel('Norm. Fiber Force (N/fiso*)');
  title('Musculotendon force at a neutral ankle angle (*)');
  ylim([0,1]);
  xlim([0,max(ankleAngularVelocity)]);

  legend show

  figure(fig_forceFiberVelocity)
  plot(measuredFiberVelocity.*-1000,...
       measuredForce./normFiberForceAlongTendonIsometric,'o-',...
       'DisplayName',legendEntry);
  hold on;
  xlabel('Fiber Shortening Speed (mm/s)');
  ylabel('Norm. Fiber Force (N/fiso*)');
  title('Musculotendon force at a neutral ankle angle (*)');
  ylim([0,1]);
  xlim([min(measuredFiberVelocity.*-1000),max(measuredFiberVelocity.*-1000)]);

  legend show
  
  figure(fig_ankleOmegaVsFiberVelocity)
  plot(measuredFiberVelocity.*-1000,...
       measuredAnkleAngularVelocity.*(-180/pi),'o-',...
       'DisplayName',legendEntry);
  hold on;
  plot(measuredFiberVelocityAlongTendon.*-1000,...
       measuredAnkleAngularVelocity.*(-180/pi),'x-',...
       'DisplayName',[legendEntry,'AT']);
     
  xlabel('Fiber Shortening Speed (mm/s)');
  ylabel('Ankle Angular Velocity (deg/s)');
  ylim([min(measuredAnkleAngularVelocity.*(-180/pi)),max(measuredAnkleAngularVelocity.*(-180/pi))]);
  %xlim([min(measuredFiberVelocity.*-1000),max(measuredFiberVelocity.*-1000)]);

  legend show
  
  figure(fig_ankleOmegaVsPennationAngle)
  plot(pennationAngle.*(180/pi),...
       measuredAnkleAngularVelocity.*(-180/pi),'o-',...
       'DisplayName',legendEntry);
  hold on;
     
  xlabel('Pennation Angle (deg)');
  ylabel('Ankle Angular Velocity (deg/s)');
  ylim([min(measuredAnkleAngularVelocity.*(-180/pi)),max(measuredAnkleAngularVelocity.*(-180/pi))]);

end



String_Neu = ['AAA' num2str(tendonStrainAtOneNormForceOverride*1000)];
Werte.(String_Neu) = measuredForce./normFiberForceAlongTendonIsometric;

%%
figure(1)
  hold on
  plot(measuredFiberVelocity.*-1000,...
       measuredForce./normFiberForceAlongTendonIsometric,'o-',...
       'DisplayName',legendEntry);
  hold on;
  plot(measuredFiberVelocity_rt.*-1000,...
       measuredForce_rt./normFiberForceAlongTendonIsometric_rt,'o-',...
       'DisplayName',legendEntry_rt);
  
figure(2)
  hold on
  plot(ankleAngularVelocity,...
       measuredForce./normFiberForceAlongTendonIsometric,'o-',...
       'DisplayName',legendEntry);
  hold on;
  plot(ankleAngularVelocity,...
       measuredForce_rt./normFiberForceAlongTendonIsometric_rt,'o-',...
       'DisplayName',legendEntry_rt);     
%%  
% figure(11)
% hold on
% plot(Hauraix_FSZ_velo,Hauraix_angular_velo,'DisplayName','Hauraix')