flag_outerLoopMode = 1;

disp('Update preload protocol to begin rotating at 95% activation');

if(flag_outerLoopMode==0)
  clc;
  close all;
  clear all;
  
  flag_simForceVelocityExpWithPreload =0;
  preloadFraction = 0;
  % 0: Shortening begins simultaneously with the development of force
  % 1: Shortening begins after 100% activation is reached.

  fractionOfFastTwitchFibers         = 0.5;
  % 0: Gives you a force-velocity curve constent with a slow-twitch fiber
  % 1: Gives you a force-velocity curve constent with a fast-twitch fiber
  
  ankleAchillesTendonMomentArm = 0.054;%/1.18; %moment arm in m  
  
  tendonStrainAtOneNormForceOverride = 0.1;
 
  flag_measurementSetting = 0;%
  %0: measueAtAnkleAngle
  %1: measureAtNormFiberLength;
  flag_useFlatActiveForceLengthCurve  = 0;
  flag_useConstantTendonStiffness     = 0;
  flag_useLinearForceVelocityCurve    = 0;
  
  flag_useTendonDamping = 0;
  normalizedTendonDamping         = 0.05*0.9;
  normalizedTendonDampingConstant = 0.05*0.1;
  
  maximumNormalizedFiberVelocity = 10; %in units of norm fiber lengths/second
  scaleLceOpt = 1;
  ankleAngleMaxPlantarFlexion   = -17;%-17 is from Holzer, -23; (Rugg et al.)   %in degrees
end

scaleLceOptStr = ['lceOptScale_',num2str(round(scaleLceOpt*100)),'p_'];


%assert( flag_useTendonDamping == 0,'Not yet implemented!');

% Running this script will perform the desired simulations and write
% the results to a series of MAT files that are stored in the data folder.
% Running the post-processing script "main_PlotMaxActivationRampShortening.m"
% will generate the desired plots.





fractionOfFastTwitchFibersStr = ...
  ['FastTwitch_',num2str(round(fractionOfFastTwitchFibers*100,0))];

preloadStr = [num2str(round(preloadFraction*100,0)),'p'];

ankleAchillesTendonMomentArmStr = sprintf('_%1.1fcm_',...
  round(ankleAchillesTendonMomentArm*100,1));
ankleAchillesTendonMomentArmStr(1,end-4)='p';

measureAtAnkleAngle  = 0;
  %from Rugg et al.  
  %ankleAngleMaxPlantarFlexion   = -17;%-17 is from Holzer, -23; (Rugg et al.)   %in degrees

  
measureAtNormFiberLength = 1;
  normFiberLengthAtMeasurement = 0.85; %A length common to all tests
  
assert( xor(measureAtAnkleAngle,measureAtNormFiberLength) == 1);

%flag_measurementSetting = measureAtNormFiberLength;
measurementSettingStr = '';

if(flag_measurementSetting == measureAtNormFiberLength)
  measurementSettingStr = 'fixedFiberLength';
end
if(flag_measurementSetting == measureAtAnkleAngle)
  measurementSettingStr = 'fixedAnkleAngle';  
end

%Note: To update the plots:
% 1. Set fractionOfFastTwitchFibers to the desired value.
% 2. Run once with          flag_simForceVelocityExpWithPreload =0
% 3. Run a second time with flag_simForceVelocityExpWithPreload =1
% 4. Then run the plot script.
%
% Warning! This overwrites the MAT files in the output folder: normally
% this is ok. If you want to put the fast-vs-slow twitch on one plot you
% will need to update the file naming convention so that 
%
%  simFv_gaslat_preload_0_dampedFiberElasticTendon.mat
%  simFv_gaslat_preload_1_dampedFiberElasticTendon.mat
%  muscleArch.mat
%  normMuscleCurves.mat
%
%


flag_generateDiagnosticPlots        = 0; %Basic/Energy/Power plots
flag_updateExistingPlots            = 0; %Only applies to the diagnostic plots



if(tendonStrainAtOneNormForceOverride <0)
  tendonStrainAtOneNormForceOverride = 0.049; %The default
end
tendonStrainAtOneNormForceOverrideStr = ...
  sprintf('_Tdn_%1.1f',round(tendonStrainAtOneNormForceOverride*100,1));
tendonStrainAtOneNormForceOverrideStr(1,end-1)='p';

if(flag_runRigidBench==1)
  tendonStrainAtOneNormForceOverrideStr = ...
    sprintf('_Tdn_%1.1f',round(0,1));
  tendonStrainAtOneNormForceOverrideStr(1,end-1)='p';
  
end

flag_useArnold2010SoleusArchitecture = 1;
muscle                             = 'gasmed';
initialHoldTime                    = 0.1;
shiftFiberForceLengthCurve         = 0.;
%%
% Denis: You can change the tendon strain when one norm force is applied by 
%        changing the above parameter. Here are some values to try
%
% 0.049     : The default value used in OpenSim 
% 0.10-0.15 : A range of values that is typical for the Achilles tendon 
%             for EMG-driven models of running          
%%

outputDataFolder = '../data/';
outputFileName   = sprintf('simFv_%s_preload_%s',...
                           muscle,preloadStr);

% flag_runRigidBench               = 0;
% flag_runClassicElasticBench      = 0;
% flag_runDampedFiberElasticBench  = 1;
% 
% if( abs(tendonStrainAtOneNormForceOverride) < 1e-3)
%   flag_runRigidBench               = 1;
%   flag_runClassicElasticBench      = 0;
%   flag_runDampedFiberElasticBench  = 0;  
% else
%   flag_runRigidBench               = 0;
%   flag_runClassicElasticBench      = 0;
%   flag_runDampedFiberElasticBench  = 1;  
% end


if(flag_generateDiagnosticPlots==1)
  if(flag_runRigidBench==1 && flag_updateExistingPlots==0 )
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
else
    figRigidTendonBasic         = [];
    figRigidTendonEnergy        = [];
    figRigidTendonPower         = [];
    figClassicElasticBasic      = [];
    figClassicElasticEnergy     = [];
    figClassicElasticPower      = [];  
    figDampedFiberElasticBasic  = [];
    figDampedFiberElasticEnergy = [];
    figDampedFiberElasticPower  = [];    
end

count=0;
countMax=21;

labelsHauraix2015 = {'w30','w90','w150','w210','w270','w330'};
Hauraix2015Fig3A   = loadDigitizedData('../dataHauraix2015/HauraixFig3A.csv',...
                                    'x','y',labelsHauraix2015,'Fig3A');
for z=1:1:length(Hauraix2015Fig3A)
  Hauraix2015Fig3A(z).x = 90-(Hauraix2015Fig3A(z).x-90);
end

angularVelocitiesHauraix2015 = [1;30;90;150;210;270;330];


if(flag_rampType==1)
  countMax = length(angularVelocitiesHauraix2015);
end

colorStart = [0,0,1];
colorEnd   = [1,0,0];
colorRecord = zeros(countMax,3);

ankleAngularVelocity = zeros(countMax,1);
fiberVelocity        = zeros(countMax,1);
normFiberForceAlongTendonIsometric = 0; %Isometric force generated by the muscle
                                        %mid-ramp
normFiberForceAlongTendonIsometric_rt = 0;

npts = 1000;



protocolName ='';
if(flag_simForceVelocityExpWithPreload==1)
  protocolName = 'Preload';
else
  protocolName = 'NoPreload';
end


detailedResultsStruct = struct(...
                         'ankleAngle',zeros(npts,countMax),...
                         'ankleAngularVelocity',zeros(npts,countMax),...
                         'simulationTime',zeros(npts,countMax),...
                         'measurementTimeAnkleAngle',zeros(1,countMax),...
                         'normFiberForceAlongTendonIsometric',zeros(1,countMax));

emptyBenchRecord = getEmptyBenchRecord(npts,countMax);

rigidTendonSimulationRecord = struct( 'muscleArchitecture',[],... 
                                      'simulationConfiguration',[],...
                                      'standardResults',[],...
                                      'detailedResults',[],... 
                                      'muscle',muscle,...
                                      'protocol', protocolName);

rigidTendonSimulationRecord.standardResults = emptyBenchRecord;
rigidTendonSimulationRecord.detailedResults = detailedResultsStruct;


classicElasticTendonSimulationRecord = struct('muscleArchitecture',[],... 
                                              'simulationConfiguration',[],...
                                              'standardResults',[],...
                                              'detailedResults',[],...
                                              'muscle',muscle,...
                                              'protocol', protocolName);

classicElasticTendonSimulationRecord.standardResults = emptyBenchRecord;
classicElasticTendonSimulationRecord.detailedResults = detailedResultsStruct;


dampedFiberElasticTendonSimulationRecord = struct('muscleArchitecture',[],... 
                                                  'simulationConfiguration',[],...
                                                  'standardResults',[],...
                                                  'detailedResults',[],...
                                                  'muscle',muscle,...
                                                  'protocol', protocolName);

dampedFiberElasticTendonSimulationRecord.standardResults = emptyBenchRecord;
dampedFiberElasticTendonSimulationRecord.detailedResults = detailedResultsStruct;



for count=1:1:countMax
    fprintf('%i of %i\n',count,countMax);
        
    n = (count-1)*10+1;
    if(flag_rampType==1)
      n = angularVelocitiesHauraix2015(count,1);
    end
    %%

    muscleAbbrArnold2010                 = muscle;
    
    flag_plotNormMuscleCurves            = 0;
    flag_updateNormMuscleCurves          = 1;
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
    
    rampStartAngle               = 0;   
    rampEndAngle                 = 0;   
    rampStr                      = '';
    rampMeasurementAngle         = 0;
    
    switch(flag_rampType)
      
      case 0
        rampStartAngle               = -15;   %degrees, -ve: dorsiflexion
        rampEndAngle                 =  15;    %degrees
        rampStr = 'ramp_Holzer_';
        rampMeasurementAngle         = 0;
      case 1
        rampStartAngle               = -20;   %degrees, -ve: dorsiflexion
        rampEndAngle                 =  35;    %degrees
        rampMeasurementAngle         = 0;
        
        rampStr = 'ramp_Hauraix2015_';
        
      otherwise
        asert(0,'flag_rampType must be 0 (Holzer) or 1 (Hauraix 2015)');
    end
    
    rampAngularVelocity          =  n;  %degrees per second.
    ankleAngularVelocity(count)  = n;
    %Parameters from the literature
      
    %From Anderson et al.
    
    %%
    %Parameters that apply to all muscles
    %%
    %maximumNormalizedFiberVelocity = 10; %in units of norm fiber lengths/second
    %vmaxStr = ['vmax_',num2str(round(maximumNormalizedFiberVelocity,0)),'_'];
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
    
    if(count==1)
      flag_updateNormMuscleCurves=1;
    else
      flag_updateNormMuscleCurves=0;
    end
    
    
    normMuscleCurves = ...
        createDefaultNormalizedMuscleCurves(muscleAbbr,...
        tendonStrainAtOneNormForceOverride,...
        fractionOfFastTwitchFibers,...
        fractionOfFastTwitchFibersStr,...
        shiftFiberForceLengthCurve,...
        outputDataFolder,...        
        flag_useFlatActiveForceLengthCurve,...
        flag_useLinearForceVelocityCurve,...
        flag_useConstantTendonStiffness,...
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
        
    lceOpt = lceOpt*scaleLceOpt;
    
    vmaxStr = ['vmax_',num2str(round(maximumNormalizedFiberVelocity,0)),'_']; 
    if(flag_useHauraixVmax==1)
      maximumNormalizedFiberVelocity = 0.308/lceOpt;
      vmaxStr = ['vmax_Hauraix_']; 
    end

         
    
    %%
    % If you want to alter the achitectural properties of the
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
    benchConfig.relTol = 1e-9;
    benchConfig.absTol = 1e-9;
    benchConfig.initialActivationState = 0;%[0:0.1:1]';
    benchConfig.color = trialColor;
    
    
    %%
    % Musculotendon Ramp-Shortening Function: Two different protocols
    %
    % With pre-load
    %  1. Musuclotendon path is at a length consistent with 15 degrees of 
    %     dorsiflexion
    %  2. Muscle is maximally activated
    %  3. The path is shortened to a length consistent with 15 degrees of 
    %     plantar flexion in certain amount of time.
    %  4. The tension of the musculotendon is sampled at the half way through
    %     the ramp
    %
    % Without pre-load
    %
    %  1. Musuclotendon path is at a length consistent with 15 degrees of 
    %     dorsiflexion
    %  2. Simultaneously the muscle is activated and the path is shortened
    %     to a length consistent with 15 degrees of plantar flexion.
    %  3. The tension of the musculotendon is sampled at the half way through
    %     the ramp
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
    % II. Change in Musculotendon Path Length
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
    % dl/dt = d omega/dt * r
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
    %  joint torque as a function of joint angle and angular velocity:close
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
    
    save([outputDataFolder,...
      'muscleArch_',fractionOfFastTwitchFibersStr,'.mat'],'muscleArch');
    
    %%
    % I. Starting Musculotendon Path Length
    %%
    
    %Solve for the length of the tendon when the fiber is at its optimal
    %length, optimal pennation angle, and is fully activated
    ftNOpt     = 1*cos(alphaOpt);
    ltNAtFiso  = normMuscleCurves.tendonForceLengthCurve.xEnd(2);
    ftNerr     = 1;
    tol        = 1e-9;
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

    lpDelta90 = (ankleAngleMaxPlantarFlexion)*(pi/180)*ankleAchillesTendonMomentArm;
    
    
    lpDeltaStart = -(rampStartAngle)*(pi/180)*ankleAchillesTendonMomentArm;
    
    lpRampStart = lpOpt + lpDelta90 + lpDeltaStart;
    
    %%
    %
    %%
    
    
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
    
    
    
    %totalSimTime = 1;
    %minHoldTime  = 0.25;
    %if( (totalSimTime-(initialHoldTime+rampTime)) < minHoldTime)
    %  totalSimTime = (initialHoldTime+rampTime)+minHoldTime;
    %end
    
    totalSimTime = initialHoldTime + rampTime;
    
    %measurementTimeAnkleAngle = initialHoldTime + 0.5*rampTime;
    measurementAngleFraction = (rampMeasurementAngle-rampStartAngle) ...
                              /(rampEndAngle-rampStartAngle);
    measurementTimeAnkleAngle = ...
      initialHoldTime + measurementAngleFraction*rampTime;
    %%
    % Construct the path function
    %%
    
    switch flag_rampType
      case 0    
        pathFcn = @(t)calcRampFunctionState(t,initialHoldTime,...
            lpRampStart,lpRampEnd,rampTime);
      case 1
        %construct a function that has the same velocity profile as 
        %Hauraix 2015 Fig 3A.
        if(count == 1)
          pathFcn = @(t)calcRampFunctionState(t,initialHoldTime,...
              lpRampStart,lpRampEnd,rampTime);          
        else
          
          odefun = @(argT,argX)calcHauraixCrankDerivative(argT,argX,...
                               Hauraix2015Fig3A(count-1));
          tspan  = [0,5];
          angle0 = 70; 
          angle1 = 122;
          eventfun = @(argT,argX)calcHauraixEvent(argT,argX,angle1);
          options = odeset('events',eventfun,'RelTol',1e-6,'AbsTol',1e-6);
          [t,y,te,ye,ie] = ode45(odefun,tspan,angle0,options);
          
          dy = zeros(size(y));
          for z=1:1:length(y)
            dy(z,1) = odefun(t(z,1),y(z,1));
          end
  
          timeSeries = t;
          angleSeries = y.*(pi/180);
          angularVelocitySeries = dy.*(pi/180);
          
          lengthSeries = lpRampStart ...
                       - (angleSeries-angleSeries(1,1)...
                         ).*ankleAchillesTendonMomentArm;
          velocitySeries =  -(angularVelocitySeries ...
                             ).*ankleAchillesTendonMomentArm;
            
          
          flag_debugHauraixFcn=0;
          if(flag_debugHauraixFcn==1)
            figHauraixFcn =figure;
            subplot(2,2,1);
              plot(t,y,'r');
              xlabel('Time');
              ylabel('Angle');
            subplot(2,2,2);
              plot(y,dy,'r','LineWidth',2);
              hold on;
              plot(Hauraix2015Fig3A(count-1).x,...
                   Hauraix2015Fig3A(count-1).y,'k');
              hold on;
              xlabel('Angle');
              ylabel('Angular Velocity');
            subplot(2,2,3);
              plot([0;rampTime],[lpRampStart;lpRampEnd],'b','LineWidth',2);
              hold on;
              plot(timeSeries,lengthSeries,'r');
              hold on;
              xlabel('Time (s)');
              ylabel('Length (m)');
            subplot(2,2,4);
              plot([0;rampTime],[1;1].*((lpRampEnd-lpRampStart)/rampTime),'b','LineWidth',2);
              hold on;
              plot(timeSeries,velocitySeries,'r');
              hold on;
              xlabel('Time (s)');
              ylabel('Lengthening Rate (m/s)');              
            here=1;
          end
          
          
          pathFcn = @(argT)calcPathFunctionState(argT,initialHoldTime,...
                            timeSeries, lengthSeries,velocitySeries);
          totalSimTime = max(timeSeries)+initialHoldTime;
          rampAngularVelocity = max(dy);
          %pathFcn = @(t)calcRampFunctionState(t,initialHoldTime,...
          %lpRampStart,lpRampEnd,rampTime);   
        
        end
        
        
        
      otherwise assert(0);
    end
    
    benchConfig.pathFcn = pathFcn;
    benchConfig.tspan   = [0, (totalSimTime)];
    
    fprintf('  %1.3f deg/s, preload %i\n',rampAngularVelocity,...
              flag_simForceVelocityExpWithPreload);
    
    %
    %
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
    if(flag_simForceVelocityExpWithPreload == 0)
      ton  = initialHoldTime;
      toff = totalSimTime;      
    end
    
    excitationFcn = @(argTime)calcStepFunction(argTime, ton, toff);    
    activationFcn = @(argEx, argAct)calcFirstOrderActivationDerivative(...
                                    argEx, argAct, tact,tdeact,tmin);
    
    benchConfig.excitationFcn = excitationFcn;
    benchConfig.activationFcn = activationFcn;
    
    %%=========================================================================
    %Rigid Tendon Model Benchmark
    %%=========================================================================
    
    if(flag_runRigidBench == 1)
        disp('  Rigid-tendon model');
        %figRigidTendonBasic  = figure;
        %figRigidTendonEnergy = figure;
        %figRigidTendonPower  = figure;
        
        rigidConfig.useFiberDamping  = 0;
        rigidConfig.useElasticTendon = 0;
        rigidConfig.damping          = 0;
        rigidConfig.iterMax          = 100;
        rigidConfig.tol              = 1e-12;
        rigidConfig.minActivation    = 0.0;
        rigidConfig.useTendonDamping = 0;
        rigidConfig.normalizedTendonDamping =normalizedTendonDamping;
        rigidConfig.normalizedTendonDampingConstant =normalizedTendonDampingConstant;
                
        tendonDampingStr = ['_TendonDamping_',...
        [num2str(round((rigidConfig.useTendonDamping ...
                       *rigidConfig.normalizedTendonDamping)*100,0)),'p']];        
        
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
        

        
        benchConfig.initialActivationState = 1;
        if(flag_simForceVelocityExpWithPreload==0)
          benchConfig.initialActivationState = 0;
        else
          benchConfig.initialActivationState = preloadFraction;
        end
        benchConfig.nSim = 1;        
        
        rigidTendonSimulationRecord.muscleArchitecture      = muscleArch;        
        rigidTendonSimulationRecord.simulationConfiguration = benchConfig;

        benchRecordRigid = ...
            runMillard2012ComputationalBenchmark(calcRigidTendonMuscleInfoFcn,...
            calcInitalRigidMuscleState ,...
            benchConfig,...
            figRigidTendonBasic,...
            figRigidTendonEnergy,...
            figRigidTendonPower);

        rigidTendonSimulationRecord.standardResults = ...
          insertResultsIntoSet( benchRecordRigid,...
                                rigidTendonSimulationRecord.standardResults,...
                                count);

        rigidTendonSimulationRecord.detailedResults.ankleAngularVelocity(:,count)= ...
          benchRecordRigid.pathVelocity./ankleAchillesTendonMomentArm;
        
        rigidTendonSimulationRecord.detailedResults.ankleAngle(:,count)=...
          + (90 + rampStartAngle)*(pi/180)  ...
          - (rigidTendonSimulationRecord.standardResults.pathLength(:,count) - lpRampStart...
             )*(1/ankleAchillesTendonMomentArm);  
        
        
        rigidTendonSimulationRecord.detailedResults.simulationTime(:,count) = ...
          ([0:(1/(npts-1)):1]').*(totalSimTime);

        measurementTime = 0;
        if(flag_measurementSetting == measureAtAnkleAngle)
          measurementTime = measurementTimeAnkleAngle;          
        end   
        
        if(flag_measurementSetting == measureAtNormFiberLength)
          measurementTime = ...
            getTimeAtNormFiberLength(rigidTendonSimulationRecord, ...
                                     normFiberLengthAtMeasurement,count);          
        end
        rigidTendonSimulationRecord.detailedResults.measurementTime(1,count) = ...
          measurementTime;

        %Get the isometric force of the muscle at the sample point
        if(count ==1)
          normFiberForceAlongTendonIsometric = ...
            interp1(...
              rigidTendonSimulationRecord.standardResults.time(:,count),...
              rigidTendonSimulationRecord.standardResults.normFiberForceAlongTendon(:,count),...
              measurementTime);
        end        
        
        rigidTendonSimulationRecord.detailedResults.normFiberForceAlongTendonIsometric(1,count) = ...
          normFiberForceAlongTendonIsometric;
        
        
        if(count==countMax)
          disp(['  saved to: ',outputDataFolder, outputFileName,...
            '_rigidTendon_',...
            scaleLceOptStr,...            
            vmaxStr,...
            rampStr,...
            fractionOfFastTwitchFibersStr,...
            tendonStrainAtOneNormForceOverrideStr,...
            tendonDampingStr,...
            ankleAchillesTendonMomentArmStr,...
            measurementSettingStr,'.mat']);
          
          save( [ outputDataFolder, outputFileName,...
                  '_rigidTendon_',...
                  scaleLceOptStr,...  
                  vmaxStr,...
                  rampStr,...
                  fractionOfFastTwitchFibersStr,...
                  tendonStrainAtOneNormForceOverrideStr,...
                  tendonDampingStr,...
                  ankleAchillesTendonMomentArmStr,...
                  measurementSettingStr,'.mat'],...
                'rigidTendonSimulationRecord');
        end
    end
    
    
    
    
    %%=========================================================================
    %Classic Elastic Tendon Model Benchmark
    %%=========================================================================
    
    if flag_runClassicElasticBench == 1
        disp('  Classic elastic-tendon model');

        
        classicElasticTendonConfig.useFiberDamping  = 0;
        classicElasticTendonConfig.useElasticTendon = 1;
        classicElasticTendonConfig.damping          = 0;
        classicElasticTendonConfig.iterMax          = 100;
        classicElasticTendonConfig.tol              = 1e-9;
        classicElasticTendonConfig.minActivation    = 0.05;
        classicElasticTendonConfig.useTendonDamping = 0; %Cannot be used here.
        classicElasticTendonConfig.normalizedTendonDamping =normalizedTendonDamping;
        classicElasticTendonConfig.normalizedTendonDamping =normalizedTendonDampingConstant;
        
        tendonDampingStr = ['_TendonDamping_',...
        [num2str(round((classicElasticTendonConfig.useTendonDamping ...
                       *classicElasticTendonConfig.normalizedTendonDamping)*100,0)),'p']];        
        

        classicElasticTendonSimulationRecord.muscleArchitecture = muscleArch;

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
          
        benchConfig.numberOfMuscleStates = 1;
        benchConfig.minimumActivation    = ...
            classicElasticTendonConfig.minActivation;
        benchConfig.name = 'CE';
        
        
        benchConfig.initialActivationState = 1;
        if(flag_simForceVelocityExpWithPreload==0)
          benchConfig.initialActivationState = 0;
        else
          benchConfig.initialActivationState = preloadFraction;
        end        
        
        benchConfig.nSim = 1;        
        
        classicElasticTendonSimulationRecord.muscleArchitecture      = muscleArch;        
        classicElasticTendonSimulationRecord.simulationConfiguration = benchConfig;


        benchRecordClassicElastic = ...
            runMillard2012ComputationalBenchmark(...
            calcClassicElasticTendonMuscleInfoFcn,...
            calcClassicElasticTendonInitialMuscleStateFcn,...
            benchConfig,...
            figClassicElasticBasic,...
            figClassicElasticEnergy,...
            figClassicElasticPower);
        
        
        classicElasticTendonSimulationRecord.standardResults = ...
          insertResultsIntoSet( benchRecordClassicElastic,...
                                classicElasticTendonSimulationRecord.standardResults,...
                                count);

        classicElasticTendonSimulationRecord.detailedResults.ankleAngularVelocity(:,count)= ...
          benchRecordClassicElastic.pathVelocity./ankleAchillesTendonMomentArm;

        classicElasticTendonSimulationRecord.detailedResults.ankleAngle(:,count)=...
          + (90 + rampStartAngle)*(pi/180)  ...
          - (classicElasticTendonSimulationRecord.standardResults.pathLength(:,count) - lpRampStart...
             )*(1/ankleAchillesTendonMomentArm);        
        
        classicElasticTendonSimulationRecord.detailedResults.simulationTime(:,count) = ...
          ([0:(1/(npts-1)):1]').*(totalSimTime);

        measurementTime = 0;
        if(flag_measurementSetting == measureAtAnkleAngle)
          measurementTime = measurementTimeAnkleAngle;          
        end           
        if(flag_measurementSetting == measureAtNormFiberLength)
          measurementTime = ...
            getTimeAtNormFiberLength(classicElasticTendonSimulationRecord, ...
                                     normFiberLengthAtMeasurement,count);          
        end
        classicElasticTendonSimulationRecord.detailedResults.measurementTime(1,count) = ...
          measurementTime;        

        %Get the isometric force of the muscle at the sample point
        if(count ==1)
          normFiberForceAlongTendonIsometric = ...
            interp1(...
              classicElasticTendonSimulationRecord.standardResults.time(:,count),...
              classicElasticTendonSimulationRecord.standardResults.normFiberForceAlongTendon(:,count),...
              measurementTime);
        end
        
        classicElasticTendonSimulationRecord.detailedResults.normFiberForceAlongTendonIsometric(1,count) = ...
          normFiberForceAlongTendonIsometric;        

        if(count==countMax)
          disp(['  saved to: ',outputDataFolder, outputFileName,...
            '_classicElasticTendon_',...
            scaleLceOptStr,...  
            vmaxStr,...
            rampStr,...
            fractionOfFastTwitchFibersStr,...
            tendonStrainAtOneNormForceOverrideStr,...
            tendonDampingStr,...
            ankleAchillesTendonMomentArmStr,...
            measurementSettingStr,'.mat']);          
          save( [outputDataFolder, outputFileName,...
            '_classicElasticTendon_',...
            scaleLceOptStr,...  
            vmaxStr,...
            rampStr,...
            fractionOfFastTwitchFibersStr,...
            tendonStrainAtOneNormForceOverrideStr,...
            tendonDampingStr,...
            ankleAchillesTendonMomentArmStr,...
            measurementSettingStr,'.mat'],...
                'classicElasticTendonSimulationRecord');
        end
        
    end
    %%=========================================================================
    %Damped Equilibrum Elastic Model Benchmark
    %%=========================================================================
    
    if flag_runDampedFiberElasticBench == 1
        disp('  Damped-fiber elastic-tendon model');
        
        dampedFiberElasticTendonConfig.useFiberDamping  = 1;
        dampedFiberElasticTendonConfig.useElasticTendon = 1;
        dampedFiberElasticTendonConfig.damping          = 0.1;
        dampedFiberElasticTendonConfig.iterMax          = 100;
        dampedFiberElasticTendonConfig.tol              = 1e-9;
        dampedFiberElasticTendonConfig.minActivation    = 0.0;
        
        dampedFiberElasticTendonConfig.useTendonDamping = flag_useTendonDamping;
        dampedFiberElasticTendonConfig.normalizedTendonDamping =normalizedTendonDamping;
        dampedFiberElasticTendonConfig.normalizedTendonDampingConstant =normalizedTendonDampingConstant;
        
        tendonDampingStr = ['_TendonDamping_',...
        [num2str(round((dampedFiberElasticTendonConfig.useTendonDamping ...
                       *dampedFiberElasticTendonConfig.normalizedTendonDamping)*100,0)),'p']];        
                
        
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
        
             
          
        benchConfig.numberOfMuscleStates = 1;
        benchConfig.minimumActivation    = ...
            dampedFiberElasticTendonConfig.minActivation;
        benchConfig.name = 'DFE';
                
        benchConfig.initialActivationState = 1;
        if(flag_simForceVelocityExpWithPreload==0)
          benchConfig.initialActivationState = 0;
        else
          benchConfig.initialActivationState = preloadFraction;
        end
        
        benchConfig.nSim = 1;
        
        dampedFiberElasticTendonSimulationRecord.muscleArchitecture      = muscleArch;        
        dampedFiberElasticTendonSimulationRecord.simulationConfiguration = benchConfig;

        benchRecordDampedFiberElasticTendon = ...
            runMillard2012ComputationalBenchmark(...
              calcDampedFiberElasticTendonMuscleInfoFcn,...
              calcDampedFiberElasticTendonInitialMuscleStateFcn,...
              benchConfig,...
              figDampedFiberElasticBasic,...
              figDampedFiberElasticEnergy,...
              figDampedFiberElasticPower);
        
        dampedFiberElasticTendonSimulationRecord.standardResults = ...
          insertResultsIntoSet( benchRecordDampedFiberElasticTendon,...
                                dampedFiberElasticTendonSimulationRecord.standardResults,...
                                count);

        dampedFiberElasticTendonSimulationRecord.detailedResults.ankleAngularVelocity(:,count)= ...
          benchRecordDampedFiberElasticTendon.pathVelocity./ankleAchillesTendonMomentArm;

        dampedFiberElasticTendonSimulationRecord.detailedResults.ankleAngle(:,count)=...
          + (90 + rampStartAngle)*(pi/180)  ...
          - (dampedFiberElasticTendonSimulationRecord.standardResults.pathLength(:,count) - lpRampStart...
             )*(1/ankleAchillesTendonMomentArm);
          
        
        dampedFiberElasticTendonSimulationRecord.detailedResults.simulationTime(:,count) = ...
          ([0:(1/(npts-1)):1]').*(totalSimTime);

        measurementTime = 0;
        if(flag_measurementSetting == measureAtAnkleAngle)
          measurementTime = measurementTimeAnkleAngle;          
        end           
        if(flag_measurementSetting == measureAtNormFiberLength)
          measurementTime = ...
            getTimeAtNormFiberLength(dampedFiberElasticTendonSimulationRecord, ...
                                     normFiberLengthAtMeasurement,count);          
        end
        dampedFiberElasticTendonSimulationRecord.detailedResults.measurementTime(1,count) = ...
          measurementTime;         
           
        %Get the isometric force of the muscle at the sample point
        if(count ==1)
          normFiberForceAlongTendonIsometric = ...
            interp1(...
              dampedFiberElasticTendonSimulationRecord.standardResults.time(:,count),...
              dampedFiberElasticTendonSimulationRecord.standardResults.normFiberForceAlongTendon(:,count),...
              measurementTime);
        end
        
        dampedFiberElasticTendonSimulationRecord.detailedResults.normFiberForceAlongTendonIsometric(1,count) = ...
          normFiberForceAlongTendonIsometric;           
        
        if(count==countMax)
          disp(['  saved to: ',outputDataFolder, outputFileName,...
            '_dampedFiberElasticTendon_',...
            scaleLceOptStr,...              
            vmaxStr,...
            rampStr,...
            fractionOfFastTwitchFibersStr,...
            tendonStrainAtOneNormForceOverrideStr,...
            tendonDampingStr,...
            ankleAchillesTendonMomentArmStr,...
            measurementSettingStr,'.mat']);          
          save( [outputDataFolder, outputFileName,...
            '_dampedFiberElasticTendon_',...
            scaleLceOptStr,...              
            vmaxStr,...
            rampStr,...
            fractionOfFastTwitchFibersStr,...
            tendonStrainAtOneNormForceOverrideStr,...
            tendonDampingStr,...
            ankleAchillesTendonMomentArmStr,...
            measurementSettingStr,'.mat'],...
                'dampedFiberElasticTendonSimulationRecord');     
        end
              
    end
    

    
end %Velocity loop

