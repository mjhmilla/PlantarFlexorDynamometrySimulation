clc;
close all;
clear all;


outputDataFolder = '../data/';

outputFolder = '../plots/';
outputFileName = 'fig_IsometricNormMusculotendonForce_ExpAndSim.pdf';

fractionOfFastTwitchFibers         = 0.5;
% 0: Gives you a force-velocity curve constent with a slow-twitch fiber
% 1: Gives you a force-velocity curve constent with a fast-twitch fiber  

fractionOfFastTwitchFibersStr = ...
  ['FastTwitch_',num2str(round(fractionOfFastTwitchFibers*100,0))];

flag_generateDiagnosticPlots = 0;

%%
% Plot configuration
%%

numberOfFiguresPerPage        = 9;
numberOfVerticalPlotRows      = 3;
numberOfHorizontalPlotColumns = 3;    
assert(numberOfVerticalPlotRows*numberOfHorizontalPlotColumns ...
         >= numberOfFiguresPerPage);

plotHorizMarginCm = 1.5;
plotVertMarginCm  = 1.5;           
pageHeight  = 29.7;
pageWidth   = 21.0;           
plotHeight  = 6;
plotWidth   = 6;

plotConfigGeneric;
%%
% Plot the experimental data
%%
expData = load('../data/Isos_for_Matt.mat');

figData=figure;

dataFields  = {'MVC_100_01',...
               'MVC_100_02',...
               'MVC_105_01',...
               'MVC_105_02',...
               'MVC_110_01',...
               'MVC_110_02',...
               'MVC_115_01',...
               'MVC_115_02',...
               'MVC_120_01',...
               'MVC_120_02',...
                'MVC_70_01',...
                'MVC_70_02',...
                'MVC_80_01',...
                'MVC_80_02',...
                'MVC_90_01',...
                'MVC_90_02',...
                'MVC_95_01',...
                'MVC_95_02'};

isometricAnkleAngles = [100, 105, 110,...
                        115, 120, 70,...
                         80,  90, 95];              
              
pIdx = Inf;              
for i=1:1:length(dataFields)
  idx = floor((i-1)/2)+1;
  
  [row,col] = find(subPlotPanelIndex==idx);          
  subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
  subplot('Position',subPlotVec);
       
  n = 0.25;
  if(pIdx ==idx)
    n=0.75;

    tqMax = max(expData.torque.(dataFields{i}));
    tqMin = mean(expData.torque.(dataFields{i})(1:100,1));
    plot(expData.timestamp_analogs.(dataFields{i}),...
         (expData.torque.(dataFields{i})-tqMin)./(tqMax-tqMin),...
         'Color',[1,1,1].*n, 'LineWidth',2);
    hold on;
    legendText = dataFields{i};

    if(pIdx == idx)
      if(idx==1)
        title(['Ankle Angle $(>90^\circ \,DF): ', legendText(1,5:(end-3)),'^\circ$']);        
      else
        title(['Ankle Angle: $', legendText(1,5:(end-3)),'^\circ$']);        
      end
    end
  end
  xlabel('Time (s)');
  ylabel('Norm. 0-1 Torque (Nm)');
  xlim([0,8]);
  ylim([-0.1,1.1]);
  box off;
  grid on;
  pIdx = idx;
end

%% 
% Simulate the same thing
%   Modified code from main_MaxActivationRampShortening20201002.m
%%

tendonStrainAtOneNormForceOverride = 0.1;
muscle          = 'gasmed';
initialHoldTime = 1;

flag_runRigidBench               = 0;
flag_runClassicElasticBench      = 0;
flag_runDampedFiberElasticBench  = 1;
flag_updateExistingPlots         = 0;

%%
%Configure the muscle
%%


countMax = length(isometricAnkleAngles);

colorStart = [0,0,1];
colorEnd   = [1,0,0];
colorRecord = zeros(countMax,3);

npts=100*8;

for count=1:1:countMax

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
  
  

  ankleAngularVelocity(count)  =  0;
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

  shiftFiberForceLengthCurve = 0.25;
  
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
        fractionOfFastTwitchFibers,...
        fractionOfFastTwitchFibersStr,...
        shiftFiberForceLengthCurve,...
        outputDataFolder,...        
        flag_updateNormMuscleCurves,...
        flag_plotNormMuscleCurves);

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


  benchConfig.npts   = npts;
  benchConfig.nSim   = 1;
  benchConfig.relTol = 1e-5;
  benchConfig.absTol = 1e-5;
  benchConfig.initialActivationState = 0;%[0:0.1:1]';
  benchConfig.color = trialColor;  

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

  ankleAngle = -(isometricAnkleAngles(1,count)-90);
  
  
  lpDeltaStart = -(ankleAngle-ankleAngleMaxPlantarFlexion ...
      )*(pi/180)*ankleAchillesTendonMomentArm;

  lpRampStart = lpOpt + lpDeltaStart;

  %%
  % II. Change in Path Length: this is an isometric test so there is no
  % change
  %%
  lpRampEnd = lpRampStart;
  lpRampMid = lpRampEnd;
  
  %%
  % III. Time required to complete the ramp
  %%

  initialHoldTime = 0;
  rampTime        = 8;
  totalSimTime    = 8;
  minHoldTime     = 0;

  measurementTime = rampTime*0.5;
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

  expIdx = count*2;
  
  %first second of data is quiet.
  tqStart = mean(expData.torque.(dataFields{expIdx})(1:1000,1));
  
  flagFoundStart= 0;
  flagFoundEnd = 0;
  
  [minTq, idxMinTq] = min(expData.torque.(dataFields{expIdx}));  
  [maxTq, idxMaxTq] = max(expData.torque.(dataFields{expIdx}));
  tqThreshold = 0.05*(maxTq-minTq);
  
  
  
  for z=1:1:length(expData.torque.(dataFields{expIdx}))
    if (flagFoundStart == 0 && ...
        expData.torque.(dataFields{expIdx})(z,1) > (tqStart+tqThreshold))
      flagFoundStart=1;
      ton=expData.timestamp_analogs.(dataFields{expIdx})(1,z);
    end
    if(z > idxMaxTq && flagFoundEnd == 0 ...
      && (maxTq-expData.torque.(dataFields{expIdx})(z,1)) > tqThreshold)
      flagFoundEnd = 1;
      toff=expData.timestamp_analogs.(dataFields{expIdx})(1,z);
    end
  end
  

  excitationFcn = @(argTime)calcStepFunction(argTime, ton, toff);    
  activationFcn = @(argEx, argAct)calcFirstOrderActivationDerivative(...
                                  argEx, argAct, tact,tdeact,tmin);

  benchConfig.excitationFcn = excitationFcn;
  benchConfig.activationFcn = activationFcn;
  
  simTime = [ benchConfig.tspan(1):...
              (benchConfig.tspan(2)-benchConfig.tspan(1))/(npts-1): ...
              benchConfig.tspan(2)];
  
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
          [],[],[]);

      normFiberForceAT = benchRecordRigid.normFiberForceAlongTendon(:,1);
      normForceMin = min(normFiberForceAT);
      normForceMax = max(normFiberForceAT);
      normForceAT01 = (normFiberForceAT-normForceMin)./(normForceMax-normForceMin);
      
      figure(figData)
      
      [row,col] = find(subPlotPanelIndex==count);          
      subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
      subplot('Position',subPlotVec);

      plot(simTime,normForceAT01,'r');
      hold on;        
        
      %save('benchRecordRigid.mat','benchRecordRigid');
      %saveas(figRigidTendonBasic,  'figRigidTendonBasic.fig','fig');
      %saveas(figRigidTendonEnergy, 'figRigidTendonEnergy.fig','fig');
      %saveas(figRigidTendonPower,  'figRigidTendonPower.fig','fig');

        
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
          benchConfig,[],[],[]);


      normFiberForceAT = benchRecordClassicElastic.normFiberForceAlongTendon(:,1);
      normForceMin = min(normFiberForceAT);
      normForceMax = max(normFiberForceAT);
      normForceAT01 = (normFiberForceAT-normForceMin)./(normForceMax-normForceMin);
      
      
      figure(figData)
      
      [row,col] = find(subPlotPanelIndex==count);          
      subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
      subplot('Position',subPlotVec);     
      
      plot(simTime,normForceAT01,'g');
      hold on;

        
        
      %save('benchRecordClassicElastic.mat','benchRecordClassicElastic');
      %saveas(figClassicElasticBasic,'figClassicElasticBasic.fig','fig');
      %saveas(figClassicElasticEnergy,'figClassicElasticEnergy.fig','fig');
      %saveas(figClassicElasticPower,'figClassicElasticPower.fig','fig');

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
          benchConfig,[],[],[]);

        
      normFiberForceAT = benchRecordDampedFiberElasticTendon.normFiberForceAlongTendon(:,1);
      normForceMin = min(normFiberForceAT);
      normForceMax = max(normFiberForceAT);
      normForceAT01 = (normFiberForceAT-normForceMin)./(normForceMax-normForceMin);      
      
      
      figure(figData)
      
      [row,col] = find(subPlotPanelIndex==count);          
      subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
      subplot('Position',subPlotVec);
      
      plot(simTime,normForceAT01,'b');
      hold on;
        
  end
  
  
end

figure(figData); 
configPlotExporter;
print('-dpdf',[outputFolder,outputFileName]); 

