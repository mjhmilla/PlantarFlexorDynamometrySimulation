clc;
close all;
clear all;

%%
%Settings
%%
flag_tendonType = 2;
% 1: 4.9%
% 2: 9.2%
% 3: rigid : 

standardTendonElasticity = 0.049;
highTendonElasticity     = 0.092;


%flag_useSubsampleOfTrials = 0;
preloadHauraixReplication = 0.15;
flag_useTendonDamping   = 1;
normalizedTendonDamping = 0.05*(2/3);

%maximumNormalizedFiberVelocity = 6; %in units of norm fiber lengths/second
vmaxStr = ['vmax_10_'];

preloadStr = [num2str(round(preloadHauraixReplication*100,0)),'p'];

if(flag_tendonType ==3)
  flag_useTendonDamping   = 0;  
  normalizedTendonDamping = 0.0; 
end
flag_rampType           = 1;
% 0: Holzer
% 1: Hauraix 2015

scaleLceOpt = 1;
scaleLceOptStr = ['lceOptScale_',num2str(round(scaleLceOpt*100)),'p_'];

switch(flag_rampType)

  case 0
    rampStr = 'ramp_Holzer_';
  case 1
    rampStr = 'ramp_Hauraix2015_';
  otherwise
    asert(0,'flag_rampType must be 0 (Holzer) or 1 (Hauraix 2015)');
end


tendonDampingStr = ['_TendonDamping_',...
        [num2str(round((flag_useTendonDamping*normalizedTendonDamping)*100,0)),'p']]; 

errorBarLineWidth  = 1;
errorBarMarkerSize = 2;
errorBarCapSize    = 0;

dataFolder      = '../data/';
outputFolder    = '../plots/';
dataHauraixFolder= '../dataHauraix2015/';

%%
%Load digitzed experimental data
%%
labelsHauraix2015 = {'w30','w90','w150','w210','w270','w330'};

Hauraix2015Fig3A   = loadDigitizedData('../dataHauraix2015/HauraixFig3A.csv',...
                                    'x','y',labelsHauraix2015,'Fig3A');
Hauraix2015Fig3B   = loadDigitizedData('../dataHauraix2015/HauraixFig3B.csv',...
                                    'x','y',labelsHauraix2015,'Fig3B');
Hauraix2015Fig3C   = loadDigitizedData('../dataHauraix2015/HauraixFig3C.csv',...
                                    'x','y',labelsHauraix2015,'Fig3C');
Hauraix2015Fig3D   = loadDigitizedData('../dataHauraix2015/HauraixFig3D.csv',...
                                    'x','y',labelsHauraix2015,'Fig3D');
Hauraix2015Fig3E   = loadDigitizedData('../dataHauraix2015/HauraixFig3E.csv',...
                                    'x','y',labelsHauraix2015,'Fig3E');
Hauraix2015Fig3F   = loadDigitizedData('../dataHauraix2015/HauraixFig3F.csv',...
                                    'x','y',labelsHauraix2015,'Fig3F');
   
%Map the ankle angle convention to the one used here                                  
scaleFactorYFig3F = 1;
for z=1:1:length(Hauraix2015Fig3A)
  Hauraix2015Fig3A(z).x = 90-(Hauraix2015Fig3A(z).x-90);
  Hauraix2015Fig3B(z).x = 90-(Hauraix2015Fig3B(z).x-90);
  Hauraix2015Fig3C(z).x = 90-(Hauraix2015Fig3C(z).x-90);
  Hauraix2015Fig3D(z).x = 90-(Hauraix2015Fig3D(z).x-90);  
  Hauraix2015Fig3E(z).x = 90-(Hauraix2015Fig3E(z).x-90);    
  Hauraix2015Fig3F(z).x = 90-(Hauraix2015Fig3F(z).x-90);      
  
  if(z==1)
    scaleFactorYFig3F = 75/max(Hauraix2015Fig3F(z).y);
  end
  Hauraix2015Fig3F(z).y = Hauraix2015Fig3F(z).y.*scaleFactorYFig3F;
end


                                  
%%
%Set the configuration of the simulation that you would like to plot
%%

scaleHauraix = 1.0;% 0.8314;

smallMomentArm    = 0.0459;
standardMomentArm = smallMomentArm*1.22;

fractionOfFastTwitchFibers          = 0.5;
ankleAchillesTendonMomentArm        = standardMomentArm;
%measurementSettingStr = '_fixedFiberLength';
measurementSettingStr = 'fixedAnkleAngle';  

preloadStr0   = ['_Preload_',num2str(round(0,0))];
preloadStr50  = ['_Preload_',num2str(round(50,0))];
preloadStr100 = ['_Preload_',num2str(round(100,0))];

%Generate the key words needed to form the name of the corresponding file
fractionOfFastTwitchFibersStr = ...
  ['FastTwitch_',num2str(round(fractionOfFastTwitchFibers*100,0))];

tendonStrainAtOneNormForceStrA = ...
  sprintf('_Tdn_%1.1f',round(standardTendonElasticity*100,1));
tendonStrainAtOneNormForceStrA(1,end-1)='p';
tendonStrainStrA = [sprintf('%1.1f',round(standardTendonElasticity*100,1)),'\%'];


tendonStrainAtOneNormForceStrB = ...
  sprintf('_Tdn_%1.1f',round(highTendonElasticity*100,1));
tendonStrainAtOneNormForceStrB(1,end-1)='p';
tendonStrainStrB = [sprintf('%1.1f',round(highTendonElasticity*100,1)),'\%'];

tendonStrainStrC = '_Tdn_0p0';

ankleAchillesTendonMomentArmStr = sprintf('_%1.1fcm_',...
  round(ankleAchillesTendonMomentArm*100,1));
ankleAchillesTendonMomentArmStr(1,end-4)='p';

tendonLabelA = ['Typical Achilles Tendon ($e^{T}_\circ=$',tendonStrainStrA,')'];
tendonLabelB = ['Compliant Achilles Tendon ($e^{T}_\circ=$',tendonStrainStrB,')']; 
tendonLabelC = ['Rigid Achilles Tendon']; 

tendonStrA = tendonStrainAtOneNormForceStrA;    
tendonStrB = tendonStrainAtOneNormForceStrB;    

tendonStr = '';
tendonStrainStr = '';
switch(flag_tendonType)
  case 1
    tendonStr = tendonStrainAtOneNormForceStrA;    
    tendonStrainStr = tendonLabelA;
  case 2
    tendonStr = tendonStrainAtOneNormForceStrB; 
    tendonStrainStr = tendonLabelB; 
  case 3
    tendonStr = '_Tdn_0p0';
    tendonStrainStr = tendonLabelC; 
        
  otherwise
    assert(0);
end


modelName = ['simFv_gasmed_preload_',preloadStr,'_dampedFiberElasticTendon_'];
if(flag_tendonType==3)
  modelName = ['simFv_gasmed_preload_',preloadStr,'_rigidTendon_'];
end
%%
%Generate the output file names
%%

outputFileName  = ['fig_Hauraix2015Fig3Reproduction_',...
                  modelName,...
                  scaleLceOptStr,...
                  vmaxStr,...
                  rampStr,...
                  fractionOfFastTwitchFibersStr,...
                  tendonStr,...
                  tendonDampingStr,... 
                  ankleAchillesTendonMomentArmStr];

expData = '../data/Own_Study.xlsx'; 

%%
%Generate the names of data sets to load
%%

% simDataSets = {...
%   ['../data/simFv_gasmed_preload_0p_dampedFiberElasticTendon_',...
%     rampStr,...
%     fractionOfFastTwitchFibersStr,...
%     tendonStrA,...
%     tendonDampingStr,... 
%     ankleAchillesTendonMomentArmStr,...
%     measurementSettingStr,'.mat'],...
%   ['../data/simFv_gasmed_preload_50p_dampedFiberElasticTendon_',...
%     rampStr,...  
%     fractionOfFastTwitchFibersStr,...
%     tendonStrA,...
%     tendonDampingStr,... 
%     ankleAchillesTendonMomentArmStr,...
%     measurementSettingStr,'.mat'],...   
%   ['../data/simFv_gasmed_preload_100p_dampedFiberElasticTendon_',...
%     rampStr,...  
%     fractionOfFastTwitchFibersStr,...
%     tendonStrA,...
%     tendonDampingStr,... 
%     ankleAchillesTendonMomentArmStr,...
%     measurementSettingStr,'.mat'],...    
%   ['../data/simFv_gasmed_preload_0p_dampedFiberElasticTendon_',...
%     rampStr,...  
%     fractionOfFastTwitchFibersStr,...
%     tendonStrB,...
%     tendonDampingStr,... 
%     ankleAchillesTendonMomentArmStr,...
%     measurementSettingStr,'.mat'],...
%   ['../data/simFv_gasmed_preload_50p_dampedFiberElasticTendon_',...
%     rampStr,...  
%     fractionOfFastTwitchFibersStr,...
%     tendonStrB,...
%     tendonDampingStr,... 
%     ankleAchillesTendonMomentArmStr,...
%     measurementSettingStr,'.mat'],...   
%   ['../data/simFv_gasmed_preload_100p_dampedFiberElasticTendon_',...
%     rampStr,...  
%     fractionOfFastTwitchFibersStr,...
%     tendonStrB,...
%     tendonDampingStr,... 
%     ankleAchillesTendonMomentArmStr,...
%     measurementSettingStr,'.mat'],...            
%   };



simDataSets = {...
  ['../data/',modelName,...
    scaleLceOptStr,...
    vmaxStr,...
    rampStr,...
    fractionOfFastTwitchFibersStr,...
    tendonStr,...
    tendonDampingStr,... 
    ankleAchillesTendonMomentArmStr,...
    measurementSettingStr,'.mat']
  };



simDataSetColor = [  1, 0, 0];
expDataSetColor = [  0, 0, 1];

simDataSetTag = {'Non-preloaded'};
trialDataSetOmega = zeros(length(simDataSetTag));

            
simDataLineHandles(6) = struct('h',[]);                
                               
 
tendonStrainName = {  ['$$e^{\mathrm{T}}_\circ$$ (',tendonStrainStr,')']};               
                    


load([dataFolder,'normMuscleCurves_',...
                  fractionOfFastTwitchFibersStr,'.mat']);
load([dataFolder,'muscleArch_',...
                  fractionOfFastTwitchFibersStr,'.mat']);

fvCurveSample = calcBezierYFcnXCurveSampleVector(...
                normMuscleCurves.fiberForceVelocityCurve, 500);
              
%%
% Subplot Configuration
%%
rampStartAngle=0;
rampEndAngle=0;
switch(flag_rampType)
  case 0
    rampStartAngle               = -15;   %degrees, -ve: dorsiflexion
    rampEndAngle                 =  15;    %degrees
  case 1
    rampStartAngle               = -20;   %degrees, -ve: dorsiflexion
    rampEndAngle                 =  35;    %degrees
  otherwise
    asert(0,'flag_rampType must be 0 (Holzer) or 1 (Hauraix 2015)');
end

scaleNormForce =  100;

ankleAchillesTendonMomentArm = 0.054; %moment arm in m
%from Rugg et al. & Holzer

omegaMaxDeg = 200;
m2mm        = 1000;
m2cm        = 100;
rad2deg     = (180/pi);
angleMin = (90+rampStartAngle);
angleMax = (90+rampEndAngle);

omegaMin = 0;
omegaMax = 200;

velMaxMM    = (omegaMaxDeg*(pi/180))*ankleAchillesTendonMomentArm*m2mm;
velMaxMM    = velMaxMM*1;




greyColor = [1,1,1].*0.5;

trialLineWidth = [0.5;0.5;0.5;0.5;0.5;0.5;0.5];
expLineWidth = [0.5;0.5;0.5;0.5;0.5;0.5;0.5].*2;
trialLineType = {'-','-','-','-','-','-','-'};

subplotXlabel = {'Ankle Angle ($$^\circ$$)',...
                 'Ankle Angle ($$^\circ$$)',...
                 'Ankle Angle ($$^\circ$$)',...
                 'Ankle Angle ($$^\circ$$)',...
                 'Ankle Angle ($$^\circ$$)',...
                 'Ankle Angle ($$^\circ$$)',...
                 'Angular Velocity ($$^\circ$$/s)'};               
               
subplotYlabel = {'Angular Velocity $$^\circ/\mathrm{s}$$',...
                 'MTU Shortening Velocity (cm/s)',...
                 'Fascicle Shortening Velocity (cm/s)',...
                 'Horizontal Shortening Velocity (cm/s)',...
                 'Tendinous tissue shortening (cm/s)',...
                 'Normalized Torque (Nm/Nm)',...
                 'Fascicle Shortening Velocity (cm/s)'};               

subplotTitle = {'A. ','B.','C.','D.','E.',...
                'F. Hauraix scaled: $$\hat{\tau}(30^\circ/s)$$ is 75\%',...
                'G. '};               
               


subplotXlim = [angleMin*0.95, angleMax*1.05;...
               angleMin*0.95, angleMax*1.05;...
               angleMin*0.95, angleMax*1.05;...
               angleMin*0.95, angleMax*1.05;...
               angleMin*0.95, angleMax*1.05;...
               angleMin*0.95, angleMax*1.05;...
               0, 330];
             
subplotYlim = [-5,340;...
               -5,30;...
               -5,30;...           
               -5,30;...
               -15,20;...
                -5,100;...
                0,25];

angleDegTicks     = [angleMin:10:angleMax];       
forceTicks        = [0:20:100];
omegaDegTicks     = [0,30,90,150,210,270,330];             
velCMTicks        = [0:5:25];             
           
forcePercentTicks = round([0:0.25:1].*100,0);

subplotTicks(12) =struct('xticks',[],'yticks',[]);


subplotTicks(1).xticks = angleDegTicks;
subplotTicks(1).yticks = omegaDegTicks;
subplotTicks(2).xticks = angleDegTicks;
subplotTicks(2).yticks = velCMTicks;
subplotTicks(3).xticks = angleDegTicks;
subplotTicks(3).yticks = velCMTicks;

subplotTicks(4).xticks = angleDegTicks;
subplotTicks(4).yticks = velCMTicks;
subplotTicks(5).xticks = angleDegTicks;
subplotTicks(5).yticks = [-25:5:25];
subplotTicks(6).xticks = angleDegTicks;
subplotTicks(6).yticks = forcePercentTicks;
subplotTicks(7).xticks = [0,30,90,150,210,270,330];
subplotTicks(7).yticks = velCMTicks;

simDataAddLegend = [1,0,0,0,0,0,0];
% simDataLegendPosition = {'SouthEast','NorthEast','SouthWest',...
%                          'NorthWest',...
%                          'SouthEast','SouthWest'};

subPlotIndexDetail = [1,2,3,4,5,6,7];                       
  
%%
% Individual Data Series Configuration
%%

flag_plotExpMeanLine = 0;





%%
% Plot configuration
%%


numberOfFiguresPerPage        = 8;
numberOfVerticalPlotRows      = 4;
numberOfHorizontalPlotColumns = 2;    
assert(numberOfVerticalPlotRows*numberOfHorizontalPlotColumns ...
         >= numberOfFiguresPerPage);

plotHorizMarginCm = 1;
plotVertMarginCm  = 1.5;           
pageHeight  = 29.7;
pageWidth   = 21.0;           
plotHeight  = 7;
plotWidth   = 7;

plotConfigGeneric;

fig_Pub    = figure;     

indexOmegaVsAngle =1;
indexVpVsAngle    =2;
indexVceVsAngle   =3;
indexVceATVsAngle =4;
indexVtVsAngle    =5;
indexTauVsAngle   =6;

indexVceATvsOmega = 7;

%%
%Add the experimental data
%%



%%
%Omega vs Angle
%%
[row,col] = find(subPlotPanelIndex==indexOmegaVsAngle);          
subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
subplot('Position',subPlotVec); 



for indexTrialSet=1:1:length(Hauraix2015Fig3A)

  normZeroToOne = (indexTrialSet-1)/(length(Hauraix2015Fig3A)-1);    
  seriesGrey    = greyColor.*0.75+expDataSetColor.*0.25;
  seriesColor   = seriesGrey.*(1-normZeroToOne) ...
                + expDataSetColor.*(normZeroToOne); 
                
  hdl = plot(  Hauraix2015Fig3A(indexTrialSet).x,...
               Hauraix2015Fig3A(indexTrialSet).y,...
            trialLineType{indexTrialSet},...
            'Color',seriesColor,...
            'LineWidth',expLineWidth(indexTrialSet,1));
  hold on;        

  %if(simDataAddLegend(1,indexTrialSet)==1)
    
  xMid = 105;
  yMid = interp1(Hauraix2015Fig3A(indexTrialSet).x,...
                 Hauraix2015Fig3A(indexTrialSet).y,...
                 xMid);

  %omegaStr = Hauraix2015Fig3A(indexTrialSet).seriesName(1,2:end);
  omegaStr = num2str(round(max( Hauraix2015Fig3A(indexTrialSet).y )));
  
  %plot(xMid,yMid,'.','Color',seriesColor);
  %hold on;

  if(indexTrialSet==length(Hauraix2015Fig3A))
    omegaStr = ['Hauraux 2015: ',omegaStr];
  end
  
  text(xMid,yMid, [omegaStr,'$$^\circ$$/s'],'HorizontalAlignment','Right',...
    'VerticalAlignment','Bottom','Color',seriesColor);
  hold on;
  %end

end    

for indexTrialSet=1:1:length(Hauraix2015Fig3B)
    
    %%
    %MTU velocity vs Angle
    %%
    [row,col] = find(subPlotPanelIndex==indexVpVsAngle);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec); 

    plot(   Hauraix2015Fig3B(indexTrialSet).x,...
            Hauraix2015Fig3B(indexTrialSet).y,...
            trialLineType{indexTrialSet},...
            'Color',seriesColor,...
            'LineWidth',expLineWidth(indexTrialSet,1));

    hold on;  
end

for indexTrialSet=1:1:length(Hauraix2015Fig3C)
    
    %%
    %Fiber velocity vs Angle
    %%
    [row,col] = find(subPlotPanelIndex==indexVceVsAngle);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec); 

    plot(   Hauraix2015Fig3C(indexTrialSet).x,...
            Hauraix2015Fig3C(indexTrialSet).y,...
            trialLineType{indexTrialSet},...
            'Color',seriesColor,...
            'LineWidth',expLineWidth(indexTrialSet,1));

    hold on;  
end

for indexTrialSet=1:1:length(Hauraix2015Fig3D)
    
    %%
    %Fiber velocity along tendon vs Angle
    %%
    [row,col] = find(subPlotPanelIndex==indexVceATVsAngle);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec); 

    plot(   Hauraix2015Fig3D(indexTrialSet).x,...
            Hauraix2015Fig3D(indexTrialSet).y,...
            trialLineType{indexTrialSet},...
            'Color',seriesColor,...
            'LineWidth',expLineWidth(indexTrialSet,1));

    hold on;  
end

for indexTrialSet=1:1:length(Hauraix2015Fig3E)
    
    %%
    %Tendon velocity  vs Angle
    %%
    [row,col] = find(subPlotPanelIndex==indexVtVsAngle);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec); 

    plot(   Hauraix2015Fig3E(indexTrialSet).x,...
            Hauraix2015Fig3E(indexTrialSet).y,...
            trialLineType{indexTrialSet},...
            'Color',seriesColor,...
            'LineWidth',expLineWidth(indexTrialSet,1));

    hold on;  
end




for indexTrialSet=1:1:length(Hauraix2015Fig3F)
    
    %%
    %Tendon velocity  vs Angle
    %%
    [row,col] = find(subPlotPanelIndex==indexTauVsAngle);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec); 

    plot(   Hauraix2015Fig3F(indexTrialSet).x,...
            Hauraix2015Fig3F(indexTrialSet).y,...
            trialLineType{indexTrialSet},...
            'Color',seriesColor,...
            'LineWidth',expLineWidth(indexTrialSet,1));

    hold on;  
end

for indexTrialSet=1:1:length(Hauraix2015Fig3E)
    
    %%
    %Omega vs fiber velocity along tendon
    %%
    [row,col] = find(subPlotPanelIndex==indexVceATvsOmega);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec); 

    angle = 90;
    vceAT = interp1(Hauraix2015Fig3D(indexTrialSet).x,...
                    Hauraix2015Fig3D(indexTrialSet).y,...
                    angle);
    omega = interp1(Hauraix2015Fig3A(indexTrialSet).x,...
                    Hauraix2015Fig3A(indexTrialSet).y,...
                    angle);
    
    plot(   omega,...
            vceAT,...
            'o',...
            'Color',seriesColor,...
            'MarkerFaceColor',seriesColor,...
            'MarkerSize',4);

    hold on;  
end


%%
% Add the simulated data
%%

for indexSim=1:1:length(simDataSets)
  
  figure(fig_Pub);

  %Load the set and remove the simulation-specific name
  data=load(simDataSets{indexSim});
  headField = fields(data);
  data = data.(headField{1});

  totalTrials = size(data.standardResults.activation,2);

  trialSet = [1:1:totalTrials];

  
  
  
  flag_seriesNameAdded = 0;
  
  normFactor = 1.0;
  for indexTrialSet=1:1:length(trialSet)
       
    indexTrial = trialSet(1,indexTrialSet);
    

    
    normZeroToOne = (indexTrialSet-1)/(length(trialSet)-1);    
    seriesGrey    = greyColor.*0.75+simDataSetColor(indexSim,:).*0.25;
    seriesColor   = seriesGrey.*(1-normZeroToOne) ...
                  + simDataSetColor(indexSim,:).*(normZeroToOne); 
                
    %indexSubplot  = indexOmegaVsAngle;        

       

    %%
    %Omega vs Angle
    %%
    [row,col] = find(subPlotPanelIndex==indexOmegaVsAngle);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec); 

    
    hdl = plot(   data.detailedResults.ankleAngle(:,indexTrial).*rad2deg,...
              -data.detailedResults.ankleAngularVelocity(:,indexTrial).*rad2deg,...
              trialLineType{indexTrialSet},...
              'Color',seriesColor,...
              'LineWidth',trialLineWidth(indexTrialSet,1));
    hold on;        

    %if(simDataAddLegend(1,indexTrialSet)==1)
      
    xLbl = NaN;
    yLbl = NaN;
    xErrMin = Inf;
    for z=1:1:length(-data.standardResults.time(:,indexTrial))
      xErr = 110-data.detailedResults.ankleAngle(z,indexTrial).*rad2deg;
      if(abs(xErr) < xErrMin )
        xErrMin=abs(xErr);
        xLbl = data.detailedResults.ankleAngle(z,indexTrial).*rad2deg;
        yLbl = -data.detailedResults.ankleAngularVelocity(z,indexTrial).*rad2deg;
      end
    end
                        
      omegaStr = num2str( round (max(abs(data.detailedResults.ankleAngularVelocity(:,indexTrial)))*180/pi,2));
      
      if(indexTrialSet == length(trialSet))
        omegaStr = ['Sim: ',omegaStr];
      end
      
      text(xLbl,yLbl, [omegaStr,'$$^\circ$$/s'],...
        'HorizontalAlignment','Left',...
        'VerticalAlignment','Bottom',...
        'Color',seriesColor);
      hold on;
    %end
    
    if(indexTrialSet==1 && indexSim == 1)
      text(subplotXlim(1,2),subplotYlim(1,2)*1.2,... 
           tendonStrainStr, 'FontSize',12,...
           'HorizontalAlignment','center');
      hold on;
    end
    
    
    
    %%
    %MTU velocity vs Angle
    %%
    [row,col] = find(subPlotPanelIndex==indexVpVsAngle);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec); 

    plot(   data.detailedResults.ankleAngle(:,indexTrial).*rad2deg,...
            -data.standardResults.pathVelocity(:,indexTrial).*m2cm,...
            trialLineType{indexTrialSet},...
            'Color',seriesColor,...
            'LineWidth',trialLineWidth(indexTrialSet,1));

    hold on;  

    %%
    %Fiber velocity vs Angle
    %%
    [row,col] = find(subPlotPanelIndex==indexVceVsAngle);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec); 

    plot(   data.detailedResults.ankleAngle(:,indexTrial).*rad2deg,...
            -data.standardResults.fiberVelocity(:,indexTrial).*m2cm,...
            trialLineType{indexTrialSet},...
            'Color',seriesColor,...
            'LineWidth',trialLineWidth(indexTrialSet,1));

    hold on;  
    
    
    %%
    %Fiber velocity along tendon vs Angle
    %%    
%     [row,col] = find(subPlotPanelIndex==indexVceATVsAngle);          
%     subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
%     subplot('Position',subPlotVec); 
% 
%     plot(   data.detailedResults.ankleAngle(:,indexTrial).*rad2deg,...
%             -data.standardResults.fiberVelocity(:,indexTrial).*m2cm,...
%             trialLineType{indexTrialSet},...
%             'Color',seriesColor,...
%             'LineWidth',trialLineWidth(indexTrialSet,1));
% 
%     hold on;      
%     
%     
%     %%
%     %lceN vs Angle
%     %%
%     indexOfPlot = subPlotIndexDetail(indexLceNRow,simDataColumn(indexSim));
%     
%     [row,col] = find(subPlotPanelIndex==indexOfPlot);          
%     subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
%     subplot('Position',subPlotVec); 
% 
%     hdl=plot(   data.detailedResults.ankleAngle(:,indexTrial).*rad2deg,...
%              data.standardResults.normFiberLength(:,indexTrial),...
%             trialLineType{indexTrialSet},...
%             'Color',seriesColor,...
%             'LineWidth',trialLineWidth(indexTrialSet,1));
% 
%     hold on;  
%     
%     simDataLineHandles(indexSim).h = [simDataLineHandles(indexSim).h,hdl];
%     here=1;

    %%
    %Tendon velocity vs Angle
    %%
    %indexOfPlot = subPlotIndexDetail(indexVtVsAngleRow,simDataColumn(indexSim));
    
    [row,col] = find(subPlotPanelIndex==indexVtVsAngle);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec); 

    plot(   data.detailedResults.ankleAngle(:,indexTrial).*rad2deg,...
            -data.standardResults.tendonVelocity(:,indexTrial).*m2cm,...
            trialLineType{indexTrialSet},...
            'Color',seriesColor,...
            'LineWidth',trialLineWidth(indexTrialSet,1));

    hold on;  

    %%
    %Path velocity vs Angle
    %%
%     [row,col] = find(subPlotPanelIndex==indexVpVsAngle);          
%     subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
%     subplot('Position',subPlotVec); 
% 
%     plot(   data.detailedResults.ankleAngle(:,indexTrial).*rad2deg,...
%             -data.standardResults.pathVelocity(:,indexTrial).*m2mm,...
%             simDataSetLineType{indexSim},...
%             'Color',seriesColor,...
%             'LineWidth',lineWidth(indexSim,1));
% 
%     hold on;  
    
    %%
    %CE Velocity Along Tendon vs Angle
    %%
    %indexOfPlot = subPlotIndexDetail(indexVceATVsAngleRow,simDataColumn(indexSim));
    
    [row,col] = find(subPlotPanelIndex==indexVceATVsAngle);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec); 

    plot(   data.detailedResults.ankleAngle(:,indexTrial).*rad2deg,...
            -data.standardResults.fiberVelocityAlongTendon(:,indexTrial).*m2cm,...
            trialLineType{indexTrialSet},...
            'Color',seriesColor,...
            'LineWidth',trialLineWidth(indexTrialSet,1));

    hold on;   
    
    %%
    %Norm torque  vs Angle
    %%
    %indexOfPlot = subPlotIndexDetail(indexTauVsAngleRow,simDataColumn(indexSim));
    
    [row,col] = find(subPlotPanelIndex==indexTauVsAngle);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec); 
    
    if(indexTrial ==1)
      normFactor = max(data.standardResults.normFiberForceAlongTendon(:,indexTrial));
    end

    plot(   data.detailedResults.ankleAngle(:,indexTrial).*rad2deg,...
            data.standardResults.normFiberForceAlongTendon(:,indexTrial).*(100/normFactor),...
            trialLineType{indexTrialSet},...
            'Color',seriesColor,...
            'LineWidth',trialLineWidth(indexTrialSet,1));

    hold on;  
    
    %%
    %Omega vs fiber velocity along tendon
    %%
    [row,col] = find(subPlotPanelIndex==indexVceATvsOmega);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec); 

    angle = 90;
    angleErrBest=Inf;
    idxAngleBest = 0;
    for z=1:1:length(data.detailedResults.ankleAngle(:,indexTrial))
      angleErr = abs(angle-data.detailedResults.ankleAngle(z,indexTrial).*rad2deg);
      if(angleErr < angleErrBest)
        angleErrBest=angleErr;
        idxAngleBest=z;
      end
    end
    

    
    plot(   -data.detailedResults.ankleAngularVelocity(idxAngleBest,indexTrial).*rad2deg,...
            -data.standardResults.fiberVelocityAlongTendon(idxAngleBest,indexTrial).*m2cm,...
            'o',...
            'Color',seriesColor,...
            'MarkerFaceColor',seriesColor,...
            'MarkerSize',4);

    hold on;     
    
  end
end



  
for indexRow=1:1:size(subPlotIndexDetail,1)

  for indexColumn = 1:1:size(subPlotIndexDetail,2)
    
    indexOfPlot = subPlotIndexDetail(indexRow,indexColumn);
  
    [row,col] = find(subPlotPanelIndex==indexOfPlot);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec);   
    
    
    
%    if(simDataAddLegend(1,indexOfPlot)==1)
%       omegaLabel1 = [num2str(round(-trialDataSetOmega(1,1),0)),'$$^\circ$$/s'];
%       omegaLabel2 = [num2str(round(-trialDataSetOmega(1,2),0)),'$$^\circ$$/s'];
%       omegaLabel3 = [num2str(round(-trialDataSetOmega(3,3),0)),'$$^\circ$$/s'];
%       omegaLabel4 = [num2str(round(-trialDataSetOmega(4,1),0)),'$$^\circ$$/s'];
%       omegaLabel5 = [num2str(round(-trialDataSetOmega(5,2),0)),'$$^\circ$$/s'];
%       omegaLabel6 = [num2str(round(-trialDataSetOmega(5,3),0)),'$$^\circ$$/s'];
%       
%       legend( [simDataLineHandles(indexColumn).h,simDataLineHandles(indexColumn+3).h],...
%         [tendonStrainName{1},': ',omegaLabel1],...
%         [tendonStrainName{2},': ',omegaLabel2],...
%         [tendonStrainName{3},': ',omegaLabel3],...
%         [tendonStrainName{4},': ',omegaLabel4],...
%         [tendonStrainName{5},': ',omegaLabel5],...
%         [tendonStrainName{6},': ',omegaLabel6],...
%       'Location','SouthEast');      
      %legend;
      %legend boxoff;
%    end


    xticks(subplotTicks(indexOfPlot).xticks);
    yticks(subplotTicks(indexOfPlot).yticks);  

    xlim(subplotXlim(indexOfPlot,:));
    ylim(subplotYlim(indexOfPlot,:));
    
    xlabel(subplotXlabel{indexOfPlot});
    ylabel(subplotYlabel{indexOfPlot});
    title(subplotTitle{indexOfPlot});
    
    box off
  end
end




figure(fig_Pub); 
configPlotExporter;
print('-dpdf',[outputFolder,outputFileName,'.pdf']); 

     
     


              
