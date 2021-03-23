clc;
close all;
clear all;

%%
%Settings
%%
flag_tendonType = 1;
standardTendonElasticity = 0.049;
highTendonElasticity     = 0.14;

flag_useSubsampleOfTrials = 1;



errorBarLineWidth  = 1;
errorBarMarkerSize = 2;
errorBarCapSize    = 0;

dataFolder      = '../data/';
outputFolder    = '../plots/';
%%
%Set the configuration of the simulation that you would like to plot
%%

scaleHauraix = 1.0;% 0.8314;
standardMomentArm = 0.054;
smallMomentArm    = standardMomentArm/1.18;

fractionOfFastTwitchFibers          = 0.5;
ankleAchillesTendonMomentArm        = standardMomentArm;
%measurementSettingStr = '_fixedFiberLength';
measurementSettingStr = '_fixedAnkleAngle';  

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


ankleAchillesTendonMomentArmStr = sprintf('_%1.1fcm_',...
  round(ankleAchillesTendonMomentArm*100,1));
ankleAchillesTendonMomentArmStr(1,end-4)='p';

tendonLabelA = ['Typical Achilles Tendon ($e^{T}_\circ=$',tendonStrainStrA,')'];
tendonLabelB = ['Compliant Achilles Tendon ($e^{T}_\circ=$',tendonStrainStrB,')']; 

tendonStrA = tendonStrainAtOneNormForceStrA;    
tendonStrB = tendonStrainAtOneNormForceStrB;    

% tendonStr = '';
% tendonStrainStr = '';
% tendonLabel = '';
% switch(flag_tendonType)
%   case 1
%     tendonStr = tendonStrainAtOneNormForceStrA;    
%     tendonLabel = ['Typical Achilles Tendon ($e^{T}_\circ=$',tendonLabelA,')'];
%     tendonStrainStr = tendonLabelA;
%   case 2
%     tendonStr = tendonStrainAtOneNormForceStrB;
%     tendonLabel = ['Compliant Achilles Tendon ($e^{T}_\circ=$',tendonLabelB,')'];    
%     tendonStrainStr = tendonLabelB;    
%   otherwise
%     assert(0);
% end

%%
%Generate the output file names
%%

outputFileName  = ['fig_Hauraix2015Fig3Reproduction_',...
                  fractionOfFastTwitchFibersStr,...
                  ankleAchillesTendonMomentArmStr,...
                  measurementSettingStr];

expData = '../data/Own_Study.xlsx'; 

%%
%Generate the names of data sets to load
%%

simDataSets = {...
  ['../data/simFv_gasmed_preload_0p_dampedFiberElasticTendon_',...
    fractionOfFastTwitchFibersStr,...
    tendonStrA,...
    ankleAchillesTendonMomentArmStr,...
    measurementSettingStr,'.mat'],...
  ['../data/simFv_gasmed_preload_50p_dampedFiberElasticTendon_',...
    fractionOfFastTwitchFibersStr,...
    tendonStrA,...
    ankleAchillesTendonMomentArmStr,...
    measurementSettingStr,'.mat'],...   
  ['../data/simFv_gasmed_preload_100p_dampedFiberElasticTendon_',...
    fractionOfFastTwitchFibersStr,...
    tendonStrA,...
    ankleAchillesTendonMomentArmStr,...
    measurementSettingStr,'.mat'],...    
  ['../data/simFv_gasmed_preload_0p_dampedFiberElasticTendon_',...
    fractionOfFastTwitchFibersStr,...
    tendonStrB,...
    ankleAchillesTendonMomentArmStr,...
    measurementSettingStr,'.mat'],...
  ['../data/simFv_gasmed_preload_50p_dampedFiberElasticTendon_',...
    fractionOfFastTwitchFibersStr,...
    tendonStrB,...
    ankleAchillesTendonMomentArmStr,...
    measurementSettingStr,'.mat'],...   
  ['../data/simFv_gasmed_preload_100p_dampedFiberElasticTendon_',...
    fractionOfFastTwitchFibersStr,...
    tendonStrB,...
    ankleAchillesTendonMomentArmStr,...
    measurementSettingStr,'.mat'],...            
  };

simDataColumn = [1,2,3,1,2,3];

simDataSetColor = [  1, 0, 0;...
                     1, 0, 0;...
                     1, 0, 0;...
                     0, 0, 1;...
                     0, 0, 1;...
                     0, 0, 1];

simDataSetTag = {'Slack', '50\% Preload','100\% Preload'};
trialDataSetTag={'Slow','Med', 'Fast'};
trialDataSetOmega = zeros(length(simDataSetTag),3);

            
simDataLineHandles(9) = struct('h',[]);                
                
tendonDataSetName = { tendonLabelA,...
                      tendonLabelA,...
                      tendonLabelA,...
                      tendonLabelB,...
                      tendonLabelB,...
                      tendonLabelB};                   
 
tendonStrainName = {  ['$$e^{\mathrm{T}}_\circ$$ (',tendonStrainStrA,')'],...
                      ['$$e^{\mathrm{T}}_\circ$$ (',tendonStrainStrA,')'],...
                      ['$$e^{\mathrm{T}}_\circ$$ (',tendonStrainStrA,')'],...
                      ['$$e^{\mathrm{T}}_\circ$$ (',tendonStrainStrB,')'],...
                      ['$$e^{\mathrm{T}}_\circ$$ (',tendonStrainStrB,')'],...
                      ['$$e^{\mathrm{T}}_\circ$$ (',tendonStrainStrB,')']};               
                    


load([dataFolder,'normMuscleCurves_',...
                  fractionOfFastTwitchFibersStr,'.mat']);
load([dataFolder,'muscleArch_',...
                  fractionOfFastTwitchFibersStr,'.mat']);

fvCurveSample = calcBezierYFcnXCurveSampleVector(...
                normMuscleCurves.fiberForceVelocityCurve, 500);
              
%%
% Subplot Configuration
%%

scaleNormForce =  100;

ankleAchillesTendonMomentArm = 0.054; %moment arm in m
%from Rugg et al.

omegaMaxDeg = 200;
m2mm        = 1000;
rad2deg     = (180/pi);
angleMin = (90-15);
angleMax = (90+15);

omegaMin = 0;
omegaMax = 200;

velMaxMM    = (omegaMaxDeg*(pi/180))*ankleAchillesTendonMomentArm*m2mm;
velMaxMM    = velMaxMM*1;

indexLceNRow          = 1;
indexVceATVsAngleRow  = 2;
indexVtVsAngleRow     = 3;
indexTauVsAngleRow    = 4;


greyColor = [1,1,1].*0.5;

trialLineWidth = [0.5;0.5;0.5];
trialLineType = {'-','--','-'};

subplotXlabel = {'Ankle Angle ($$^\circ$$)',...
                 'Ankle Angle ($$^\circ$$)',...
                 'Ankle Angle ($$^\circ$$)',...
                 'Ankle Angle ($$^\circ$$)',...
                 'Ankle Angle ($$^\circ$$)',...
                 'Ankle Angle ($$^\circ$$)',...
                 'Ankle Angle ($$^\circ$$)',...
                 'Ankle Angle ($$^\circ$$)',...
                 'Ankle Angle ($$^\circ$$)',...
                 'Ankle Angle ($$^\circ$$)',...
                 'Ankle Angle ($$^\circ$$)',...
                 'Ankle Angle ($$^\circ$$)'};               
               
subplotYlabel = {'Normalized Length ($$\ell^{\mathrm{CE}}/\ell^{\mathrm{CE}}_\circ$$)',...
                 'Normalized Length ($$\ell^{\mathrm{CE}}/\ell^{\mathrm{CE}}_\circ$$)',...
                 'Normalized Length ($$\ell^{\mathrm{CE}}/\ell^{\mathrm{CE}}_\circ$$)',...
                 'Velocity (mm/s)',...
                 'Velocity (mm/s)',...
                 'Velocity (mm/s)',...
                 'Velocity (mm/s)',...
                 'Velocity (mm/s)',...
                 'Velocity (mm/s)',...
                 'Normalized Force (\%)',...
                 'Normalized Force (\%)',...
                 'Normalized Force (\%)'};               

subplotTitle = {'CE Length: 100\% Preload',...
                'CE Length: 50\% Preload',...
                'CE Length: Slack',...
                'CE Velocity (AT): 100\% Preload',...
                'CE Velocity (AT): 50\% Preload',...
                'CE Velocity (AT): Slack',...
                'Tendon Velocity: 100\% Preload',...
                'Tendon Velocity: 50\% Preload',...
                'Tendon Velocity: Slack',...                
                'Ankle Torque: 100\% Preload',...
                'Ankle Torque: 50\% Preload',...
                'Ankle Torque: Slack'};               
               


subplotXlim = [angleMin*0.95, angleMax*1.05;...
               angleMin*0.95, angleMax*1.05;...
               angleMin*0.95, angleMax*1.05;...
               angleMin*0.95, angleMax*1.05;...
               angleMin*0.95, angleMax*1.05;...
               angleMin*0.95, angleMax*1.05;...
               angleMin*0.95, angleMax*1.05;...
               angleMin*0.95, angleMax*1.05;...
               angleMin*0.95, angleMax*1.05;...
               angleMin*0.95, angleMax*1.05;...
               angleMin*0.95, angleMax*1.05;...
               angleMin*0.95, angleMax*1.05];
             
subplotYlim = [0,1.55;...
               0,1.55;...
               0,1.55;...           
               -5,290;...
               -5,290;...
               -5,290;...
               -170,200;...
               -170,200;...
               -170,200;...                   
               -5,105;...
               -5,105;...
               -5,105];

angleDegTicks     = [75, 90, 105];       
vceATTicks        = [0:50:300];
vtTicks           = [-200:50:200];
forceTicks        = [0:20:100];
lceNTicks         = [0:0.25:1.5];
omegaDegTicks     = round([0:0.125:1].*omegaMaxDeg,0);             
velMMTicks        = round([0:0.125:1].*velMaxMM,0);             
           
forcePercentTicks = round([0:0.25:1].*100,0);

subplotTicks(12) =struct('xticks',[],'yticks',[]);


subplotTicks(1).xticks = angleDegTicks;
subplotTicks(1).yticks = lceNTicks;
subplotTicks(2).xticks = angleDegTicks;
subplotTicks(2).yticks = lceNTicks;
subplotTicks(3).xticks = angleDegTicks;
subplotTicks(3).yticks = lceNTicks;

subplotTicks(4).xticks = angleDegTicks;
subplotTicks(4).yticks = round(vceATTicks,0);
subplotTicks(5).xticks = angleDegTicks;
subplotTicks(5).yticks = round(vceATTicks,0);
subplotTicks(6).xticks = angleDegTicks;
subplotTicks(6).yticks = round(vceATTicks,0);

subplotTicks(7).xticks = angleDegTicks;
subplotTicks(7).yticks = round(vtTicks,0);
subplotTicks(8).xticks = angleDegTicks;
subplotTicks(8).yticks = round(vtTicks,0);
subplotTicks(9).xticks = angleDegTicks;
subplotTicks(9).yticks = round(vtTicks,0);


subplotTicks(10).xticks = angleDegTicks;
subplotTicks(10).yticks = forceTicks;
subplotTicks(11).xticks = angleDegTicks;
subplotTicks(11).yticks = forceTicks;
subplotTicks(12).xticks = angleDegTicks;
subplotTicks(12).yticks = forceTicks;

simDataAddLegend = [1,0,0,0,0,0,0,0,0,0,0,0];
% simDataLegendPosition = {'SouthEast','NorthEast','SouthWest',...
%                          'NorthWest',...
%                          'SouthEast','SouthWest'};

subPlotIndexDetail = [3 2 1; 6 5 4; 9 8 7; 12 11 10];                       
  
%%
% Individual Data Series Configuration
%%

flag_plotExpMeanLine = 0;





%%
% Plot configuration
%%


numberOfFiguresPerPage        = 12;
numberOfVerticalPlotRows      = 4;
numberOfHorizontalPlotColumns = 3;    
assert(numberOfVerticalPlotRows*numberOfHorizontalPlotColumns ...
         >= numberOfFiguresPerPage);

plotHorizMarginCm = 1;
plotVertMarginCm  = 1.5;           
pageHeight  = 29.7;
pageWidth   = 21.0;           
plotHeight  = 6;
plotWidth   = 6;

plotConfigGeneric;

fig_Pub    = figure;              

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
  if(flag_useSubsampleOfTrials ==1)
    trialSet = [1, round((totalTrials-1)*0.5) ,totalTrials];
  end
  
  
  
  flag_seriesNameAdded = 0;
  
  normFactor = 1.0;
  for indexTrialSet=1:1:length(trialSet)
       
    indexTrial = trialSet(1,indexTrialSet);
    
    trialDataSetOmega(indexSim, indexTrialSet) = ...
      min(data.detailedResults.ankleAngularVelocity(:,indexTrial).*rad2deg);
    
    normZeroToOne = (indexTrialSet-1)/(length(trialSet)-1);    
    seriesGrey    = greyColor.*0.75+simDataSetColor(indexSim,:).*0.25;
    seriesColor   = seriesGrey.*(1-normZeroToOne) ...
                  + simDataSetColor(indexSim,:).*(normZeroToOne); 
                
    %indexSubplot  = indexOmegaVsAngle;        

       

    %%
    %Omega vs Angle
    %%
%     [row,col] = find(subPlotPanelIndex==indexOmegaVsAngle);          
%     subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
%     subplot('Position',subPlotVec); 
% 
%     
%     hdl = plot(   data.detailedResults.ankleAngle(:,indexTrial).*rad2deg,...
%               -data.detailedResults.ankleAngularVelocity(:,indexTrial).*rad2deg,...
%               simDataSetLineType{indexSim},...
%               'Color',seriesColor,...
%               'LineWidth',lineWidth(indexSim,1));
%     hold on;        
%     if(indexTrial==trialSet(1,end))
%       simDataLineHandles(indexSim).h = hdl;
%     end
            
    

    %%
    %Contractile element velocity vs Angle
    %%
%     [row,col] = find(subPlotPanelIndex==indexVceVsAngle);          
%     subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
%     subplot('Position',subPlotVec); 
% 
%     plot(   data.detailedResults.ankleAngle(:,indexTrial).*rad2deg,...
%             -data.standardResults.fiberVelocity(:,indexTrial).*m2mm,...
%             simDataSetLineType{indexSim},...
%             'Color',seriesColor,...
%             'LineWidth',lineWidth(indexSim,1));
% 
%     hold on;  

    %%
    %lceN vs Angle
    %%
    indexOfPlot = subPlotIndexDetail(indexLceNRow,simDataColumn(indexSim));
    
    [row,col] = find(subPlotPanelIndex==indexOfPlot);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec); 

    hdl=plot(   data.detailedResults.ankleAngle(:,indexTrial).*rad2deg,...
             data.standardResults.normFiberLength(:,indexTrial),...
            trialLineType{indexTrialSet},...
            'Color',seriesColor,...
            'LineWidth',trialLineWidth(indexTrialSet,1));

    hold on;  
    
    simDataLineHandles(indexSim).h = [simDataLineHandles(indexSim).h,hdl];
    here=1;

    %%
    %Tendon velocity vs Angle
    %%
    indexOfPlot = subPlotIndexDetail(indexVtVsAngleRow,simDataColumn(indexSim));
    
    [row,col] = find(subPlotPanelIndex==indexOfPlot);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec); 

    plot(   data.detailedResults.ankleAngle(:,indexTrial).*rad2deg,...
            -data.standardResults.tendonVelocity(:,indexTrial).*m2mm,...
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
    indexOfPlot = subPlotIndexDetail(indexVceATVsAngleRow,simDataColumn(indexSim));
    
    [row,col] = find(subPlotPanelIndex==indexOfPlot);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec); 

    plot(   data.detailedResults.ankleAngle(:,indexTrial).*rad2deg,...
            -data.standardResults.fiberVelocityAlongTendon(:,indexTrial).*m2mm,...
            trialLineType{indexTrialSet},...
            'Color',seriesColor,...
            'LineWidth',trialLineWidth(indexTrialSet,1));

    hold on;   
    
    %%
    %Norm torque  vs Angle
    %%
    indexOfPlot = subPlotIndexDetail(indexTauVsAngleRow,simDataColumn(indexSim));
    
    [row,col] = find(subPlotPanelIndex==indexOfPlot);          
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
    
  end
end



  
for indexRow=1:1:size(subPlotIndexDetail,1)

  for indexColumn = 1:1:size(subPlotIndexDetail,2)
    
    indexOfPlot = subPlotIndexDetail(indexRow,indexColumn);
  
    [row,col] = find(subPlotPanelIndex==indexOfPlot);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subPlotVec);   
    
    
    
    if(simDataAddLegend(1,indexOfPlot)==1)
      omegaLabel1 = [num2str(round(-trialDataSetOmega(indexColumn,1),0)),'$$^\circ$$/s'];
      omegaLabel2 = [num2str(round(-trialDataSetOmega(indexColumn,2),0)),'$$^\circ$$/s'];
      omegaLabel3 = [num2str(round(-trialDataSetOmega(indexColumn,3),0)),'$$^\circ$$/s'];
      omegaLabel4 = [num2str(round(-trialDataSetOmega(indexColumn+3,1),0)),'$$^\circ$$/s'];
      omegaLabel5 = [num2str(round(-trialDataSetOmega(indexColumn+3,2),0)),'$$^\circ$$/s'];
      omegaLabel6 = [num2str(round(-trialDataSetOmega(indexColumn+3,3),0)),'$$^\circ$$/s'];
      
      legend( [simDataLineHandles(indexColumn).h,simDataLineHandles(indexColumn+3).h],...
        [tendonStrainName{1},': ',omegaLabel1],...
        [tendonStrainName{2},': ',omegaLabel2],...
        [tendonStrainName{3},': ',omegaLabel3],...
        [tendonStrainName{4},': ',omegaLabel4],...
        [tendonStrainName{5},': ',omegaLabel5],...
        [tendonStrainName{6},': ',omegaLabel6],...
        'Location','SouthEast');      
      legend boxoff;
    end


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

     
     


              
