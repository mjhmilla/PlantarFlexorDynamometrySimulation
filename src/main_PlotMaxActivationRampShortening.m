clc;
close all;
clear all;

outputFolder = '../plots/';
outputFileName = 'fig_ForceVelocity_Simulation_Vs_Experiment.pdf';


expData = '../data/Own_Study.xlsx'; 



%This was called `path' before: this is also the name of a Matlab's
%search path: do not make variables with the same name. Matlab should 
%give you an error, but it doesn't, it just silently lets you overwrite
%an important internal varible

simDataSets = {'../data/simFv_gaslat_preload_1_dampedFiberElasticTendon.mat',...
              '../data/simFv_gaslat_preload_0_dampedFiberElasticTendon.mat',...               
              '../data/simFv_gaslat_preload_1_rigidTendon.mat',...               
               '../data/simFv_gaslat_preload_0_rigidTendon.mat'};
             
simDataSetName = {'Sim.: ET+Preload',...
                  'Sim.: ET',...                  
                  'Sim.: RT+Preload',...                 
                  'Sim.: RT'};      
                
             
simDataSetColor  = [255, 165, 0;...
                    205, 92, 92;...                    
                    102, 255, 255;...                                          
                    148,  0, 211]./255;
               
simDataSetLineType = {'-','-','-','-'};                  

simDataSetLineWidth= [2.5,0.75,2.5,0.75];

simDataLegendPosition = {'NorthEast','NorthEast','NorthEast'};
%%
% Plot configuration
%%

numberOfFiguresPerPage        = 3;
numberOfVerticalPlotRows      = 3;
numberOfHorizontalPlotColumns = 1;    
assert(numberOfVerticalPlotRows*numberOfHorizontalPlotColumns ...
         >= numberOfFiguresPerPage);

plotHorizMarginCm = 2;
plotVertMarginCm  = 2;           
pageHeight  = 29.7;
pageWidth   = 21.0;           
plotHeight  = 9;
plotWidth   = 16;

plotConfigGeneric;

fig_Fv = figure;

%%
% Tabulated experimental data
%%

Norm_Muscle_Force_Individual = xlsread(expData,'Muscle_Force_Individual');
Vicon_angular_Velocity = xlsread(expData,'Vicon_angular_Velocity');
FSZ_Velocity_Individual = xlsread(expData,'Velocity_Individual');
Isomed_angular_Velocity = xlsread(expData,'Isomed_angular_Velocity');

Chino_Joint_Velo = [0 51.642 98.657 134.254 188.657];                       % Data taken from Chino et al., 2008
Chino_Joint_Velo_STD = [0 5.038 9.739 20.821 16.119];

Chino_force = [1 .868 .768 .697 .603];
Chino_force_STD = [0 .143 .203 .242 .178];

Chino_Fascicle_Velo = [80.168 54.246 42.107 22.128 0];
Chino_Fascicle_Velo = flip(Chino_Fascicle_Velo);
Chino_Fascicle_Velo_STD = [23.52 19.221 12.266 5.564 0];
Chino_Fascicle_Velo_STD = flip(Chino_Fascicle_Velo_STD);

Hauraix_Joint_Velo = [0   29.4442   89.0488  147.9360  203.2319...
                      249.9101  277.9173];                                  % Data taken from Hauraix et al., 2015
Hauraix_Joint_Velo = Hauraix_Joint_Velo.*0.8314;                            % account for the difference in angular ankle joint velocity and dynamometer velocity
Hauraix_Force = [431.76 380.33 262.52 198.54 159.07 139.33 129.77];
Hauraix_Force_norm = Hauraix_Force./Hauraix_Force(1);
Hauraix_Fascicle_Velo = [0 10.61 46.06 72.85 98.53 120.58 131.19];

Hauraix13_Fascicle_Velo = [56 85 113 131 153 168];                          % Data taken from Hauraix et al., 2013
Hauraix13_Joint_Velo = [30 90 150 210 270 330];

figure(fig_Fv);

%%
% R^2 ankle angular velocity vs. fascicle velocity
%%
[row,col] = find(subPlotPanelIndex==1);          
 subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
 subplot('Position',subPlotVec);    

  plot(Chino_Fascicle_Velo,Chino_Joint_Velo,'r','DisplayName','Chino 2008');
  hold on; 
  plot(Hauraix_Fascicle_Velo,Hauraix_Joint_Velo,'b','DisplayName','Hauraix 2015');
  hold on
  plot(Hauraix13_Fascicle_Velo,Hauraix13_Joint_Velo,'--b','DisplayName','Hauraix 13');
  hold on;
 
  m = errorbar( mean(FSZ_Velocity_Individual,2),...
                mean(-Vicon_angular_Velocity,2),...
                std(-Vicon_angular_Velocity,0,2),...
                std(-Vicon_angular_Velocity,0,2),...
                std(FSZ_Velocity_Individual,0,2),...
                std(FSZ_Velocity_Individual,0,2),'ko',...
                'DisplayName','Exp.');  
  hold on
  m.Marker = 'o';
  %    m.Color = 'k';
  %    m.LineStyle = '';    

  [p] = polyfit( mean(FSZ_Velocity_Individual,2),...
                 mean(-Vicon_angular_Velocity,2),1);
  xxx = linspace(min(mean(FSZ_Velocity_Individual,2)),...
                 max(mean(FSZ_Velocity_Individual,2)));

  yyy = polyval(p,xxx);

  ft = fittype( 'poly1' );
  [fit_A,a] = fit(mean(FSZ_Velocity_Individual,2),...
                  mean(-Vicon_angular_Velocity,2),ft);

  plot(xxx,yyy,'k:')
  hold on;

  box off;
  
  
  txt = ['R$^2$ of linear fit:',num2str(a.rsquare)];
  text(10,250,txt,'Interpreter','latex');
  box off
  
  xlabel('Velocity (mm/s)','Interpreter','latex');
  ylabel('Angular Velocity ($^\circ$/s)','Interpreter','latex');   
  title('Fascicle velocity vs. Ankle joint angular velocity');
  
  
%%
% Tendon force vs. fascicle velocity
%%

  [row,col] = find(subPlotPanelIndex==2);          
   subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
   subplot('Position',subPlotVec); 

  hold on
  k = errorbar(Chino_Fascicle_Velo,Chino_force,...
      Chino_force_STD,Chino_force_STD,...
      Chino_Fascicle_Velo_STD,Chino_Fascicle_Velo_STD,...
      'DisplayName','Chino 2008');

  k.Marker = '^';
  k.Color = 'r';
  k.LineStyle = ':';

  plot(Hauraix_Fascicle_Velo,Hauraix_Force_norm,':sb',...
      'DisplayName','Hauraix 2015');

  hold on
  e = errorbar( mean(FSZ_Velocity_Individual,2),...  
                mean(Norm_Muscle_Force_Individual,2),...
                std(Norm_Muscle_Force_Individual,0,2),...
                std(Norm_Muscle_Force_Individual,0,2),...
                std(FSZ_Velocity_Individual,0,2),...
                std(FSZ_Velocity_Individual,0,2),...
                'DisplayName','Exp.');

  e.Marker = 'o';
  e.Color = 'k';
  e.LineStyle = ':';

  ylim([0 1.1])
  box off;
  



  xlabel('Fascicle Velocity (mm/s)')
  ylabel('Norm. Force (\%)')
  title('Fascicle force vs. Fascicle velocity');
  
%%
% Tendon force vs. fascicle velocity
%%

  [row,col] = find(subPlotPanelIndex==3);          
   subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
   subplot('Position',subPlotVec); 


  k = errorbar(Chino_Joint_Velo,Chino_force,...
      Chino_force_STD,Chino_force_STD,...
      Chino_Joint_Velo_STD,Chino_Joint_Velo_STD,...
      'DisplayName', 'Chino 2008');
  hold on

  k.Marker = '^';
  k.Color = 'r';
  k.LineStyle = ':';
  % 
  plot(Hauraix_Joint_Velo,Hauraix_Force_norm,':sb','DisplayName','Hauraix 2015');
  hold on;

  e = errorbar( mean(-Vicon_angular_Velocity,2),...
                mean(Norm_Muscle_Force_Individual,2),...
                std(Norm_Muscle_Force_Individual,0,2),...
                std(Norm_Muscle_Force_Individual,0,2),...
                std(-Vicon_angular_Velocity,0,2),...
                std(-Vicon_angular_Velocity,0,2),...
                'DisplayName','Exp.');

  e.Marker = 'o';
  e.Color = 'k';
  e.LineStyle = ':';    

  ylim([0 1.1]);  
  box off;
  

  
        
  %legend('Hauraix et al., 2015','current study')
  xlabel('Angular velocity ($^\circ$/s)')
  ylabel('Norm. Force (\%)')
  title('Fascicle force vs. Ankle joint angular velocity');
  




%%
% Add the simulated data
%%
  
for i=1:1:length(simDataSets)
  
  %Load the set and remove the simulation-specific name
  data=load(simDataSets{i});
  headField = fields(data);
  data = data.(headField{1});

  trials = size(data.standardResults.activation,2);
  measurementLength                   = zeros(trials,1);
  measuredForce                       = zeros(trials,1);
  measuredFiberVelocity               = zeros(trials,1);
  measuredFiberVelocityAlongTendon    = zeros(trials,1);
  measuredAnkleAngularVelocity        = zeros(trials,1);
  pennationAngle                      = zeros(trials,1);
  
  
  for j=1:1:trials
    %Use interpolation to evaluate the data at the time of the measurement


      measuredForce(j,1) = ...
        interp1(data.detailedResults.simulationTime(:,j), ...
                data.standardResults.normFiberForceAlongTendon(:,j), ...
                data.detailedResults.measurementTime(1,j));

      measuredForce(j,1) = measuredForce(j,1) ...
        ./data.detailedResults.normFiberForceAlongTendonIsometric(1,j);
              
      measuredFiberVelocity(j,1) = ...
        interp1( data.detailedResults.simulationTime(:,j), ...
                 data.standardResults.fiberVelocity(:,j), ...
                 data.detailedResults.measurementTime(1,j));
                                         
      measuredFiberVelocityAlongTendon(j,1) = ...
        interp1(data.detailedResults.simulationTime(:,j), ...
                data.standardResults.fiberVelocityAlongTendon(:,j), ...
                data.detailedResults.measurementTime(1,j));                                         
 
      measuredAnkleAngularVelocity(j,1) = ...
        interp1(data.detailedResults.simulationTime(:,j), ...
                data.detailedResults.ankleAngularVelocity(:,j), ...
                data.detailedResults.measurementTime(1,j));   
                                              
      pennationAngle(j,1) = ...
        interp1(data.detailedResults.simulationTime(:,j), ...
                data.standardResults.pennationAngle(:,j), ...
                data.detailedResults.measurementTime(1,j));   
                                              

  end
 
  %%
  % R^2 ankle angular velocity vs. fascicle velocity
  %%
  m2mm    = 1000;
  rad2deg = 180/pi; 
  
  [row,col] = find(subPlotPanelIndex==1);          
  subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
  subplot('Position',subPlotVec);    
  
    
  plot(measuredFiberVelocity.*(-m2mm), ...
        measuredAnkleAngularVelocity.*(-rad2deg),...
        simDataSetLineType{i},...        
        'Color', simDataSetColor(i,:),...
        'LineWidth',simDataSetLineWidth(1,i),...        
        'DisplayName',simDataSetName{i});
  hold on;
   
    
  %%
  % Tendon force vs. fascicle velocity
  %%

  [row,col] = find(subPlotPanelIndex==2);          
  subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
  subplot('Position',subPlotVec); 
   
  plot(measuredFiberVelocity.*(-m2mm), ...
       measuredForce, ...
       simDataSetLineType{i},...  
        'Color',simDataSetColor(i,:),...        
        'LineWidth',simDataSetLineWidth(1,i),...        
        'DisplayName',simDataSetName{i});
  hold on;


   
  %%
  % Tendon force vs. fascicle velocity
  %%

  [row,col] = find(subPlotPanelIndex==3);          
  subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
  subplot('Position',subPlotVec);    
  
  plot( measuredAnkleAngularVelocity.*(-rad2deg),...
        measuredForce,...
        simDataSetLineType{i},...
        'Color',simDataSetColor(i,:),...        
        'LineWidth',simDataSetLineWidth(1,i),...
        'DisplayName',simDataSetName{i});
  hold on;
  

  
end

for i=1:1:3
  [row,col] = find(subPlotPanelIndex==i);          
  subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
  subplot('Position',subPlotVec);    
  legend('Location',simDataLegendPosition{i});
  %legend boxoff;
end

figure(fig_Fv); 
configPlotExporter;
print('-dpdf',[outputFolder,outputFileName]); 

