clc;
close all;
clear all;

flag_plotMainFigure              = 1;
flag_plotPreloadComparisonFigure = 1;
flag_useTendonDamping            = 1;

%Plot 0% and 100% preload

if(flag_plotMainFigure==1)
  flag_preloadSetting = 1;
  
%   flag_useMaganarisCEArchitecture = 0;
%   for flag_tendonType = 1:1:3  
%     close all;    
%     main_PlotMaxActivationRampShortening;    
%   end
  flag_useMaganarisCEArchitecture = 1;
  for flag_tendonType = 1:1:3     
    close all;    
    main_PlotMaxActivationRampShortening;    
  end
  
end

%Plot 0%, 50% and 100% preload
if(flag_plotPreloadComparisonFigure==1)
  flag_preloadSetting = 2;
  
%   flag_useMaganarisCEArchitecture = 0;
%   for flag_tendonType = 1:1:3    
%      close all;
%      main_PlotMaxActivationRampShortening;    
%   end
  flag_useMaganarisCEArchitecture = 1;
  for flag_tendonType = 1:1:3    
     close all;    
     main_PlotMaxActivationRampShortening;    
  end  
end