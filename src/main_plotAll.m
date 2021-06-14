clc;
close all;
clear all;

%Plot 0% and 100% preload
flag_preloadSetting = 1;
for flag_tendonType = 1:1:3      
  main_PlotMaxActivationRampShortening;    
end

%Plot 0%, 50% and 100% preload

% flag_preloadSetting = 2;
% for flag_tendonType = 1:1:3    
%   main_PlotMaxActivationRampShortening;    
% end
