% close all
% clear all
% clc
% 
% Velo_Vicon = xlsread('H:\Matlab_Tools\$$_Code\KGR\Own_Study.xlsx','Vicon_angular_Velocity');
% Velo_Isomed = xlsread('H:\Matlab_Tools\$$_Code\KGR\Own_Study.xlsx','Isomed_angular_Velocity');
% Velo_fascicle = xlsread('H:\Matlab_Tools\$$_Code\KGR\Own_Study.xlsx','Velocity_Individual');
% Data_Length = xlsread('H:\Matlab_Tools\$$_Code\KGR\neu_l�nge.csv');
% Preset_Isomed = [0:20:200];
% 
% 
% [fitobject,gof] = fit(mean(-Velo_Vicon')',mean(Velo_fascicle')','poly1')
% 
% n=1;
% for i = 0:180
%     FIT(n) = polyval([fitobject.p1 fitobject.p2],i);
%     XXX (n) = i;
%     n = n + 1;
% end
% 
% txt_fit=['R� = ' num2str(round(gof.rsquare,3))];


% save('H:\Matlab_Tools\$$_Code\KGR\Length_Velocity_Plot.mat');
load('Length_Velocity_Plot.mat');

Optimal_fascicle_length = 42; %(Holzer et al., 2020 - 25� DF (Hessel et al., 2019))


%% Create Plot
figure(1)
subplot(1,2,1)
%e = errorbar(Preset_Isomed,mean(Data_Length),std(Data_Length))
e = errorbar(mean(-Velo_Vicon'),mean(Velo_fascicle'),std(Velo_fascicle'),std(Velo_fascicle'),std(Velo_Vicon'),std(Velo_Vicon'),'o')
hold on
plot(XXX,FIT,'k')
text(40,80,txt_fit,'HorizontalAlignment','Center','Fontsize',20)
text(5,90,'A','HorizontalAlignment','Center','Fontsize',40)

e.Marker = 'o';
e.Color = 'k';
e.LineStyle = 'none';

xticks([0:20:200])
yticks([0:10:100])
axis([-10 210 0 100])

xlabel('Angular Velocity [�/s]')
ylabel('GM Fascicle Velocity [mm/s]')



subplot(1,2,2)
k = errorbar(Preset_Isomed,mean(Data_Length),std(Data_Length))
hold on
text(5,54,'B','HorizontalAlignment','Center','Fontsize',40)

k.Marker = 'o';
k.Color = 'k';
k.LineStyle = '-';

xticks([0:20:200])
yticks([0:5:60])
axis([-10 210 0 60])

xlabel('Preset Dynamometer Velocity [�/s]')
ylabel('GM Fascicle Length at 0� [mm]')
set(gcf,'Position',[117,287,1610,611]);

%%
figure(2)
subplot(1,2,1)
%e = errorbar(Preset_Isomed,mean(Data_Length),std(Data_Length))
e = errorbar(mean(-Velo_Vicon'),...
             mean(Velo_fascicle'),...
             std(Velo_fascicle'),...
             std(Velo_fascicle'),...
             std(Velo_Vicon'),...
             std(Velo_Vicon'),'o')
hold on
plot(XXX,FIT,'k')
text(40,80,txt_fit,'HorizontalAlignment','Center','Fontsize',20)
text(5,100,'A','HorizontalAlignment','Center','Fontsize',40)

e.Marker = 'o';
e.Color = 'k';
e.LineStyle = 'none';

xticks([0:20:200])
yticks([0:10:110])
axis([-10 210 0 110])

xlabel('Angular Velocity [�/s]')
ylabel('GM Fascicle Velocity [mm/s]')



subplot(1,2,2)
k = errorbar(Preset_Isomed,mean(Data_Length)./Optimal_fascicle_length,std(Data_Length)./Optimal_fascicle_length)
hold on
text(5,1,'B','HorizontalAlignment','Center','Fontsize',40)

k.Marker = 'o';
k.Color = 'k';
k.LineStyle = '-';

xticks([0:20:200])
yticks([0:.1:1.1])
axis([-10 210 0 1.1])

xlabel('Preset Dynamometer Velocity [�/s]')
ylabel('Normalized GM Fascicle Length')
set(gcf,'Position',[117,287,1610,611]);

