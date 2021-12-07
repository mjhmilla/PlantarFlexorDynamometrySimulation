close all;
clear all;
clc;

ankleAngleMVC = -17;

ankleAngle = [-20;-10;0;10;20;30];
pennationAngleA = [30 4; 36 5; 41 5; 45 6; 48 7; 50 7];
pennationAngleB = [32 5; 36 5; 43 6; 43 6; 48 6; 53 7];

fascicleLengthA = [39 6; 34 5; 33 5; 29 4; 28 4; 25 3];
fascicleLengthB = [39 6; 36 5; 31 4; 30 4; 25 4; 24 4];

force = [931, 95; 810, 79; 577, 52; 500, 48; 331, 35; 222,20];


%From Table 1 of
%
%Maganaris CN. Force‐length characteristics of the in vivo human 
%gastrocnemius muscle. Clinical Anatomy: The Official Journal of the 
%American Association of Clinical Anatomists and the British Association 
%of Clinical Anatomists. 2003 May;16(3):215-23.


fig_pennation = figure;

%Pennation angles
lbPennationAngleA = pennationAngleA(:,1) - pennationAngleA(:,2);
ubPennationAngleA = pennationAngleA(:,1) + pennationAngleA(:,2);

lbPennationAngleB = pennationAngleB(:,1) - pennationAngleB(:,2);
ubPennationAngleB = pennationAngleB(:,1) + pennationAngleB(:,2);

meanLBPennationAngle  = (lbPennationAngleA+lbPennationAngleB).*0.5;
meanUBPennationAngle  = (ubPennationAngleA+ubPennationAngleB).*0.5;
meanPennationAngle    = (pennationAngleA(:,1)+pennationAngleB(:,1)).*0.5;

pennationAngleMVC = interp1(ankleAngle,meanPennationAngle,ankleAngleMVC,'linear','extrap');

%Fascicle Lengths
lbFascicleLengthA = fascicleLengthA(:,1)-fascicleLengthA(:,2);
ubFascicleLengthA = fascicleLengthA(:,1)+fascicleLengthA(:,2);

lbFascicleLengthB = fascicleLengthB(:,1)-fascicleLengthB(:,2);
ubFascicleLengthB = fascicleLengthB(:,1)+fascicleLengthB(:,2);

meanLBFascicleLength = (lbFascicleLengthA + lbFascicleLengthB).*0.5;
meanUBFascicleLength = (ubFascicleLengthA + ubFascicleLengthB).*0.5;
meanFascicleLength   = (fascicleLengthA(:,1)+fascicleLengthB(:,1)).*0.5;

pennationAngleMVC = interp1(ankleAngle,meanPennationAngle,ankleAngleMVC,'linear','extrap');

fascicleLengthMVC = interp1(ankleAngle,meanFascicleLength,ankleAngleMVC,'linear','extrap');

forceMVC = interp1(ankleAngle, force(:,1), ankleAngleMVC,'linear','extrap');

fprintf('%1.3f\tPennation angle at MVC\n',pennationAngleMVC);
fprintf('%1.3f\tFascicle length at MVC\n',fascicleLengthMVC);
fprintf('%1.3f\tForce at MVC\n',forceMVC);

subplot(2,2,1);
  
  %fill(bndA(:,1),bndA(:,2),[1,1,1].*0.9+[1,0,0].*0.1,'FaceAlpha',0.5,'EdgeColor','none');
  %hold on
  seriesA=plot(ankleAngle,pennationAngleA(:,1),'r');
  hold on;
  plot(ankleAngle,lbPennationAngleA,'Color',[1,1,1].*0.5+[1,0,0].*0.5);
  hold on;
  plot(ankleAngle,ubPennationAngleA,'Color',[1,1,1].*0.5+[1,0,0].*0.5);
  hold on;

  %fill(bndB(:,1),bndB(:,2),[1,1,1].*0.9+[0,0,1].*0.1,'FaceAlpha',0.5,'EdgeColor','none');
  %hold on
  seriesB=plot(ankleAngle,pennationAngleB(:,1),'b');
  hold on;
  plot(ankleAngle,lbPennationAngleB,'Color',[1,1,1].*0.5+[0,0,1].*0.5);
  hold on;
  plot(ankleAngle,ubPennationAngleB,'Color',[1,1,1].*0.5+[0,0,1].*0.5);
  hold on;

  legend([seriesA,seriesB],{'Location A','Location B'},'Location','SouthEast');
  box off;
  
  xlabel('Ankle Angle');
  ylabel('Pennation Angle (deg)');
  title('Maganaris 2003 Table 1');

subplot(2,2,2);
  seriesC=plot(ankleAngle,meanPennationAngle(:,1),'k');
  hold on;
  plot(ankleAngle,meanLBPennationAngle,'Color',[1,1,1].*0.5);
  hold on;
  plot(ankleAngle,meanUBPennationAngle,'Color',[1,1,1].*0.5);
  hold on;  
  
  seriesMVC=plot(ankleAngleMVC,pennationAngleMVC,'*k');
  hold on;  
  legend([seriesC,seriesMVC],{'Mean(A+B)','MVC'},'Location','SouthEast');
  box off;
  
  xlabel('Ankle Angle');
  ylabel('Pennation Angle (deg)');
  title('Maganaris 2003 Table 1: Mean A and B');
  
subplot(2,2,3);
  seriesFLA = plot(ankleAngle,fascicleLengthA(:,1),'r');
  hold on;
  plot(ankleAngle,lbFascicleLengthA(:,1),'Color',[1,1,1].*0.5+[1,0,0].*0.5);
  hold on;
  plot(ankleAngle,ubFascicleLengthA(:,1),'Color',[1,1,1].*0.5+[1,0,0].*0.5);
  hold on;
  
  seriesFLB = plot(ankleAngle,fascicleLengthB(:,1),'b');
  hold on;
  plot(ankleAngle,lbFascicleLengthB(:,1),'Color',[1,1,1].*0.5+[0,0,1].*0.5);
  hold on;
  plot(ankleAngle,ubFascicleLengthB(:,1),'Color',[1,1,1].*0.5+[0,0,1].*0.5);
  hold on;
  
  legend([seriesFLA,seriesFLB],{'Location A','Location B'},'Location','SouthEast');
  box off;
  
  xlabel('Ankle Angle');
  ylabel('Fascicle Length (mm)');
  title('Maganaris 2003 Table 1');
    
subplot(2,2,4);  
  seriesFLC = plot(ankleAngle, meanFascicleLength(:,1),'k');
  hold on;
  plot(ankleAngle,meanLBFascicleLength,'Color',[1,1,1].*0.5);
  hold on;
  plot(ankleAngle,meanUBFascicleLength,'Color',[1,1,1].*0.5);
  hold on;  
  
  seriesFLMVC = plot(ankleAngleMVC,fascicleLengthMVC,'*k');
  hold on;  
  legend([seriesFLC,seriesFLMVC],{'Mean(A+B)','MVC'},'Location','SouthEast');
  box off;
  
  xlabel('Ankle Angle');
  ylabel('Fascicle Length (mm)');
  title('Maganaris 2003 Table 1: Mean A and B');
  
  
  