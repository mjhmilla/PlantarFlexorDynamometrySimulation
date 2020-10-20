function dataHull=calcDataHull(x,y,yneg,ypos,xneg,xpos)
           


ellipseQuarterPts = 25;
ellipsePts = 4*ellipseQuarterPts;
                      

angles = [0:(1/(ellipsePts-1)):1]'.*(2*pi);

circlePts = [cos(angles),sin(angles)];

dataHull(length(x)-1) = struct('x',[],'y',[]);

%Add the points that define the ellipse that encloses the data
idx = 1;
idxB = 1;

dataEllipse = zeros(ellipsePts,2);
dataEllipsePrev = zeros(ellipsePts,2);

flag_debug=0;
if(flag_debug==1)
  test=figure;
end

for i=1:1:length(x)
  
  dataEllipsePrev = dataEllipse;
  
  idx = 1;%(i-1)*ellipsePts + 1;
  %First quarter
  idxA = idx;
  idxB = idxA+ellipseQuarterPts-1;
  idxCA = 1;
  idxCB = idxCA+ellipseQuarterPts-1;
  dataEllipse(idxA:idxB,1) = x(i,1) + circlePts(idxCA:idxCB,1).*abs(xpos(i,1));
  dataEllipse(idxA:idxB,2) = y(i,1) + circlePts(idxCA:idxCB,2).*abs(ypos(i,1));
  
  %Second quarter
  idxA = idxB+1;
  idxB = idxA+ellipseQuarterPts-1;
  idxCA = idxCB+1;
  idxCB = idxCA+ellipseQuarterPts-1;  
  dataEllipse(idxA:idxB,1) = x(i,1) + circlePts(idxCA:idxCB,1).*abs(xneg(i,1));  
  dataEllipse(idxA:idxB,2) = y(i,1) + circlePts(idxCA:idxCB,2).*abs(ypos(i,1));  

  %Third quarter
  idxA = idxB+1;
  idxB = idxA+ellipseQuarterPts-1;
  idxCA = idxCB+1;
  idxCB = idxCA+ellipseQuarterPts-1;    
  dataEllipse(idxA:idxB,1) = x(i,1) + circlePts(idxCA:idxCB,1).*abs(xneg(i,1));  
  dataEllipse(idxA:idxB,2) = y(i,1) + circlePts(idxCA:idxCB,2).*abs(yneg(i,1));  

  %Fourth quarter
  idxA = idxB+1;
  idxB = idxA+ellipseQuarterPts-1;
  idxCA = idxCB+1;
  idxCB = idxCA+ellipseQuarterPts-1;    
  dataEllipse(idxA:idxB,1) = x(i,1) + circlePts(idxCA:idxCB,1).*abs(xpos(i,1));  
  dataEllipse(idxA:idxB,2) = y(i,1) + circlePts(idxCA:idxCB,2).*abs(yneg(i,1));  
     
  if(i > 1)
    dataPair=[dataEllipsePrev;dataEllipse];
    
    k = convhull([dataPair(:,1)],...
                 [dataPair(:,2)]);
                 
    dataHull(i-1).x = dataPair(k,1);
    dataHull(i-1).y = dataPair(k,2);
    
    dataHull(i,1).xCenter = [x(i-1,1);x(i,1)];
    dataHull(i,1).yCenter = [y(i-1,1);y(i,1)];
    

    
    if(flag_debug==1)
      figure(test);
      clf;
      fill(dataHull(i-1).x,dataHull(i-1).y,[1,1,1].*0.5);
      hold on;
      plot(dataEllipsePrev(:,1),dataEllipsePrev(:,2),'c');
      hold on;
      plot(dataEllipse(:,1),dataEllipse(:,2),'b');
      hold on;
      plot([x(i,1)],[y(i,1)],'ok');
      hold on;
      plot([x(i,1)-xneg(i,1);x(i,1)+xpos(i,1)],[y(i,1)],'-k');
      hold on;
      plot([x(i,1)-xneg(i,1);x(i,1)+xpos(i,1)],[y(i,1);y(i,1)],'-k');
      hold on;
      plot([x(i,1);x(i,1)],[y(i,1)-yneg(i,1);y(i,1)+ypos(i,1)],'-k');
      hold on;
      xlabel('X');
      ylabel('Y');
      
    end
  end
  
end