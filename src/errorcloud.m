function figH = errorcloud(figH,x,y,yneg,ypos,xneg,xpos,lineType,...
                           fillColor,lineColor,markerFaceColor,lineWidth,...
                           displayName)
                
dataHull=calcDataHull(x,y,yneg,ypos,xneg,xpos);

figure(figH);


for i=1:1:length(dataHull)
  fill(dataHull(i).x,dataHull(i).y,fillColor,'LineStyle','none',...
       'HandleVisibility','off');
  hold on;
end

% plot(x,y, lineType,'Color',lineColor,'MarkerFaceColor',markerFaceColor,...
%           'LineWidth',lineWidth,'DisplayName',displayName);
%         
% hold on;
