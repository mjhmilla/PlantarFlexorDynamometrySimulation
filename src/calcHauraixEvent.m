function [value,isterminal,direction] = calcHauraixEvent(t,x,endAngle)

value =0;
isterminal = 1;
direction = 1;
if(x >= endAngle)
  value=1;
end