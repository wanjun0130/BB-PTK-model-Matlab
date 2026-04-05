% this is the full ODEs for 8 compartment lung model
% Jingyu Yu @ 12/20/2008
function [dy] = LungModelOde(t,x,m,g)
dy = m*x+g;
end
