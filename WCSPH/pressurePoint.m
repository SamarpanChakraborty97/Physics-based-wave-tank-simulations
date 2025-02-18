% create a particle for measurement of pressure at a specified location

function s = pressurePoint(s)

height = 6 * 10^(-1); % height of water column
dx = height/50; % particle spacing = height/50
h = 2 * dx;
% width = 2* height;
% dxc = 0.2 * dx;


x = 5.366 * height - 2 * h;
y = 0.19 * height;
pos = [x y];

nMeasure = length(pos(:,1));
s.comp.nMeasuredLocations = nMeasure;
s.measuredLocations.pos   = pos;
s.measuredLocations.vel   = zeros(nMeasure,3);
s.measuredLocations.mass  = zeros(nMeasure,1);   %mass                         [kg/m^2]
s.measuredLocations.smoothingLength = s.free.smoothingLength*ones(nMeasure,1);  %smoothingLength - mean, variance [m]
s.measuredLocations.color = s.measured.color*ones(nMeasure,1);
