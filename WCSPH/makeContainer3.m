%make container1 makes a conical container with an opening a the bottom



function s = makeContainer3(s);

%tempWallx = s.free.xSpan(2);

%make a fixed container
dx = 5 * 10^(-3);
%h = 2 * dx;

length1 = 7.036;
%depth = 0.4;

dxc = 0.8 * dx;

%y-value of wall 1
y1 = [0.35:-dxc:0]';
%x-values of wall 1
x1 = -0.5 * ones(length(y1),1);

x2 = [-0.5+dxc:dxc:length1]';
%x2 = [dx:dx:0.5-dx]';
y2 = zeros(length(x2),1);

x3 = [length1+dxc : dxc: length1 + 3.4]';
x31 = [dxc : dxc: 3.4]';
y3 = tand(3.75) .* x31;

% y5 = [0.4:-dxc:dxc]';
% x5 = zeros(length(y5),1);
% 
% y6 = [4-dxc:-dxc:1.6+dxc]';
% x6 = zeros(length(y6),1);

% x7 = [-0.5+dxc:dxc:length1+8]';
% y7 = 4*ones(length(x7),1);
% % 
% y8 = [y3(end)+dxc:dxc:4-dxc]';
% x8 = x3(end)*ones(length(y8),1);

y4 = [0.35:-dxc:0]';
x4 = zeros(length(y4),1);
% 
% y5 = [1.2-dxc:-dxc:dxc]';
% x5 = (tempWallx+2*dx)*ones(length(y5),1);
% %now put them all together


pos = [[x1 y1];[x2 y2];[x3 y3];[x4 y4]];
pos2 = [[x1 y1];[x2 y2];[x3 y3]];

% r2=find(pos(:,1)==tempWallx+2*dx & pos(:,2)==dxc)
% r1 = find(pos(:,1)==tempWallx+2*dx & pos(:,2)==1.2-dxc)



nC = length(pos(:,1));

Piston_particle_start = length(pos2(:,1))+1
Piston_particle_end = length(pos(:,1))

s.comp.nConstrainedParticles = nC;

s.containerParticles.pos   = pos;
s.containerParticles.vel   = zeros(nC,3);
s.containerParticles.mass  = s.free.mass*ones(nC,1);   %mass                         [kg/m^2]
s.containerParticles.smoothingLength = s.free.smoothingLength*ones(nC,1);  %smoothingLength - mean, variance [m]
s.containerParticles.color = s.cons.color*ones(nC,1);


% plot(pos(:,1),pos(:,2),'k.')
% hold on 
% axis equal
% axis([-0.2 3.5 -0.1 2.1])