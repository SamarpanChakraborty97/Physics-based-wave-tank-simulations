%make container1 makes a conical container with an opening a the bottom



function s = makeContainer2(s)


%make a fixed container
dx = 0.8*s.free.smoothingLength;

%y-value of wall 1
y1 = [.5:-dx:0]';

%x-values of wall 1
x1 = [0*zeros(length(y1),1)];

x2 = [dx:dx:0.8-dx]';
%x2 = [dx:dx:0.5-dx]';
y2 = zeros(length(x2),1);

y3 = [0:dx:0.5]';
x3 = (x2(end)+dx)*ones(length(y3),1);

x4 = [x3(end)-dx:-dx:dx]';
y4 = 0.5*ones(length(x4),1);

%now put them all together


pos = [[x1 y1];[x2 y2];[x3 y3];[x4 y4]];



nC = length(pos(:,1));

s.comp.nConstrainedParticles = nC;

s.containerParticles.pos   = pos;
s.containerParticles.vel   = zeros(nC,3);
s.containerParticles.mass  = s.free.mass*ones(nC,1);   %mass                         [kg/m^2]
s.containerParticles.smoothingLength = s.free.smoothingLength*ones(nC,1);  %smoothingLength - mean, variance [m]
s.containerParticles.color = s.cons.color*ones(nC,1);


% plot(pos(:,1),pos(:,2),'k.')
% hold on 
% axis equal
% axis([-2 6 -2 2])
 








