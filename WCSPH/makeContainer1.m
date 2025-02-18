%make container1 makes a conical container with an opening a the bottom



function s = makeContainer1(s);



%make a fixed container
dx = 0.5*s.free.smoothingLength;

%y-value of wall 1
y1 = [dx:dx:6]';

%x-values of wall 2
x2 = [0:dx:10]';

%y-values of wall 3
y3 = [dx:dx:6]';


%now put them all together
pos = [[0*ones(length(y1),1) y1];[x2 0*ones(length(x2),1)];[x2(end)*ones(length(y3),1) y3];];



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
 








