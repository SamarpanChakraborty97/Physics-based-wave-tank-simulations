function outputTextFile(s)

%output # free particles
%free particle position X
%free particle position Y
%free particle position Z
%free particle mass
%free particle radius
%free parameters

%output # constrained particles
%constrained particle position X
%constrained particle position Y
%constrained particle position Z
%constrained particle mass
%constrained particle radius
%constrained parameters
%constrained particles movement

%computational parameters

cd(s.path.pathToCase);
FID = fopen(s.path.nameOfInputFile,'w');
%FID = 1;

%free particles
nFree = length(s.freeParticles.pos(:,1));
fprintf(FID,'%10.0u\n',nFree);                        %#free
fprintf(FID,'%16.16f\n',s.freeParticles.pos(:,1));    %free x-pos
fprintf(FID,'%16.16f\n',s.freeParticles.pos(:,2));    %free y-pos
fprintf(FID,'%16.16f\n',s.freeParticles.pos(:,3));    %free z-pos
%make a rountine to write density for version 002
fprintf(FID,'%16.16f\n',s.freeParticles.mass(:));     %free radius
fprintf(FID,'%16.16f\n',s.freeParticles.smoothingLength(:));    %free radius
fprintf(FID,'%10.0u\n',7*ones(nFree,1));                      %free color

fprintf(FID,'%16.16f\n',s.freeParticles.vf);          %free particles max velocity
fprintf(FID,'%16.16f\n',s.freeParticles.nu);          %free particles nu
fprintf(FID,'%16.16f\n',s.freeParticles.delta);          %free particles delta-SPH term 
fprintf(FID,'%16.16f\n',s.freeParticles.muAlpha);          %free particles alpha viscosity
fprintf(FID,'%16.16f\n',s.freeParticles.muBeta);      %free particles beta viscosity
fprintf(FID,'%16.16f\n',s.freeParticles.epsilon);     %free particles epsilon
fprintf(FID,'%16.16f\n',s.freeParticles.rRef);        %reference density


%constrained particles
nContainer = length(s.containerParticles.pos(:,1));
fprintf(FID,'%10.0u\n',nContainer);                        %#container
fprintf(FID,'%16.16f\n',s.containerParticles.pos(:,1)); %container x-pos
fprintf(FID,'%16.16f\n',s.containerParticles.pos(:,2)); %container y-pos
fprintf(FID,'%16.16f\n',s.containerParticles.pos(:,2)); %container z-pos
%make a rountine to write density for version 002
fprintf(FID,'%16.16f\n',s.containerParticles.mass(:));    %container mass
fprintf(FID,'%16.16f\n',s.containerParticles.smoothingLength(:));    %container smoothingLengh
fprintf(FID,'%10.0u\n',2*ones(nContainer,1));    %constrained color
%fprintf(FID,'%10.0u\n',nFree+(1:nContainer)');    %constrained color


% fprintf(FID,'%10u\n',3);                      %Xfunction 3 - p1*sin(2*pi*t*p2+p3)
% %fprintf(FID,'%10u\n',2);
% fprintf(FID,'%10u\n',1);                      %p1
% fprintf(FID,'%10u\n',1.5);                      %p2
% fprintf(FID,'%10u\n',0);                      %p3
% fprintf(FID,'%10u\n',10);                      %Yfunction 10 - do nothing
% fprintf(FID,'%10u\n',2);                        %range1
% fprintf(FID,'%10u\n',312);                      %range2
% 
% fprintf(FID,'%10u\n',3);                       %Xfunction 3 - p1*sin(2*pi*t*p2+p3)
% fprintf(FID,'%10u\n',1);                      %p1
% fprintf(FID,'%10u\n',1.5);                      %p2
% fprintf(FID,'%10u\n',pi);                    %p3
% fprintf(FID,'%10u\n',10);                      %Yfunction 10 - do nothing
% fprintf(FID,'%10u\n',566);                       %range1
% fprintf(FID,'%10u\n',695);                      %range2
% 
fprintf(FID,'%10u\n',10);                      %Xfunction 10 - do nothing - takes no parameters
fprintf(FID,'%10u\n',10);                      %Yfunction 10 - do nothing
fprintf(FID,'%10u\n',10);                      %Zfunction 10 - do nothing
fprintf(FID,'%10u\n',1);                       %range1
fprintf(FID,'%10u\n',10);                      %range2
% %

% %fprintf(FID,'%10u\n',3);                       %Xfunction 3 - p1*sin(2*pi*t*p2+p3)
% fprintf(FID,'%10u\n',10);                       %X function 10 - do nothing
% fprintf(FID,'%10u\n',10);                       %Y function 10 - do nothing
% %fprintf(FID,'%10u\n',0.8);                      %p2
% %fprintf(FID,'%10u\n',0);                      %p3
% fprintf(FID,'%10u\n',1);                        %Z function 1 - unit step
% fprintf(FID,'%10u\n',0.1900);
% fprintf(FID,'%10u\n',0.3000);
% %fprintf(FID,'%10u\n',2);
% %fprintf(FID,'%10u\n',312);
% fprintf(FID,'%10u\n',3340);                        %range1
% fprintf(FID,'%10u\n',3738);                      %range2


%computational parameters
fprintf(FID,'%10u\n',0);                      %no more functions
fprintf(FID,'%16.16f\n',s.comp.gravity);      %gravity
fprintf(FID,'%10.0u\n',s.comp.nT);            %# time steps
fprintf(FID,'%16.16f\n',s.comp.dt);           %dt
fprintf(FID,'%10.0u\n',s.comp.storageStride); %storage Stride


fclose(FID);