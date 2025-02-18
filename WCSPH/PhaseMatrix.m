N = 6;
PhaseMat = zeros(N,320);
k = -10+20*rand(1,1);
for i=1:6
    PhaseMat(i,:) = -10+20.*rand(1,320);
    PhaseMat(i,150:170) = k;
end
writematrix(PhaseMat,'Phases.txt','Delimiter',' ')