format short

F = 4;
L = 40;
h = 0.8;
D = 1.4;
fl = 0.4;
fp = 0.91;
fH = 2.0;
Ab = 0.015;
% kp = 1.97;
% Lp = 3.13;
xb = 14;
% tb = 60;

Lnew = L/F
hnew = h/F
Dnew = D/F
Tp = 1/fp;
Tl = 1/fl;
Th = 1/fH;
Tlnew = Tl/sqrt(F);
Thnew = Th/sqrt(F);
Tpnew = Tp/sqrt(F);
fpnew = 1/Tpnew
flnew = 1/Tlnew
fhnew = 1/Thnew
Abnew = Ab/F
% Lpnew = Lp/F
% kpnew = 2*pi/Lpnew
xbnew = xb/F
% tbnew = tb/sqrt(F)






