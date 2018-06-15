clear ll

lambda=0.780;

r=125/2;    %micron, for Nufern 630-HP ~2micron

N=4;

m=2;    %4 for QWP, 2 for HWP
%R=5;    %cm

%R_micron=R*10*1000;


R=0.133*2*pi*(r^2)/lambda*N*m;

%delta_n=0.133*(r/R_micron)^2;  %for silicon fiber

%delta_phi=delta_n*(2*pi*R_micron)/lambda;
