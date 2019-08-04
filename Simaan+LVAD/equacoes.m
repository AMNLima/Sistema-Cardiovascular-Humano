function xdot = equacoes(x,RR,z0)

A_ecg = [1.2, -5, 30, -7.5, 0.75];
B_ecg = [0.25, 0.1, 0.1, 0.1, 0.4];
TH = [-1/3, -1/12, 0, 1/12, 1/2]*pi;
w = 2*pi/RR;
k = 0;
alpha = 1-sqrt(x(1)^2 + x(2)^2);
th = atan2(x(2),x(1));
for ii = 1:5
   k = k - (A_ecg(ii)*(th-TH(ii))*exp(-(th-TH(ii))^2/(2*B_ecg(ii)^2)));
end

%% EQUATIONS
xdot(1) = alpha*x(1)      -w*x(2)    +0*x(3);
xdot(2) =     w*x(1)  +alpha*x(2)    +0*x(3);
xdot(3) =     0*x(1)      +0*x(2)    -1*x(3) + k + z0;














