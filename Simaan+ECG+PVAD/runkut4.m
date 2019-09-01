function [xf, qi, qo] = runkut4(Ts,x,A,B)

% Equação do distema
xdot  = A*x + B;
kx1  = Ts*xdot;
x1  = x  + 0.5*kx1;

xdot = A*x1 + B;
kx2 = Ts*xdot;
x1 = x + 0.5*kx2;

xdot = A*x1 + B;
kx3 = Ts*xdot;
x1 = x + kx3;

xdot = A*x1 + B;
kx4 = Ts*xdot;

q_aux = (kx1 + 2*kx2 + 2*kx3 + kx4)/6;
qi = q_aux(6);
qo = q_aux(7);
xf = x + (kx1 + 2*kx2 + 2*kx3 + kx4)/6;