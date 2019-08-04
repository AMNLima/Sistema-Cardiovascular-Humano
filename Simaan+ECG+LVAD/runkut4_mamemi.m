function [xdot1, xf] = runkut4(Ts,x,RR,z0)

xdot1 = equacoes(x,RR,z0);
kx1 = Ts*xdot1;
x1 = x + 0.5*kx1';

xdot = equacoes(x1,RR,z0);
kx2 = Ts*xdot;
x1 = x + 0.5*kx2';

xdot = equacoes(x1,RR,z0);
kx3 = Ts*xdot;
x1 = x + kx3';

xdot = equacoes(x1,RR,z0);
kx4 = Ts*xdot;

xf = x + (kx1 + 2*kx2 + 2*kx3 + kx4)'/6;

