function xf = runkut4(Ts,x,A,B,p,u,w)

% Equação do distema
xdot  = A*x + B*u + p*w;
kx1  = Ts*xdot;
x1  = x  + 0.5*kx1;

xdot = A*x1 + B*u + p*w;
kx2 = Ts*xdot;
x1 = x + 0.5*kx2;

xdot = A*x1 + B*u + p*w;
kx3 = Ts*xdot;
x1 = x + kx3;

xdot = A*x1 + B*u + p*w;
kx4 = Ts*xdot;

xf = x + (kx1 + 2*kx2 + 2*kx3 + kx4)/6;