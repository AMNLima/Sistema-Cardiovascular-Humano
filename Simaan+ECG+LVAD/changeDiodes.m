function [A, B, p] = changeDiodes(Pao, Pae, Pve, E)

if(Pae > Pve)
    Da = 0;
    Dm = 1;
elseif Pve > Pao
    Da = 1;
    Dm = 0;
else
    Da = 0;
    Dm = 0;
end

if Pve > 1
    Rk = 0;
else
    Rk = -3.5*(Pve-1);
end

% Resistances
Rs = 1; % Sistemic Vascular
Rm = 0.0050; % Mitral Valve
Ra = 0.0010; % Aortic Valve
Rc = 0.0398; % Characteristic
Ri = 0.0677; % Inlet Resistance of Cannulae
Ro = 0.0677; % Outlet Resistance of Cannulae

% Compliances
Cr = 4.400; % Left atrial
Cs = 1.330; % Systemic
Ca = 0.080; % Aortic
Cae = Cr;
Cao = Ca;

% Inertance
Ls = 0.0005; % Inertance of blood in Aorta
L = Ls;
Li = 0.0127;
Lo = 0.0127;

% Pressure Difference Parameters
beta0 = 0.17070;
beta1 = 0.02177;
beta2 = -9.3e-5;

V0 = 10; % Volume do ventrículo em ml

%x = [Vve Pao Qa Ps Pae Qb];
A = [ -((Dm/Rm) + (Da/Ra))*E      Da/Ra                     0           0           Dm/Rm                         -1;
      Da/(Ra*Cao)*E               -(Da/(Ra*Cao))            (-1/Cao)    0           0                             1/Cao;
      0                           1/L                       -Rc/L       -1/L        0                             0;
      0                           0                         1/Cs        -1/(Rs*Cs)  1/(Rs*Cs)                     0;
      Dm/(Rm*Cae)*E               0                         0           1/(Rs*Cae)  (-1/Cae)*((1/Rs) + (Dm/Rm))   0;
      E/(Li + Lo + beta1)         -1/(Li + Lo + beta1)      0           0           0                             -(beta0 + Ri + Rk + Ro)/(Li + Lo + beta1)];

%x = [Vve Pao Qa Ps Pae Qb]';
B = [   ((Dm/Rm) + (Da/Ra))*V0;
        (-Da*V0)/(Ra*Cao);
        0;
        0;  
        -(Dm*V0)/(Rm*Cae);
        -V0/(Li + Lo + beta1)];
    
p = [0;
    0;
    0;
    0;
    0;
    -beta2/(Li + Lo + beta1)];
end

