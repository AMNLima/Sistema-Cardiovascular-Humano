function [A, B] = changeDiodes(Pao, Pae, Pve, E, Pc, alpha, Pex, Vve, Qo)

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

if(Pve > Pc)
        Di = 1;
        Do = 0.03;
    elseif (Pc > Pao)
        Di = 0.08;
        Do = 1;
    else
        Di = 0.08;
        Do = 0.03;
end


% Resistances
Rs = 0.8738; % Sistemic Vascular
Rm = 0.0050; % Mitral Valve
Ra = 0.0010; % Aortic Valve
Rc = 0.0398; % Characteristic
Ri = 0.15; % Inlet Resistance of Cannulae
Rp = 0.05;
Rd = 0.01;

% Compliances
Cr = 4.000; % Left atrial
Cs = 2.896; % Systemic
Ca = 0.080; % Aortic
Cae = Cr;
Cao = Ca;
Cp = 2;
Cd = 4;

% Inertance
Ls = 0.001025; % Inertance of blood in Aorta
L = Ls;
Li = 0.0854;
Lo = 0.0087;
Lp = 0.0033;

% Pressure Difference Parameters
% beta0 = 0.17070;
% beta1 = 0.02177;
% beta2 = -9.3e-5;
betai = Di/(Li + Di * Lp);
betao = Do/(Lo + Do * Lp);
gamma = betai * Lp * betao;

Ro = 0.05 + 0.00015 * abs(Qo);
Ri_ = Ri + exp(-0.25*Vve);
omegai = Ri_ / (Li + Di * Lp);
omegao = Ro / (Lo + Do * Lp);

V0 = 15; % Volume do ventrículo em ml
Vd = 107;

phi = (1-gamma*Lp);

% x = [Vve, Pae, Qa, Pao, Ps, Qi, Qo, Pd, Vc]
A = [-(Dm/Rm + Da/Ra) * E, Dm/Rm, 0, Da/Ra, 0, -1, 0, 0, 0;
    (Dm * E)/(Rm * Cae), -1/(Cae) * (1/Rs + Dm/Rm), 0, 0, 1/(Rs * Cae), 0, 0, 0, 0;
    0, 0, -Rc/Ls, 1/Ls, -1/Ls, 0, 0, 0, 0;
    (Da * E)/(Ra * Cao), 0, -1/Cao, -Da/(Ra * Cao), 0, 0, 1/Cao, 0, 0;
    0, 1/(Rs * Cs), 1/Cs, 0, -1/(Rs * Cs), 0, 0, 0, 0;
    betai * E / phi, 0, 0, - gamma / phi, 0, - (alpha * (betai - gamma) + (betai * Rp + omegai + gamma * Rp))/(phi), (alpha * (betai - gamma) + (betai * Rp - gamma * Rp - betai * Lp * omegao))/(phi), -(betai - gamma)/(phi), -(betai - gamma)/(Cp * (phi));
    0, 0, 0, -betao / (phi), 0, (alpha * (betao - gamma) + (betao * Rp - betao * Lp * omegai - gamma * Rp))/(1 - gamma * Lp), -(alpha * (betao - gamma) + (betao * Rp + omegao - gamma * Rp))/(phi), (betao - gamma) / (1 - gamma * Lp), (betao - gamma)/(Cp * (phi));
    0, 0, 0, 0, 0, 0, 0, -1/(Rd * Cd), 0;
    0, 0, 0, 0, 0, 1, -1, 0, 0
    ];

% x = [Vve, Pae, Qa, Pao, Ps, Qi, Qo, Pd, Vc]
B = [(Dm/Rm + Da/Ra)*E*V0;
    -(Dm*E*V0)/(Rm*Cae);
    0;
    -(Da*E*V0)/(Ra*Cao);
    0;
    ((betai - gamma)/Cp) * Vd - betai * E * V0;
    ((betao - gamma)/Cp) * Vd;
    Pex / (Rd * Cd);
    0
    ];
end

