%%
clear all
close all
clc

%Simulation Time;
start_t = 0;
Ts      = 0.0001;
end_t   = 8; % end of simulation time
T = start_t:Ts:end_t;
N = length(T);

HR = 60;      % Heart rate;
RR = 60/HR;    % Duration between R waves
tmax = 0.2 + 0.1555*RR;

ts_VAD = 0.50; % Duration time of VAD systole
t_eject = (RR*ts_VAD)/Ts; % Ejection time for VAD in fill-to-empty operation
t_eject_c = t_eject;


xecg = zeros(1,N);
yecg = zeros(1,N);
zecg = zeros(1,N);
maxi        = zeros(1,N);
mini        = zeros(1,N);
h = zeros(1,N);
a           = zeros(1,N);
n           = zeros(1,N);
g           = zeros(1,N);
r           = zeros(1,N);
v           = zeros(1,length(T));
w           = zeros(1,length(T));
R_dtct      = zeros(1,length(T));
z0          = zeros(1,length(T));

xecg(1) = -1;

deltaS = zeros(1,N);
deltaPVAD = zeros(1,N);
gammad = zeros(1,N);

deltaS(1) = 1;
%Initialize the x-matrix (state variables)
% x = [    x(1)     x(2)     x(3) ]';
  x = [ xecg(1)  yecg(1)  zecg(1) ]';

ip = 0;
lg = 0;
ii = 1;
e  = zeros(1,N);

% R detector variables
sigma = 2;
delta = 2;
beta = 15;
atv = 0;

for i = 1:N-1
[xdot, x] = runkut4_mamemi(Ts,x,RR,z0(i));

z0(i+1) = 0.15*sin(2*pi*(60/(12+randn))*T(i+1));
HR = 60 + 2*randn;
RR = 60/HR;

xecg(i+1) = x(1);
yecg(i+1) = x(2);
zecg(i+1) = x(3);

% R-wave detection
if zecg(i+1) > maxi(i) 
    maxi(i+1) = maxi(i) + sigma*delta;    
elseif zecg(i+1) <= maxi(i)
    maxi(i+1) = maxi(i) - delta;
end

if zecg(i+1) < mini(i) 
    mini(i+1) = mini(i) - sigma*delta;    
elseif zecg(i+1) >= mini(i)
    mini(i+1) = mini(i) + delta;
end

h(i+1) = zecg(i+1) - (maxi(i+1)+mini(i+1))/2;

a(i+1) = maxi(i+1)-mini(i+1);

if a(i+1) <= h(i+1)
    n(i+1) = sign(h(i+1)*(abs(h(i+1))-a(i+1)));
else
    n(i+1) = 0;
end

if i > beta
    if (n(i)>0) && (n(i)>n(i-beta)) && (n(i)>n(i+beta))   
        g(i) = n(i) - max(n(i-beta),n(i+beta));
    elseif (n(i)<0) && (n(i)<n(i-beta)) && (n(i)<n(i+beta))   
        g(i) = n(i) + min(n(i-beta),n(i+beta));
    else
        g(i) = 0;
    end

    if(g(i)>g(i-1) && g(i)>g(i+1))
        r(i) = g(i);
    else
        r(i) = 0;
    end

    if(g(i)<g(i-1) && g(i)<g(i+1))
        v(i) = g(i);
    else
        v(i) = 0;
    end
end
    if r(i)>0
        w(i) = r(i);
    end
    
    if v(i)<0
        w(i) = -v(i);
    end

if(w(i) == 1)
    atv = 1;
end

if ((atv ==1) && (zecg(i)>=0.03))
    atv = 2;
end

R_dtct(i) = 0;
if ((atv ==2) && (zecg(i)<=0.03))
    j = 1;
    aux_E = 1;
    atv = 0;
    R_dtct(i) = 1;
end

if (i > ip)
    fprintf('Executing ... \t%d %%\r',lg);
    lg = lg + 10;
    ip = ip + (N-1)/10;
end
end
%%

% Definição do tempo
Ts = 0.0001;
t = 0:Ts:8;

r_wave_times = find(max(zecg)*R_dtct) * Ts;

% Função Elastança
HR = 60;
Emax = .8;
Emin = 0.05;
T = 60/HR;
Tmax = 0.2 + 0.15*T;

ejection_delay = 0.25 * T;
ejection_time = T / 2;

tn = zeros(1, length(t));
En = zeros(1, length(t));
E = zeros(1, length(t));

for i=1:length(t)-1
    [~, index] = min(abs(t(i) - r_wave_times));
    t_c = r_wave_times(index);
    tn(i) = (mod(t(i), t_c))/Tmax;
    
    En(i) = 1.55 * (((tn(i)/0.7)^1.9)/(1+(tn(i)/0.7)^1.9)) * (1/(1+(tn(i)/1.17)^21.9));   
    E(i) = (t(i) > r_wave_times(1))*(Emax - Emin) * En(i) + Emin;
end

% figure(1)
% plot(t, E, 'b', t, 40*max(zecg)*R_dtct, 'r')
%% Equação do sistema

omega = 12000 + 100.*t;

Vve = zeros(1, length(t));
Pao = zeros(1, length(t));
Qa = zeros(1, length(t));
Ps = zeros(1, length(t));
Pae = zeros(1, length(t));
Pve = zeros(1, length(t));
Qo = zeros(1, length(t));
Qi = zeros(1, length(t));
Vc = zeros(1, length(t));
Pd = zeros(1, length(t));
Pex = zeros(1, length(t));
Pc = zeros(1, length(t));
omegai = zeros(1, length(t));
omegao = zeros(1, length(t));
Pe = 180;
Pf = 0;

% x = [Vve, Pae, Qa, Pao, Ps, Qi, Qo, Pd, Vc]
x = [140 5 0 90 90 0 0 0 0]'; % Valores iniciais
V0 = 15;
alpha = 0.15;

ejection_counter = 0;
dQi = 0;
dQo = 0;
Rp = 0.05;
Lp = 0.0033;
for i=1:length(t)-1
    Pve(i) = E(i)*(Vve(i)-V0); % Pressão no ventrículo esquerdo
    
    if ejection_counter > 0
        Pex(i) = Pe;
        ejection_counter = ejection_counter - 1;
    elseif ismember(i, (r_wave_times / Ts) + (ejection_delay / Ts))
        Pex(i) = Pe;
        ejection_counter = ejection_time / Ts;
    else 
        Pex(i) = Pf;
    end
    
    Pair = Pd(i) + alpha * (Qi(i) - Qo(i));
    Pc(i) = Pair + Rp * (Qi(i) - Qo(i)) + Lp * (dQi - dQo);
    
    [A, B] = changeDiodes(Pao(i), Pae(i), Pve(i), E(i), Pc(i), alpha, Pex(i), Vve(i), Qo(i)); % Função diodo
        
    %w = (omega(i)*2*pi/60)^2;
    
    [x, qi, qo] = runkut4(Ts, x, A, B); % Runge-Kuta 4 ordem
    
    dQi = qi;
    dQo = qo;
    % x = [Vve, Pae, Qa, Pao, Ps, Qi, Qo, Pd, Vc]
    Vve(i+1) = x(1);
    Pae(i+1) = x(2);
    Qa(i+1) = x(3);
    Pao(i+1) = x(4);
    Ps(i+1) = x(5);
    Qi(i+1) = x(6);
    Qo(i+1) = x(7);
    Pd(i+1) = x(8);
    Vc(i+1) = x(9);
end

%% Plota tudo

figure(1)

subplot(4, 1, 1);
plot(t, Pao, 'r', t, Pae, 'g', t, Pve, 'b')
legend('Pao','Pae', 'Pve')
grid on
title('Simula��oo do modelo 0D do sistema cardiovascular humano')
ylabel('Press�o (mmHg)')
xlabel('tempo (s)')

subplot(4, 1, 2);
plot(t, Qa, 'm')
grid on
ylabel('Fluxo (m/s)')
xlabel('tempo (s)')
title('Fluxo A�rtico')

subplot(4, 1, 3);
plot(t, Vve, 'r')
grid on
ylabel('Volume (ml)')
xlabel('tempo (s)')
title('Volume no ventr�culo esquerdo')

subplot(4, 1, 4);
plot(t,zecg,'b',t,max(zecg)*R_dtct,'k')
legend('QRS', 'Onda R detectada')
grid on
ylabel('ECG sint�tico')
xlabel('tempo(s)')
title('Sinal de ECG sint�tico com detec��o de onde R (complexo QRS)')