% Parametros
clc
clear

% Definição do tempo
Ts = 0.0001;
t = 0:Ts:5;

% Função Elastança
HR = 75;
Emax = 2;
Emin = 0.06;
T = 60/HR;
Tmax = 0.2 + 0.15*T;

tn = (mod(t,T))/Tmax;
En = 1.55 * (((tn./0.7).^1.9)./(1+(tn./0.7).^1.9)) .* (1./(1+(tn./1.17).^21.9));
E = (Emax - Emin).*En + Emin;

%plot(E)
%% Equação do sistema

Vve = zeros(1, length(t));
Pao = zeros(1, length(t));
Qa = zeros(1, length(t));
Ps = zeros(1, length(t));
Pae = zeros(1, length(t));
Pve = zeros(1, length(t));

%x = [Vve Pao Qa Ps Pae];
x = [140 90 0 90 5]'; % Valores iniciais
V0 = 10;

for i=1:length(t)-1
    Pve(i) = E(i)*(Vve(i)-V0); % Pressão no ventrículo esquerdo
    
    [A, B] = changeDiodes(Pao(i), Pae(i), Pve(i), E(i)); % Função diodo
    
    x = runkut4(Ts, x, A, B, E(i)); % Runge-Kuta 4 ordem
    
    Vve(i+1) = x(1);
    Pao(i+1) = x(2);
    Qa(i+1) = x(3);
    Ps(i+1) = x(4);
    Pae(i+1) = x(5);
end

%% Plota tudo

figure(1)

subplot(3, 1, 1);
plot(t, Pao, 'r', t, Pae, 'g', t, Pve, 'b')
legend('Pao','Pae', 'Pve')
grid on
title('Simulação do modelo 0D do sistema cardiovascular humano')
ylabel('Pressão (mmHg)')
xlabel('tempo (s)')

subplot(3, 1, 2);
plot(t, Qa, 'm')
grid on
ylabel('Fluxo (m/s)')
xlabel('tempo (s)')
title('Fluxo Aórtico')

subplot(3, 1, 3);
plot(t, Vve, 'r')
grid on
ylabel('Volume (ml)')
xlabel('tempo (s)')
title('Volume no ventrículo esquerdo')