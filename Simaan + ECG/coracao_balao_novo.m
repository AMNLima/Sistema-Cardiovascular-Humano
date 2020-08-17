clear all
close all
clc

%Simulation Time;
start_t = 0;
Ts      = 0.0001;
end_t   = 4; % end of simulation time
T = start_t:Ts:end_t;
N = length(T);

HR = 98;      % Heart rate;
RR = 60/HR;    % Duration between R waves

xecg = zeros(1,N);
yecg = zeros(1,N);
zecg = zeros(1,N);
z0 = zeros(1,length(T));

xecg_P = zeros(1,N);
yecg_P = zeros(1,N);
zecg_P = zeros(1,N);
z0_P = zeros(1, length(T));

xecg_QRS = zeros(1,N);
yecg_QRS = zeros(1,N);
zecg_QRS = zeros(1,N);
z0_QRS = zeros(1,length(T));

xecg_T = zeros(1,N);
yecg_T = zeros(1,N);
zecg_T = zeros(1,N);
z0_T = zeros(1,length(T));

xecg(1) = -1;
xecg_P(1) = -1;
xecg_QRS(1) = -1;
xecg_T(1) = -1;

%Initialize the x-matrix (state variables)
% x = [    x(1)     x(2)     x(3) ]';
  x = [ xecg(1)  yecg(1)  zecg(1) ]';
  x_P = [ xecg_P(1) yecg_P(1) zecg_P(1)]';
  x_QRS = [ xecg_QRS(1) yecg_QRS(1) zecg_QRS(1)]';
  x_T = [ xecg_T(1) yecg_T(1) zecg_T(1)]';
  
for i = 1:N-1
    [xdot, x] = runkut4_PQRST(Ts,x,RR,z0(i), 'ALL');
    
    [xdot_P, x_P] = runkut4_PQRST(Ts,x_P,RR,z0_P(i), 'P');
    
    [xdot_P, x_QRS] = runkut4_PQRST(Ts,x_QRS,RR,z0_QRS(i), 'QRS');
    
    [xdot_T, x_T] = runkut4_PQRST(Ts,x_T,RR,z0_T(i), 'T');
    
    z0(i+1) = 0.15*sin(2*pi*(60/(12+randn))*T(i+1));
    z0_P(i+1) = z0(i+1);
    z0_QRS(i+1) = z0(1+1);
    z0_T(i+1) = z0(i+1);
    
    HR = 60 + 2*randn;
    RR = 60/HR;
    
    xecg(i+1) = x(1);
    yecg(i+1) = x(2);
    zecg(i+1) = x(3);
    
    xecg_P(i+1) = x_P(1);
    yecg_P(i+1) = x_P(2);
    zecg_P(i+1) = x_P(3);
    
    xecg_QRS(i+1) = x_QRS(1);
    yecg_QRS(i+1) = x_QRS(2);
    zecg_QRS(i+1) = x_QRS(3);
    
    xecg_T(i+1) = x_T(1);
    yecg_T(i+1) = x_T(2);
    zecg_T(i+1) = x_T(3);
end

% Definição do tempo
Ts = 0.0001;
t = 0:Ts:4;


%% Plota tudo

figure(1)

plot(t, zecg, 'b', t, zecg_P, 'r',t, zecg_QRS, 'g', t, zecg_T, 'k')
legend('ALL', 'P', 'QRS', 'T')
grid on
ylabel('ECG sintético')
xlabel('tempo(s)')

%%
figure(2)
subplot(2,1,1)
plot(t, zecg_P+zecg_QRS+zecg_T, 'r', t, zecg, 'b')
subplot(2,1,2)
plot(t, ((zecg_P+zecg_QRS+zecg_T)-zecg).^2)


