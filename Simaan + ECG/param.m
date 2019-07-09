%Parameters
%%

global Li Lo Lp Rp Ra Vo Cao Cla Cs Rs Rm Rc Ls Rd Cd Cp Vd_vad

%cardiovascular system model parameters (taken from Simaan2009 adopted from Breitenstein);
Rs = 1.0000; % Rs-system resistance;(0.83-normal,weak; 1.4-severly weak without pump; 0.83-severly weak with pump)(mmHg.sec/mL)
Rm = 0.0050; % Rm-mitral valve open;(mmHg.sec/mL)
Ra = 0.0060; % Ra-aortic valve open;(mmHg.sec/mL)
Rc = 0.0398; % Rc-characteristic resistance;(mmHg.sec/mL)
Ls = 0.0005; % Ls-inertance of blood in aorta;(mmHg.sec^2/mL)

Cs  = 1.33; % Cs-system compliance;(mL/mmHg)
Cao = 0.08; % Aortic Complinace (ml/mmHg)
Cla  = 4.4; %Cr-pulmonary compliance;(mL/mmHg)

Vo = 10; % x-intercept of linear EDPVR
%%
%Thoratec Parameters
Max_vol = 100;  % (ml) Maximum VAD Volume
Vd_vad  = 107; % (ml)
Vd_vad_v(1) = Vd_vad;

V_total = 370; % total blood volume; (280-normal)(mL) + Vvad(1) + Vao(1)

%Human Values for Cannula
Tcannula_length_out = 32; %cm
Tcannula_length_in = 24; %cm
Tcannula_area = (0.9^2)*pi; %cm^2
rho_blood = 1; %gm/ml
Thoratec_length = 10; %cm
Thoratec_area = 5; %cm

%INITIAL VALUES OF INLET AND OUTLET RESISTANCE FOR THORATEC CANNULA
Ri = 0.15;
Ro = 0.05;

%EXPRESSION OF FINDING INDUCTANCE VALUES IN THE MODEL
Li = 1.75*(9/4)*(rho_blood*Tcannula_length_in)/(1359.5*Tcannula_area); %Rideout text
Lo = (rho_blood*Tcannula_length_out)/(1359.5*Tcannula_area); %Rideout text

Lp = (9/4)*(rho_blood*Thoratec_length)/(1359.5*Thoratec_area);
% Lp = 0; % A ausência desse indutor não altera o comportamento do PVAD

Rp = 0.05; % Resistor in the VAD (mmHg*s/ml)
% Rp = 0; % A ausência desse resistor não altera o comportamento do PVAD
Cp = 2.0;  % ml/mmHg as estimated from PV experiments

Rd = 0.01; % Resistance of Drive line air;
Cd = 4.0; % Complinace of Drive line air;
