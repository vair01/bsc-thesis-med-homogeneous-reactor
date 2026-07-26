clear; close all; clc
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultaxesfontsize', 16)

% Vair Piova Andrea, June 2024
% Design of a homogeneous reactor for the production of radionuclides for 
% medical purposes, fuel cycle.

%% Data
% Reactor data
pow=50e3;   % (W)
diameter=1e2;   % (cm)
height=diameter;    % (cm)
T_core=60+273.15;   % (K)
T_0=20+273.15;  % (K)

% Geometric calculations
vol=pi*(diameter^2)/4*height;    % (cm^3)
dens_pot=pow*1000/vol;      % (W/cm^3)

% Constants
N_av=6.022e23;    % Avogadro's number (#atoms/mol)
MM_water=2*1.00794+15.9994;    % (g/mol) => Appendix A, Murray
MM_abs=2*157.25+3*15.9994;   % (g/mol) => Appendix A, Murray

% Density of water @ 60°C
rho_m=0.983;   % (g/cm^3)

% Molybdenum properties
rho_Mo=10.280;  % (g/cm^3)
MM_Mo98=98; % (g/mol)
MM_Mo99=99; % (g/mol)
half_life_Mo99=65.976*3600; % (s)
lambda_Mo99=log(2)/half_life_Mo99;  % (1/s)
yield_235=6.132/100;    % Mo-99 production rate from U-235 fission
yield_239=6.185/100;    % Mo-99 production rate from Pu-239 fission

% Mo-98 pellet => cylinder with d=2cm, h=2cm
d_pellet=2;   % (cm)
h_pellet=2;   % (cm)
vol_pellet=pi*d_pellet^2/4*h_pellet;  % (cm^3)
N_Mo98=rho_Mo*N_av/MM_Mo98; % (#atoms/cm^3)

% Calculation of the average microscopic cross-sections for transmutation:
% Sigma @ E0
sigE0_f_235=586.7371;  % (barn)
sigE0_abs_235=sigE0_f_235+99.39421;  % (barn)
sigE0_f_239=747.353;    % (barn)
sigE0_c_238=2.683199;   % (barn)
sigE0_c_Mo98=0.130015;  % (barn)
sigE0_c_Mo99=8.00243;   % (barn)

% Westcott's corrective factors, from Janis
g_U235=0.9755;
g_U238=1.0049;
g_Pu239=1.05203;
g_Mo98=1.00087;
g_Mo99=1.00038;

sigma_f_235=sqrt(pi)/2*(T_0/T_core)^(1/2)*sigE0_f_235*g_U235; % (barn)
sigma_f_239=sqrt(pi)/2*(T_0/T_core)^(1/2)*sigE0_f_239*g_Pu239;    % (barn)
sigma_c_238=sqrt(pi)/2*(T_0/T_core)^(1/2)*sigE0_c_238*g_U238; % (barn)
sigma_abs_235=sqrt(pi)/2*(T_0/T_core)^(1/2)*sigE0_abs_235*g_U235; % (barn)
sigma_c_Mo98=sqrt(pi)/2*(T_0/T_core)^(1/2)*sigE0_c_Mo98*g_Mo98; % (barn)
sigma_c_Mo99=sqrt(pi)/2*(T_0/T_core)^(1/2)*sigE0_c_Mo99*g_Mo99; % (barn)

%% Time discretisation
months=6;
t_end=months*30*24*3600;   % (s)
dt=3600;    % one step every hour
tt=(0:dt:t_end);
nt=length(tt);

% Flux functions dependent on U235 and Pu239
flusso=@(N235,N239) (pow/(3.2044e-11*(sigma_f_235*1e-24*N235+ ...
    sigma_f_239*1e-24*N239)*vol)); %(n/cm^2/s)
% Flux at reactor center
rr=0;
hh=0;
loc_flux=@(N235,N239) 3.63*pow/(3.2044e-11*( ...
    sigma_f_235*1e-24*N235+sigma_f_239*1e-24*N239)*vol)*besselj( ...
    0,2.405*rr/(diameter/2))*cos(pi*hh/height);

%% FUEL CYCLE without removal (=> Forward Euler method)
% Initial conditions
% Enrichment, initial moderation ratio and Nf
mod_ratio=1.2253;  % Nm/Nf (-)
enr=3/100;  % (-)

Nm=rho_m*N_av/MM_water; % (#molecules/cm^3)
Nf=Nm/mod_ratio; % (#molecules/cm^3)

% Initial concentrations
N235_old=Nf*enr;    % (#atoms/cm^3)
N238_old=Nf*(1-enr);    % (#atoms/cm^3)
N239_old=0; % (#atoms/cm^3)
N239_new=N239_old;  % (#atoms/cm^3)
NMo98_old=N_Mo98;   % (#atoms/cm^3)
NMo99_old=0;    % (#atoms/cm^3)
NMo99fuel_old=0;    % (#atoms/cm^3)

% Initialising bisection
maxiter=1e4;
inf=1e15;
sup=1e19;
toll=1e-10;

fun_keff=@(Nabs) criticality(Nabs,Nf,enr,N239_new,height,diameter,T_core);
[Nabs_new]=bisection(fun_keff,inf,sup,maxiter,toll);

% Evaluation of the average flux and flux at the center of the reactor
flux=flusso(N235_old,N239_old); % (n/cm^2/s)
flux_center=loc_flux(N235_old,N239_old);  % (n/cm^2/s)

% Vector pre-allocation
k_eff=zeros(size(tt));
U235=zeros(size(tt));
U238=zeros(size(tt));
Pu239=zeros(size(tt));
N_Gdox=zeros(size(tt));
Mo99=zeros(size(tt));
Mo99_fuel=zeros(size(tt));
flux_vec=zeros(size(tt));
flux_vec_center=zeros(size(tt));

% Variables @ t=0
U235(1)=N235_old;   % (#atoms/cm^3)
U238(1)=N238_old;   % (#atoms/cm^3)
Pu239(1)=N239_old;  % (#atoms/cm^3)
N_Gdox(1)=Nabs_new; % (#atoms/cm^3)
Mo99(1)=NMo99_old;  % (#atoms/cm^3)
Mo99_fuel(1)=NMo99fuel_old;     % (#atoms/cm^3)
flux_vec(1)=flux; % (n/cm^2/s)
flux_vec_center(1)=flux_center;   % (n/cm^2/s)

% Check that k_eff is equal to 1 (=> criticality)
k_eff(1)= 1 - criticality(Nabs_new,Nf,enr,N239_old,height,diameter,T_core); % (-)

for ii=2:nt

    % Transmutation eq.
    N235_new=N235_old*(1-dt*sigma_abs_235*1e-24*flux);  % (#atoms/cm^3)
    N238_new=N238_old*(1-dt*sigma_c_238*1e-24*flux);    % (#atoms/cm^3)
    N239_new=N239_old*(1-dt*sigma_f_239*1e-24*flux)+ ...
        sigma_c_238*N238_old*1e-24*flux*dt; % (#atoms/cm^3)
    % Production of Mo-99 by transmutation of Mo-98 in the pellet
    NMo98_new=NMo98_old*(1-dt*sigma_c_Mo98*1e-24*flux_center);  % (#atoms/cm^3)
    NMo99_new=NMo99_old*(1-dt*lambda_Mo99-dt*sigma_c_Mo99*1e-24*flux)+ ...
        dt*sigma_c_Mo98*1e-24*NMo98_old*flux_center; % (#atoms/cm^3)   
    % Production of Mo-99 directly from U-235 fission
    NMo99fuel_new=NMo99fuel_old*(1-dt*lambda_Mo99-dt*sigma_c_Mo99*1e-24*flux)+ ...
        dt*sigma_f_235*1e-24*N235_old*flux*yield_235+ ...
        dt*sigma_f_239*1e-24*N239_old*flux*yield_239; % (#atoms/cm^3)
    
    Nf=N235_new+N238_new;   % (#atoms/cm^3)
    enr=N235_new/Nf;    % (-)
    mod_ratio=Nm/Nf;
    MM_uranyl=235*enr+238*(1-enr)+8*15.9994+2*14.0067;    % (g/mol) => UO2(NO3)2

    % The bisection method is used to calculate Nabs and obtain k_eff = 1
    fun_keff=@(Nabs) criticality(Nabs,Nf,enr,N239_new,height,diameter,T_core);
    [Nabs_new]=bisection(fun_keff,inf,sup,maxiter,toll);
   
    % Check that k_eff is equal to 1
    k_eff(ii)= 1 - criticality(Nabs_new,Nf,enr,N239_old,height,diameter,T_core);    % (-)

    % Flux evaluation
    flux=flusso(N235_new,N239_new); %(n/cm^2/s)
    flux_center=loc_flux(N235_new,N239_new);  %(n/cm^2/s)
    
    % Concentrations update
    N235_old=N235_new;  % (#atoms/cm^3)
    N238_old=N238_new;  % (#atoms/cm^3)
    N239_old=N239_new;  % (#atoms/cm^3)
    NMo98_old=NMo98_new;    % (#atoms/cm^3)
    NMo99_old=NMo99_new;    % (#atoms/cm^3) 
    NMo99fuel_old=NMo99fuel_new;    % (#atoms/cm^3)

    % Store results in respective vectors
    U235(ii)=N235_new;  % (#atoms/cm^3)
    U238(ii)=N238_new;  % (#atoms/cm^3)
    Pu239(ii)=N239_new; % (#atoms/cm^3)
    N_Gdox(ii)=Nabs_new;    % (#atoms/cm^3)
    Mo99(ii)=NMo99_new; % (#atoms/cm^3)
    Mo99_fuel(ii)=NMo99fuel_new;    % (#atoms/cm^3)
    flux_vec(ii)=flux;    % (n/cm^2/s)
    flux_vec_center(ii)=flux_center;  % (n/cm^2/s)

end

%% Plots (First Cycle)
figure(1)
subplot(2,2,1)
plot(tt/(24*3600),U235)
ylabel('N_{U235} (atoms/cm^3)')
xlabel('Time (days)')
xlim([0 180])
grid on
box on
subplot(2,2,2)
plot(tt/(24*3600),U238)
ylabel('N_{U238} (atoms/cm^3)')
xlabel('Time (days)')
xlim([0 180])
grid on
box on
subplot(2,2,3)
plot(tt/(24*3600),Pu239)
ylabel('N_{Pu239} (atoms/cm^3)')
xlabel('Time (days)')
xlim([0 180])
grid on
box on
subplot(2,2,4)
plot(tt/(3600*24),round(k_eff,5))
ylabel('k_{eff} (-)')
xlabel('Time (days)')
xlim([0 180])
grid on
box on
sgtitle('Temporal evolution of isotopic concentrations and k_{eff}', 'FontWeight', 'bold')

figure(2)
plot(tt/(3600*24),N_Gdox)
title('Change of Gd_2O_3 concentration in time')
ylabel('N_{Gd_2O_3} (molecules/cm^3)')
xlabel('Time (days)')
xlim([0 180])
grid on
box on

figure(3)
plot(tt/(3600*24),Mo99)
title('Activity of Mo-99 (transmutation of Mo-98 in the pellet)')
ylabel('N_{Mo-99} (atoms/cm^3)')
xlabel('Time (days)')
xlim([0 180])
grid on
box on

figure(4)
plot(tt/(3600*24),Mo99_fuel)
title('Activity of Mo-99 (U-235 fission)')
ylabel('N_{Mo-99} (atoms/cm^3)')
xlabel('Time (days)')
xlim([0 180])
grid on
box on

disp('Press any button to continue to the cycle with Mo-99 removal after 7 days')
pause; close all; clc

%% FUEL CYCLE with Mo-99 removal after 7 days (=> Forward Euler method)
% Initial conditions
% Initial enrichment and Nf
enr=3/100;  % (-)

Nm=rho_m*N_av/MM_water; % (#molecules/cm^3)
Nf=Nm/mod_ratio; % (#molecules/cm^3)

% Cycle counter
n=1;

% Initial concentrations
N235_old=Nf*enr;    % (#atoms/cm^3)
N238_old=Nf*(1-enr);    % (#atoms/cm^3)
N239_old=0; % (#atoms/cm^3)
N239_new=N239_old;  % (#atoms/cm^3)
NMo98_old=N_Mo98;   % (#atoms/cm^3)
NMo99_old=0;    % (#atoms/cm^3)
NMo99fuel_old=0;    % (#atoms/cm^3)

fun_keff=@(Nabs) criticality(Nabs,Nf,enr,N239_new,height,diameter,T_core);
[Nabs_new]=bisection(fun_keff,inf,sup,maxiter,toll);

flux=flusso(N235_old,N239_old); % (n/cm^2/s)
flux_center=loc_flux(N235_old,N239_old);  % (n/cm^2/s)

% Vector pre-allocation
k_eff=zeros(size(tt));
U235=zeros(size(tt));
U238=zeros(size(tt));
Pu239=zeros(size(tt));
N_Gdox=zeros(size(tt));
Mo99=zeros(size(tt));
Mo99_fuel=zeros(size(tt));
flux_vec=zeros(size(tt));
flux_vec_center=zeros(size(tt));

% Variables @ t=0
U235(1)=N235_old;   % (#atoms/cm^3)
U238(1)=N238_old;   % (#atoms/cm^3)
Pu239(1)=N239_old;  % (#atoms/cm^3)
N_Gdox(1)=Nabs_new; % (#atoms/cm^3)
Mo99(1)=NMo99_old;  % (#atoms/cm^3)
Mo99_fuel(1)=NMo99fuel_old;     % (#atoms/cm^3)
flux_vec(1)=flux; % (n/cm^2/s)
flux_vec_center(1)=flux_center;   % (n/cm^2/s)

% Correctly initialize k_eff(1) at time zero for the second cycle
k_eff(1)= 1 - criticality(Nabs_new,Nf,enr,N239_old,height,diameter,T_core); % (-)

for ii=2:nt

    % Transmutation eq.
    N235_new=N235_old*(1-dt*sigma_abs_235*1e-24*flux);  % (#atoms/cm^3)
    N238_new=N238_old*(1-dt*sigma_c_238*1e-24*flux);    % (#atoms/cm^3)
    N239_new=N239_old*(1-dt*sigma_f_239*1e-24*flux)+sigma_c_238*N238_old*1e-24*flux*dt; % (#atoms/cm^3)
    % Production of Mo-99 by transmutation of Mo-98 in the pellet
    NMo98_new=NMo98_old*(1-dt*sigma_c_Mo98*1e-24*flux_center);  % (#atoms/cm^3)
    NMo99_new=NMo99_old*(1-dt*lambda_Mo99-dt*sigma_c_Mo99*1e-24*flux)+ ...
        dt*sigma_c_Mo98*1e-24*NMo98_old*flux_center; % (#atoms/cm^3)   
    % Production of Mo-99 directly from U-235 fission
    NMo99fuel_new=NMo99fuel_old*(1-dt*lambda_Mo99-dt*sigma_c_Mo99*1e-24*flux)+ ...
        dt*sigma_f_235*1e-24*N235_old*flux*yield_235+ ...
        dt*sigma_f_239*1e-24*N239_old*flux*yield_239; % (#atoms/cm^3)
    
    if ii==168*n    % every 7 days
        activity = NMo99fuel_old*vol*lambda_Mo99/3.7e10;
        NMo98_new=N_Mo98;   % (#atoms/cm^3)
        NMo99_new=0;    % (#atoms/cm^3)
        NMo99fuel_new=0;
        n=n+1;
    end

    Nf=N235_new+N238_new;   % (#atoms/cm^3)
    enr=N235_new/Nf;    % (-)
    MM_uranyl=235*enr+238*(1-enr)+8*15.9994+2*14.0067;    % (g/mol) => UO2(NO3)2

    % The bisection method is used to calculate Nabs and obtain k_eff = 1
    fun_keff=@(Nabs) criticality(Nabs,Nf,enr,N239_new,height,diameter,T_core);
    [Nabs_new]=bisection(fun_keff,inf,sup,maxiter,toll);
   
    % Check that k_eff is equal to 1
    k_eff(ii)= 1 - criticality(Nabs_new,Nf,enr,N239_old,height,diameter,T_core);    % (-)

    % Flux evaluation
    flux=flusso(N235_new,N239_new); %(n/cm^2/s)
    flux_center=loc_flux(N235_new,N239_new);  %(n/cm^2/s)
    
    % Concentrations update
    N235_old=N235_new;  % (#atoms/cm^3)
    N238_old=N238_new;  % (#atoms/cm^3)
    N239_old=N239_new;  % (#atoms/cm^3)
    NMo98_old=NMo98_new;    % (#atoms/cm^3)
    NMo99_old=NMo99_new;    % (#atoms/cm^3) 
    NMo99fuel_old=NMo99fuel_new;    % (#atoms/cm^3)

    % Store results in respective vectors
    U235(ii)=N235_new;  % (#atoms/cm^3)
    U238(ii)=N238_new;  % (#atoms/cm^3)
    Pu239(ii)=N239_new; % (#atoms/cm^3)
    N_Gdox(ii)=Nabs_new;    % (#atoms/cm^3)
    Mo99(ii)=NMo99_new; % (#atoms/cm^3)
    Mo99_fuel(ii)=NMo99fuel_new;    % (#atoms/cm^3)
    flux_vec(ii)=flux;    % (n/cm^2/s)
    flux_vec_center(ii)=flux_center;  % (n/cm^2/s)

end

%% Plots (Second Cycle)
figure(1)
subplot(2,2,1)
plot(tt/(24*3600),U235)
ylabel('N_{U235} (atoms/cm^3)')
xlabel('Time (days)')
xlim([0 180])
grid on
box on
subplot(2,2,2)
plot(tt/(24*3600),U238)
ylabel('N_{U238} (atoms/cm^3)')
xlabel('Time (days)')
xlim([0 180])
grid on
box on
subplot(2,2,3)
plot(tt/(24*3600),Pu239)
ylabel('N_{Pu239} (atoms/cm^3)')
xlabel('Time (days)')
xlim([0 180])
grid on
box on
subplot(2,2,4)
plot(tt/(3600*24),round(k_eff,5))
ylabel('k_{eff}')
xlabel('Time (days)')
xlim([0 180])
grid on
box on
sgtitle('Temporal evolution of isotopic concentrations and k_{eff}', 'FontWeight', 'bold')

figure(2)
plot(tt/(3600*24),N_Gdox)
title('Change of Gd_2O_3 concentration in time')
ylabel('N_{Gd_2O_3} (molecules/cm^3)')
xlabel('Time (days)')
xlim([0 180])
grid on
box on

figure(3)
plot(tt/(3600*24),Mo99*vol_pellet*lambda_Mo99/3.7e10)
title('Activity of Mo-99 (transmutation of Mo-98 in the pellet)')
ylabel('A_{Mo-99} (Ci)')
xlabel('Time (days)')
xlim([0 180])
grid on
box on
figure(4)
plot(tt/(3600*24),Mo99_fuel*vol*lambda_Mo99/3.7e10)
title('Activity of Mo-99 (U-235 fission)')
ylabel('A_{Mo-99} (Ci)')
xlabel('Time (days)')
xlim([0 180])
grid on
box on

%% Evaluation of 6-day activity and fraction of global production
tempo=6*24*3600;    % (s)
activity_6day=activity*exp(-lambda_Mo99*tempo)

perc_global_prod=activity_6day/9000*100

%% Function that runs the bisection method (non-linear eq.)
function [x,k] = bisection(f,a,b,kmax,tol)

if f(a)*f(b) > 0
    error('Change the bounds for the bisection method calculation');
end

fa=f(a);
x=(a+b)/2;
fx=f(x);
fr=fx;

for k = 1:kmax
    if abs(fx)/abs(fr) <= tol
        break
    else
        if fa*fx <0
            b=x;
        else
            a = x;
            fa = fx;
        end
    end
x = (a+b)/2;
fx=f(x);
end
end


%% Function that calculates criticality
function fun_keff=criticality(Nabs,Nf,enr,Pu239,alt,diam,Tcore)

% Constants
N_av=6.022e23;  % (#atoms/mol)
T_0=20+273.15;  % (K)
nu_U235=2.45; % average number of neutrons produced from U235 fission (-)
nu_Pu239=2.88;  % average number of neutrons produced from Pu239 fission (-)
lethargy_mod=0.920; % moderator lethargy, Table 3.1 Lamarsh (-)
MM_abs=2*157.25+3*15.9994;   % (g/mol) => Appendix A, Murray
MM_Pu239=239.052163;   % (g/mol) => Appendix A, Murray

% Water properties (@ 60°C)
rho_m=0.9832;   % (g/cm^3)
MM_water=2*1.00794+15.9994;    % (g/mol) => Appendix A, Murray

% Amount of Gd-157 in nature
enr_Gd=15.68/100;   % (-)

% Westcott's corrective factors, from Janis 
g_f_U235=0.97963;
g_c_U235=0.99246;
g_f_U238=1.00152;
g_c_U238=1.00193;
g_f_Pu239=1.01667;
g_c_Pu239=1.01673;
g_c_O=1.00042;
g_c_N=1.00017;
g_c_H=1.00034;
g_c_Gd157=0.85333;

% Various cross-sections @ E0
% cross section U235
sigE0_f_235=586.7371;  % (barn)
sigE0_c_235=99.39421;   % (barn)
% cross section U238
sigE0_f_238=1.85134e-5; % (barn)
sigE0_c_238=2.683199;   % (barn)
% cross section Pu239
sigE0_f_239=747.353;    % (barn)
sigE0_c_239=270.124;    % (barn)
% cross section N
sigE0_abs_N=0.075014;    % (barn)
% cross section O
sigE0_abs_O=1.69912e-4; % (barn)
% cross section H
sigE0_abs_H=0.3327;    % (barn)
% cross section Gd157
sigE0_abs_Gd157=252912; % (barn)

% Sigma scattering @ 1 eV (epithermal zone)
sigma_s_235=12.6; % (barn)
sigma_s_238=9.1;  % (barn)
sigma_s_O=3.8; % (barn)
sigma_s_N=9.9;   % (barn)
sigma_s_H=20.7;    % (barn)
sigma_s_Gd157=13.5183;  % (barn)
sigma_s_Pu239=10.1592;  % (barn)

% The cross-sections (functions of T) are corrected
% => Hp: the flux follows a Gaussian distribution
% Average sigma for fission
sigma_f_235=sqrt(pi)/2*(T_0/Tcore)^(1/2)*sigE0_f_235*g_f_U235; % (barn)
sigma_f_239=sqrt(pi)/2*(T_0/Tcore)^(1/2)*sigE0_f_239*g_f_Pu239;    % (barn)
sigma_f_238=sqrt(pi)/2*(T_0/Tcore)^(1/2)*sigE0_f_238*g_f_U238;    % (barn)
% Average sigma for capture
sigma_c_235=sqrt(pi)/2*(T_0/Tcore)^(1/2)*sigE0_c_235*g_c_U235; % (barn)
sigma_c_239=sqrt(pi)/2*(T_0/Tcore)^(1/2)*sigE0_c_239*g_c_Pu239; % (barn)
sigma_c_238=sqrt(pi)/2*(T_0/Tcore)^(1/2)*sigE0_c_238*g_c_U238; % (barn)
% Average sigma for absorption (fission + capture)
sigma_abs_235=sigma_c_235+sigma_f_235; % (barn)
sigma_abs_239=sigma_c_239+sigma_f_239; % (barn)
sigma_abs_238=sigma_c_238+sigma_f_238; % (barn)
sigma_abs_N=sqrt(pi)/2*(T_0/Tcore)^(1/2)*sigE0_abs_N*g_c_N;    % (barn)
sigma_abs_O=sqrt(pi)/2*(T_0/Tcore)^(1/2)*sigE0_abs_O*g_c_O;    % (barn)
sigma_abs_H=sqrt(pi)/2*(T_0/Tcore)^(1/2)*sigE0_abs_H*g_c_H;    % (barn)
sigma_abs_Gd157=sqrt(pi)/2*(T_0/Tcore)^(1/2)*sigE0_abs_Gd157*g_c_Gd157;    % (barn)

sigma_abs_mod=sigma_abs_H*2+sigma_abs_O;    % (barn)
sigma_abs_abs=sigma_abs_Gd157*2*enr_Gd+sigma_abs_O*3;   % (barn)

MM_uranyl=235*enr+238*(1-enr)+8*15.9994+2*14.0067;    % (g/mol) => UO2(NO3)2

% Criticality condition
keff=1; % (-) => must remain constant over time
Nm=rho_m*N_av/MM_water; % (#molecules/cm^3)
mod_ratio=Nm/Nf; % (-)

% 4 factors evaluation
% eta => reproduction factor
sigma_abs_fuel=sigma_abs_235*enr+sigma_abs_238*(1-enr)+sigma_abs_O*8+sigma_abs_N*2;
eta=(nu_U235*sigma_f_235*Nf*enr+nu_Pu239*sigma_f_239*Pu239)/ ...
(sigma_abs_fuel*Nf+sigma_abs_239*Pu239);

% f => thermal utilization factor
f=@(Nabs)(sigma_abs_fuel*Nf+sigma_abs_239*Pu239)/(sigma_abs_fuel*Nf+ ...
    sigma_abs_mod*Nm+sigma_abs_abs*Nabs+sigma_abs_239*Pu239);

% p => resonance escape probability
Na=Nf*(1-enr);
sigma_p= @(Nabs) (sigma_s_235*Nf*enr+sigma_s_238*Nf*(1-enr)+ ...
    8*sigma_s_O*Nf+2*sigma_s_N*Nf+sigma_s_Pu239*Pu239+2*sigma_s_H*Nm+ ...
    sigma_s_O*Nm+2*Nabs*enr_Gd*sigma_s_Gd157+3*Nabs*sigma_s_O);    %cm^-1
I=@(Nabs) 2.73*(sigma_p(Nabs)/Na).^0.486; % (barn) => typical a,c constants for U238
lethargy_fuel=1-((MM_uranyl-1)^2/(2*MM_uranyl))*log((MM_uranyl+1)/(MM_uranyl-1));
lethargy_abs=1-((MM_abs-1)^2/(2*MM_abs))*log((MM_abs+1)/(MM_abs-1));
lethargy_Pu=1-((MM_Pu239-1)^2/(2*MM_Pu239))*log((MM_Pu239+1)/(MM_Pu239-1));
% The sigma values in the epithermal zone are calculated
sigma_s_mod=sigma_s_H*2+sigma_s_O;  % (barn)
sigma_s_fuel=sigma_s_235*enr+sigma_s_238*(1-enr)+8*sigma_s_O+2*sigma_s_N;   % (barn)
sigma_s_abs=2*sigma_s_Gd157*enr_Gd+3*sigma_s_O; % (barn)
lethargy= @(Nabs) (lethargy_mod*sigma_s_mod*mod_ratio+ ...
    lethargy_abs*sigma_s_abs.*Nabs/Nf+lethargy_fuel*sigma_s_fuel+ ...
    lethargy_Pu*sigma_s_Pu239*Pu239/Nf)./(sigma_s_mod*mod_ratio+ ...
    sigma_s_fuel+sigma_s_Pu239*Pu239/Nf+sigma_s_abs.*Nabs/Nf);
p=@(Nabs) exp(-(Na.*I(Nabs))./(lethargy(Nabs).*sigma_p(Nabs)));

% epsilon => fast fission factor
epsilon=(1+0.690*(1-enr)*1/mod_ratio)/(1+0.563*(1-enr)*1/mod_ratio);

% The probabilities of non-leakage (thermal and fast) are calculated using 
% a simplified method
tau_th=26;  % (cm^2) => Murray's Table 4.4
Lm=2.85;   % (cm) => Murray's Table 4.4
j0=2.40483;
B_2=(pi/alt)^2+(j0/(diam/2))^2;
L_2=@(Nabs) Lm^2*(1-f(Nabs));
P_nl=@(Nabs) (1./(1+L_2(Nabs).*B_2))*exp(-B_2*tau_th);

% Final function to be included in the bisection method
fun_keff=(keff-eta*f(Nabs)*p(Nabs)*epsilon*P_nl(Nabs));
end