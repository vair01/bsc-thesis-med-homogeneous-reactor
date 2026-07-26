clear; close all; clc
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultaxesfontsize', 20)

% Vair Piova Andrea, June 2024
% Design of a homogeneous reactor for the production of radionuclides for 
% medical purposes, initial evaluation, absence of neutron poison.

%% Data
% Reactor data
pow=50e3;   % (W)
T_core=60+273.15; % (K)
diam=1e2;   % (cm)
height=diam;    % (cm)

% Geometrical calculations
vol=pi*(diam^2)/4*height;    % (cm^3)
pow_dens=pow/vol;      % (kW/litro)

% Constants
T_0=20+273.15;  % (K)
nu=2.45; % average number of neutrons produced from U235 fission (-)
N_av=6.022e23;    % Avogadro's number (#atoms/mol)

% Water properties (@ 60°C)
MM_water=2*1.00794+15.9994;    % molar mass (g/mol)
rho_m=0.9832;   % density (g/cm^3)
lethargy_mod=0.920; % moderator lethargy, water (-), table 3.1 Lamarsh

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

% various cross sections @ E0 = 0.0253 eV
% cross section U235
sigE0_f_235=586.7371;  % (barn)
sigE0_c_235=99.39421;   % (barn)
sigE0_abs_235=sigE0_f_235+sigE0_c_235;  % (barn)
sigE0_tot_235=698.218;  % (barn)
% cross section U238
sigE0_f_238=1.85134e-5; % (barn)
sigE0_c_238=2.683199;   % (barn)
sigE0_abs_238=sigE0_c_238+sigE0_f_238;   % (barn)
sigE0_tot_238=11.9229;   % (barn)
% cross section Pu239
sigE0_f_239=747.353;    % (barn)
sigE0_c_239=270.124;    % (barn)
sigE0_abs_239=sigE0_f_239+sigE0_c_239;  % (barn)
% cross section N
sigE0_abs_N=0.075014;    % (barn)
sigE0_tot_N=12.171;    % (barn)
% cross section O
sigE0_abs_O=1.69912e-4; % (barn)
sigE0_tot_O=3.9137; % (barn)
% cross section H
sigE0_abs_H=0.3327;    % (barn)
sigE0_tot_H=30.4139;    % (barn)
% cross section Gd157
sigE0_abs_Gd157=252912; % (barn)

% scattering cross section @ 1 eV (in the epithermal zone)
sigma_s_235=12.6; % (barn)
sigma_s_238=9.1;  % (barn)
sigma_s_O=3.8; % (barn)
sigma_s_N=9.9;   % (barn)
sigma_s_H=20.7;    % (barn)

sigma_s_mod=sigma_s_H*2+sigma_s_O;  % (barn)

% Temperature correction for temperature-dependent cross sections
% => Hp: flux follows a Gaussian distribution
% Average sigma for fission
sigma_f_235=sqrt(pi)/2*(T_0/T_core)^(1/2)*sigE0_f_235*g_f_U235; % (barn)
sigma_f_239=sqrt(pi)/2*(T_0/T_core)^(1/2)*sigE0_f_239*g_f_Pu239;    % (barn)
sigma_f_238=sqrt(pi)/2*(T_0/T_core)^(1/2)*sigE0_f_238*g_f_U238;    % (barn)
% Average sigma for capture
sigma_c_235=sqrt(pi)/2*(T_0/T_core)^(1/2)*sigE0_c_235*g_c_U235; % (barn)
sigma_c_239=sqrt(pi)/2*(T_0/T_core)^(1/2)*sigE0_c_239*g_c_Pu239; % (barn)
sigma_c_238=sqrt(pi)/2*(T_0/T_core)^(1/2)*sigE0_c_238*g_c_U238; % (barn)
% Average sigma for absorption (fission + capture)
sigma_abs_235=sigma_c_235+sigma_f_235; % (barn)
sigma_abs_239=sigma_c_239+sigma_f_239; % (barn)
sigma_abs_238=sigma_c_238+sigma_f_238; % (barn)
sigma_abs_N=sqrt(pi)/2*(T_0/T_core)^(1/2)*sigE0_abs_N*g_c_N;    % (barn)
sigma_abs_O=sqrt(pi)/2*(T_0/T_core)^(1/2)*sigE0_abs_O*g_c_O;    % (barn)
sigma_abs_H=sqrt(pi)/2*(T_0/T_core)^(1/2)*sigE0_abs_H*g_c_H;    % (barn)
sigma_abs_Gd157=sqrt(pi)/2*(T_0/T_core)^(1/2)*sigE0_abs_Gd157*g_c_Gd157;    % (barn)

sigma_abs_mod=sigma_abs_H*2+sigma_abs_O;    % (barn)

% Variable moderation ratio and enrichment vectors
mod_ratio_var=linspace(0.005,100,1e5);
enr_var=linspace(1,5,5)/100;

% Vector pre-allocation
k_inf=zeros(size(mod_ratio_var));
k_eff=zeros(size(mod_ratio_var));
eta=zeros(size(mod_ratio_var));
f=zeros(size(mod_ratio_var));
p=zeros(size(mod_ratio_var));
epsilon=zeros(size(mod_ratio_var));
P_nl=zeros(size(mod_ratio_var));

%% Double for loop
for jj=1:length(enr_var)
    enr=enr_var(jj);

    for ii=1:length(mod_ratio_var)
        mod_ratio=mod_ratio_var(ii);

        MM_uranyl=235*enr+238*(1-enr)+8*15.9994+2*14.0067;    % Uranyl nitrate: UO2(NO3)2 (g/mol)
        Nm=rho_m*N_av/MM_water; % (#molecules/cm^3)
        Nf=Nm/mod_ratio; % (#molecules/cm^3)
    
        % 4 factors evaluation
        % eta => reproduction factor (-)      
        sigma_abs_fuel=sigma_abs_235*enr+sigma_abs_238*(1-enr)+sigma_abs_O*8+sigma_abs_N*2;
        eta(ii)=(nu*sigma_f_235*enr)/sigma_abs_fuel;
        
        % f => thermal utilization factor (-)
        f(ii)=sigma_abs_fuel/(sigma_abs_fuel+sigma_abs_mod*mod_ratio);
        
        % p => resonance escape probability (-)
        Na=Nf*(1-enr);  % resonance absorber density, U238 (#atoms/cm^3)
        sigma_p=(sigma_s_235*Nf*enr+sigma_s_238*Nf*(1-enr)+8*sigma_s_O*Nf+2*sigma_s_N*Nf+2*sigma_s_H*Nm+sigma_s_O*Nm);    % (barn*#molecules/cm^3)
        I=2.73*(sigma_p/Na)^0.486; % resonance integral, typical a,c constants for U238 (barn)
        lethargy_fuel=1-((MM_uranyl-1)^2/(2*MM_uranyl))*log((MM_uranyl+1)/(MM_uranyl-1));   % (-)

        sigma_s_fuel=sigma_s_235*enr+sigma_s_238*(1-enr)+8*sigma_s_O+2*sigma_s_N; % (barn)
        lethargy=(lethargy_mod*sigma_s_mod*mod_ratio+lethargy_fuel*sigma_s_fuel)/(sigma_s_mod*mod_ratio+sigma_s_fuel);    % (-)
        p(ii)=exp(-(Na*I)/(lethargy*sigma_p));
        
        % epsilon => fast fission factor (-)
        epsilon(ii)=(1+0.690*(1-enr)*(1/mod_ratio))/(1+0.563*(1-enr)*(1/mod_ratio));
    
        % Get k_inf (-)
        k_inf(ii)=eta(ii)*f(ii)*p(ii)*epsilon(ii);

        % Calculate P_nl to find k_eff
        tau_th=26;  % Fermi neutron age (cm^2), taken from Lamarsh Table 3.1
        Lm=2.85;   % Diffusion length of the moderator only (cm)
        j0=2.40483; % First root of the Bessel function
        B_2=(pi/height)^2+(j0/(diam/2))^2; % Geometric buckling (1/cm^2)
        L_2=Lm^2*(1-f(ii)); % Approximation (cm^2)
        P_nl(ii)=(1/(1+L_2*B_2))*exp(-B_2*tau_th);
        k_eff(ii)=k_inf(ii)*P_nl(ii);
    end
    
    %% Plots
    % Plot of the 4 factors
    figure(1)
    plot(mod_ratio_var,eta,'displayname',['enr=',num2str(enr*100),'%'])
    title('Effect of moderation ratio \Theta and enrichment on \eta')
    xlabel('\Theta (-)')
    ylabel('\eta (-)')
    legend('-dynamiclegend')
    grid on
    hold on

    figure(2)
    plot(mod_ratio_var,f,'displayname',['enr=',num2str(enr*100),'%'])
    title('Effect of moderation ratio \Theta and enrichment on f')
    xlabel('\Theta (-)')
    ylabel('f (-)')
    legend('-dynamiclegend')
    ylim([0 1])
    grid on
    hold on

    figure(3)
    plot(mod_ratio_var,p,'displayname',['enr=',num2str(enr*100),'%'])
    title('Effect of moderation ratio \Theta and enrichment on p')
    xlabel('\Theta (-)')
    ylabel('p (-)')
    legend('-dynamiclegend')
    grid on
    hold on

    figure(4)
    plot(mod_ratio_var,epsilon,'displayname',['enr=',num2str(enr*100),'%'])
    title('Effect of moderation ratio \Theta and enrichment on \epsilon')
    xlabel('\Theta (-)')
    ylabel('\epsilon (-)')
    legend('-dynamiclegend')
    grid on
    hold on

    % Plot k_inf and k_eff as a function of moderation ratio and enrichment
    figure(5)
    plot(mod_ratio_var,k_inf,'displayname',['enr=',num2str(enr*100),'%'])
    title('Infinite multiplication factor as a function of the enrichment and the moderation ratio')
    xlabel('\Theta (-)')
    ylabel('k_{\infty} (-)')
    legend('-dynamiclegend')
    grid on
    box on
    hold on
    
    figure(6)
    plot(mod_ratio_var,k_eff,'--','displayname',['enr=',num2str(enr*100),'%'])
    title('Effective multiplication factor as a function of the enrichment and the moderation ratio')
    xlabel('\Theta (-)')
    ylabel('k_{eff} (-)')
    legend('-dynamiclegend')
    grid on
    box on
    hold on

end

%% Post Processing
figure(6)
yline(1.03,'LineWidth',1.5,'displayname','initial k_{eff}')