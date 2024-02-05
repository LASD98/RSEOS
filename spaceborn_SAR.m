clear
close all
clc

%% HW 2-1

%% Data

% given
pho_x=10; % azimuth resolution (processing a full synthetic aperture) [m]
f0=435e6; % carrier frequency [Hz]
B=6e6; % transmitted bandwidth (maximum allowed by ITU regulations for frequency allocation) [Hz]
H=666000; % orbital height [m]
theta=28/180*pi;  % Reference incidence angle (o?-nadir) [°] 
T=290; % Reference noise temperature [K]
AGBmin=250;
AGBmax=500; % Reference range of AGB [tons per hectare]
sigma_hh_dB=-8; % Reference value for signal backscatter in HH polarization [dB]
             % (assuming the signal is dominated by trunkground double bounce scattering)
c = 3e8; % propagation velocity [m/s]

% obtained
lambda=c/f0; % wavelength [m]
pho_r=c/2/B; % range resolution [m] 
r=H/cos(theta); % reference satellite-target distance [m]
             
%% AGB estimation law

% AGB = A*sigma_hv^(p)

% sigma_hv is the backscatter coefficient of vegetation in HV polarization, 
% A and p are empirically derived constants

A=2.4e6;
p=2.6;

%% 1. Derive the relative accuracy to be required on sigma_hv to ensure an accuracy of 20% on AGB

% acc_x = dx/x
% AGB = A*sigma_hv^(p) -> dAGB = A*p*sigma_HV^(p-1)*dsigma_hv
% -> dAGB/AGB = (A*p*sigma_HV^(p-1)*dsigma_hv)/(A*sigma_hv^(p))
% -> acc_AGB = p*acc_sigma_hv = 20% -> 
acc_sigma_hv = 0.2/p;

%% 2. Determine the size of the averaging window to meet the requirement derived at point 1

% acc_sigma_hv = 1/sqrt(Leq); Leq = equivalent number of independent pixels
                                  % in the estimation window
Leq=(1/acc_sigma_hv)^2;

% Leq = L*sampling/resolution; L = real number of pixels in the estimation window
% we can assume L = Leq
L=Leq;

rcA=pho_r*pho_x; % resolution cell area [m^2]

% as we have 169 pixels, we can assume we have a 13x13 averaging window

Lx=sqrt(L)*pho_x;
Lr=sqrt(L)*pho_r;
avwA=Lx*Lr; % averaging window area [m^2]

%% 3. Determine a), b), c) as necessary to ensure a Signal to Noise Ratio 
% of at least 10 dB at HV polarization across the whole range of AGB values

SNR_dB=10;

sigma_hv_min=(AGBmin/A)^(1/p);
sigma_hv_max=(AGBmax/A)^(1/p);

ill_angle=lambda/2/pho_x;
Aill=r*pho_r*ill_angle;

% We have to make hypotesis on some antenna parameters to calculate the
% desired quantities
d=12; % antenna diameter [m] (real)
Ramb=50000; % ambiguous range [m]
RCS=sigma_hv_min*Aill; % radar cross section [m^2] 
RCSmax=sigma_hv_max*Aill;
atmL_dB=3; % 2 way atmospheric losses [dB]
F_dB = 3; % [dB] noise figure 
duty_cycle=0.1; % duty cycle (real)

% (a) Pulse Repetition Interval (PRI) 

PRI=2/c*Ramb;
PRF=1/PRI; % pulse repetition frequency

% (b) Duration of the transmitted pulses (=Tobs)

Tobs=PRI*duty_cycle;

% (c) Transmitted power for each transmitted pulse

A=pi*(d/2)^2; % antenna area [m^2]
G=4*pi/lambda^2*A; % antenna gain
K = 1.38e-23; % Boltzmann constant 
F = 10^(F_dB/10); % noise figure 
N0 = K*T*F; % noise spectral density 
atmL=10^(atmL_dB/10); % atmospheric losses
SNR=10^(SNR_dB/10); % required SNR

% Pr = Pt*(1/((4*pi*R^2)^2)*G*A*RCS/L) 
% Er = Pr*Toss = Pt*(Toss/((4*pi*R^2)^2)*G*A*RCS/L) 
% SNR = Er/N0 = Pt*(Toss/((4*pi*R^2)^2)*G*A*RCS/L/N0) 
Pt_W = SNR/(Tobs/((4*pi*r^2)^2)*G*A*RCS/atmL/N0); % transmitted power [Watt]
Pt_dBW = 10*log10(Pt_W); % [dBW = dB w.r.t. 1 W]
SNRmax = Pt_W*(Tobs/((4*pi*r^2)^2)*G*A*RCSmax/atmL/N0);
SNRmax_dB = 10*log10(SNRmax);

%% 4. Discuss whether your response at points 1) and 2) needs to be modi?ed when accounting for thermal noise



%% 5. Discuss the choice of the focusing processor (is 1D azimuth matched filtering enough?)



%% 6. One of BIOMASS secondary goals is to map terrain topography using SAR Interferometry. 
% By virtue of the penetration capabilities of P-Band waves, one can assume that the signal in HH polarization 
% is mostly contributed by trunk-ground double bounce scattering from several tree trunks within any resolution cell. 
% Accordingly, the phase center can be (at least approximately) assumed to be directly associated with terrain topography:

% (a) Evaluate the amount of phase noise in HH interferograms for a height of ambiguity z_amb = 80 m 

zamb=80;
kz=2*pi/zamb;
sigma_hh=10^(sigma_hh_dB/10);

sigma_delta_phi=sigma_hh*kz;
sigma_delta_phi_dB = 10*log10(sigma_delta_phi);

%Bn=lambda*r*sin(theta)/(2*zamb); %perpendicular baseline

% (b) Determine the size of the averaging window required to estimate topography to within an accuracy of 5 m

deltaR=5;

SNR6=1/sigma_delta_phi;

gamma=(SNR6/(1+SNR6))*exp(1i*4*pi*deltaR/lambda);

gamma=norm(gamma);

L6=(1-gamma^2)/(2*sigma_delta_phi^2*gamma^2);
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                