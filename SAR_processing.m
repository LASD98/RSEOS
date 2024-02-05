clearvars -except SAR_data
% clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HW.2-2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
    SAR_data = Generate_raw_SAR_data('proton.jpg');
end
D = SAR_data.D;
t_ax = SAR_data.t_ax;
xa = SAR_data.xa;
g = SAR_data.g;
r_ax = SAR_data.r_ax;
x_ax = SAR_data.x_ax;
f0 = SAR_data.f0;
As = SAR_data.As;
B = SAR_data.B;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pho_r = 3
% 
% dxa = 2.5000
% 
% SAR_data = 
%   struct with fields:
% 
%                    D: [535×3003 double]
%        D_description: 'raw data matrix [fast time, sensor position]'
% 
%                 t_ax: [1×535 double]
%     t_ax_description: 'fast times [s]'
% 
%                   xa: [1×3003 double]
%       xa_description: 'sensor positions along orbit [m]'
% 
%                    g: [1×201 double]
%        g_description: 'transmitted pulse'
% 
%                 r_ax: [1×535 double]
%     r_ax_description: 'range positions of the targets [m]'
% 
%                 x_ax: [1×334 double]
%     x_ax_description: 'azimuth positions of the targets [m]'
% 
%                   f0: 5.0000e+09
%       f0_description: 'carrier frequency [Hz]'
% 
%                   As: 7005
%       As_description: 'Synthetic Aperture length [m]'
% 
%                    B: 50000000
%        B_description: 'Bandwidth [Hz]'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Radar and Synthetic Aperture Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 3e8; % wave propagation velocity [m/s]
f0 = f0; % carrier frequency [Hz]
lambda = c/f0; % wavelength [m]
B = B; % bandwidth [Hz]
pho_r = c/2/B; % range resolution [m]
r_min = t_ax(1)*c/2;
r_max = t_ax(end)*c/2;
As = As;
Lx = lambda/As*r_max; % maximum synthetic aperture
dx_ant = Lx/2*0.5; % spatial sampling of the synthetic aperture
Nx = length(xa);
dt = 1/2/B;
Nr = length(r_ax);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matched Filtering (Range Compression)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vec = (length(t_ax)-length(g))/2;
g = [zeros(1,vec),g,zeros(1,vec)];
Drc = zeros(Nr,Nx);
conj_flip_g = conj(fliplr(g))';
for i = 1:Nx 
   Drc(:,i) = conv2(D(:,i),conj_flip_g,'same')*dt;
end

figure
subplot(1,2,1)
imagesc(xa,r_ax,abs(D)), title('Raw Data')
xlabel('azimuth [m]'), ylabel('slant range [m]')
subplot(1,2,2)
imagesc(xa,r_ax,abs(Drc)), title('Range Compressed Data')
xlabel('azimuth [m]'), ylabel('slant range [m]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOCUSING BY 1D MATCHED FILTERING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
Dfr = zeros(Nr,Nx);
for k = 1:Nr
    R = sqrt(xa.^2 + r_ax(k).^2);
    filter = +exp(1i*4*pi/lambda*R); 
    Dfr(k,:) = conv2(Drc(k,:),filter,'same')*dt;  
end
t1DMF=toc;
figure
subplot(1,3,1), imagesc(xa,r_ax,abs(Dfr))
title('Focused data by 1D Matched filter')
xlabel('azimuth [m]'), ylabel('slant range [m]')
subplot(1,3,2), imagesc(xa,r_ax,abs(Dfr))
title('Focused data by 1D Matched filter (Zoomed In)')
xlabel('azimuth [m]'), ylabel('slant range [m]'), xlim([-255 255])
subplot(1,3,3), imagesc(xa,r_ax,abs(Dfr))
title('Along-track resolution')
xlabel('azimuth [m]'), ylabel('slant range [m]'), xlim([-4 4])

figure
subplot(1,2,1), imagesc(x_ax,r_ax,abs(Dfr(:,1402:1602)))
title('Focused data by 1D Matched filter')
xlabel('target azimuth [m]'), ylabel('slant range [m]'), xlim([-255 255])
subplot(1,2,2), imagesc(x_ax,r_ax,abs(Dfr(:,1402:1602)))
title('Along-track resolution')
xlabel('target azimuth [m]'), ylabel('slant range [m]'), xlim([-4 4])

figure
imagesc(x_ax,r_ax,abs(Dfr(:,1402:1602)))
title('resolution')
xlabel('target azimuth [m]'), ylabel('slant range [m]'), xlim([-4 4]), ylim([7.003e5 7.003e5+6])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOCUSING BY TIME DOMAIN BACK PROJECTION (TDBP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[R,X] = ndgrid(r_ax, xa);
wbar = waitbar(0,'Backprojecting');
Ns = round(As/2/dx_ant); 
S = 0;
tic
for n = 1:Nx
    waitbar(n/Nx,wbar)
    Rn_rx = sqrt( (X-xa(n)).^2 + R.^2);
    % to speed-up the computation, only the samples within a synthetic 
    % aperture are computed
    ind_x = max(1,n-Ns):min(Nx,n+Ns);
    Rn_rx = Rn_rx(:,ind_x);
    Sn = zeros(Nr,Nx);
    Sn(:,ind_x) = interp1(r_ax,Drc(:,n),Rn_rx).*exp(+1i*4*pi/lambda*Rn_rx); 
    S = S + Sn;
end
close(wbar)
tTDBP=toc;
figure
imagesc(xa,r_ax,abs(S)), title('Focused image by TDBP')
xlabel('azimuth [m]'), ylabel('slant range [m]')
figure
imagesc(xa,r_ax,abs(S))
title('Focused image by TDBP (Zoomed In)')
xlabel('azimuth [m]'), ylabel('slant range [m]'), xlim([-255 255])

figure
imagesc(xa,r_ax,abs(S))
title('resolution')
xlabel('target azimuth [m]'), ylabel('slant range [m]'), xlim([-4 4]), ylim([7.003e5 7.003e5+6])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


