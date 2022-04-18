function [ dtot_su, etot_su ] = Sat_com_Kaband( xyz_BS, xyz_gat, xyz_sat, xyz_UAV, fc_s, BW_su, BW_gs, R, Pg, Ps, thetaT, Gt_smax, Ne, ...
    Gr_s, Gt_g, Gr_u, Re, theta_0, ns, N0, packet_size, hb, hu, hs, nb, Pr_intrf, n_smp, Eb, Db, Eq_g, Dq_g, Eq_s, Dq_s, c )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

d_gs = sqrt(sum((xyz_gat(1,:) - xyz_sat)'.^2));  % 3D distance of gateway and sat
r_gs = sqrt(sum((xyz_gat(1,[1 2]) - xyz_sat(1,[1 2]))'.^2)); 
theta_gs = atan((hs-hb)./r_gs);

d_su = sqrt(sum((xyz_sat - xyz_UAV(1,:))'.^2));  % 3D distance of UAV and sat
r_su = sqrt(sum((xyz_sat(1,[1 2]) - xyz_UAV(1,[1 2]))'.^2)); 
theta_su = atan((hs-hu)./r_su);

Dp_gs = d_gs/c;
Dp_su = d_su/c;

% probability of link unavailability 
dmax_gs = sqrt((Re*sin(theta_0))^2 + hs^2 + 2*Re*hs) - Re*sin(theta_0);
dmax_su = sqrt(((Re + hu)*sin(theta_0))^2 + (hs - hu)^2 + 2*(Re + hu)*(hs - hu)) - (Re + hu)*sin(theta_0);

Pvis_gs = 1-(1-(dmax_gs^2 - hs^2)/(4*Re*(Re + hs)))^ns;
Pvis_su = 1-(1-(dmax_su^2 - (hs - hu)^2)/(4*(Re + hu)*(Re + hs)))^ns;

El_gs = 1 - Pvis_gs;
El_su = 1 - Pvis_su;

R_gs = R; R_su = R;

% Shadow fading and clutter loss for suburban scenario in Ka-band
if theta_gs <= 10*pi/180
    SF = [1.9*randn(1,0.782*n_smp), 10.7*randn(1,0.218*n_smp)];
    CL = 0.218*29.5;
    K_gs = 10^(10/10);
elseif theta_gs > 10*pi/180 && theta_gs <= 20*pi/180
    SF = [1.6*randn(1,0.869*n_smp), 10.0*randn(1,0.131*n_smp)];
    CL = 0.131*24.6;
    K_gs = 10^(10/10);
elseif theta_gs > 20*pi/180 && theta_gs <= 30*pi/180
    SF = [1.9*randn(1,0.919*n_smp), 11.2*randn(1,0.081*n_smp)];
    CL = 0.081*21.9;
    K_gs = 10^(12/10);
elseif theta_gs > 30*pi/180 && theta_gs <= 40*pi/180
    SF = [2.3*randn(1,0.929*n_smp), 11.6*randn(1,0.071*n_smp)];
    CL = 0.071*20.0;
    K_gs = 10^(15/10);
elseif theta_gs > 40*pi/180 && theta_gs <= 50*pi/180
    SF = [2.7*randn(1,0.935*n_smp), 11.8*randn(1,0.065*n_smp)];
    CL = 0.065*18.7;
    K_gs = 10^(18/10);
elseif theta_gs > 50*pi/180 && theta_gs <= 60*pi/180
    SF = [3.1*randn(1,0.94*n_smp), 10.8*randn(1,0.06*n_smp)];
    CL = 0.06*17.8;
    K_gs = 10^(21/10);
elseif theta_gs > 60*pi/180 && theta_gs <= 70*pi/180
    SF = [3*randn(1,0.949*n_smp), 10.8*randn(1,0.051*n_smp)];
    CL = 0.051*17.2;
    K_gs = 10^(24/10);
elseif theta_gs > 70*pi/180 && theta_gs <= 80*pi/180
    SF = [3.6*randn(1,0.998*n_smp), 10.8*randn(1,0.002*n_smp)];
    CL = 0.002*16.9;
    K_gs = 10^(27/10);
else
    SF = [0.4*randn(1,0.998*n_smp), 10.8*randn(1,0.002*n_smp)];
    CL = 0.002*16.8;
    K_gs = 10^(30/10);
end
                
PL_gs_dB = 32.45 + 20*log10(fc_s) + 20*log10(d_gs) + SF + CL;
PL_su_dB = 32.45 + 20*log10(fc_s) + 20*log10(d_su);

PL_gs = 10.^(PL_gs_dB/10);
PL_su = 10.^(PL_su_dB/10);

% rician channel generation
if theta_su <= 10*pi/180
    K_su = 10^(20/10);
elseif theta_su > 10*pi/180 && theta_su <= 20*pi/180
    K_su = 10^(20/10);
elseif theta_su > 20*pi/180 && theta_su <= 30*pi/180
    K_su = 10^(21/10);
elseif theta_su > 30*pi/180 && theta_su <= 40*pi/180
    K_su = 10^(22/10);
elseif theta_su > 40*pi/180 && theta_su <= 50*pi/180
    K_su = 10^(23/10);
elseif theta_su > 50*pi/180 && theta_su <= 60*pi/180
    K_su = 10^(24/10);
elseif theta_su > 60*pi/180 && theta_su <= 70*pi/180
    K_su = 10^(26/10);
elseif theta_su > 70*pi/180 && theta_su <= 80*pi/180
    K_su = 10^(28/10);
else
    K_su = 10^(30/10);
end

sigma_gs = 1/sqrt(2*(K_gs+1));
sigma_su = 1/sqrt(2*(K_su+1));

rho_gs = sqrt(K_gs/(K_gs+1));
rho_su = sqrt(K_su/(K_su+1));

h_gs = rho_gs + sigma_gs*(randn(1,n_smp) + 1i*randn(1,n_smp));
h_su = rho_su + sigma_su*(randn(1,n_smp) + 1i*randn(1,n_smp));

g_gs = h_gs.*conj(h_gs);
g_su = h_su.*conj(h_su);

% Gateway directional antenna
% thetaZ_gs = pi/2 - theta_gs;
% ge = 10^(ge_max/10)*sin(thetaZ_gs).^2;
% gA = sin(Ne*pi*(cos(thetaZ_gs)-cos(thetaT))/2).^2./(Ne*sin(pi*(cos(thetaZ_gs)-cos(thetaT))/2).^2);
% Gt_g = ge.*gA;

% Satellite directional antenna
for i = 1:nb
    d_sc = sqrt(sum((xyz_sat - [xyz_BS(i,[1,2]),0])'.^2));  % 3D distance of Sat and cell i 
    d_uc = sqrt(sum((xyz_UAV(1,:) - [xyz_BS(i,[1,2]),0])'.^2));  % 3D distance of UAV and cell i 
    thetaB_su = acos((d_su^2 + d_sc^2 - d_uc^2)/(2*d_su*d_sc));  % angle between UAV and boresight of cell i
    if thetaB_su == 0
        Gt_s(i) = 1*10^(Gt_smax/10);
    else
        Gt_s(i) = 10^(Gt_smax/10)*4*abs(besselj(1,20*pi*sin(thetaB_su))./(20*pi*sin(thetaB_su))).^2;  % linear
    end
end

snr_gs = (Pg*Gt_g*Gr_s*g_gs./PL_gs)/(BW_gs*N0);
V_gs = 1-(1+snr_gs).^(-2);

Dt_gs = packet_size/R_gs;
Dt_su = packet_size/R_su;

et_gs = qfunc((BW_gs*log2(1+snr_gs)-R_gs)*log(2)./sqrt(BW_gs*V_gs/Dt_gs));

IBI = rand(1,nb);
I3 = IBI <= Pr_intrf; 
for i = 1:nb
    rx_pow_su(i,:) = Ps*Gt_s(i)*Gr_u*g_su./PL_su;
    Sat_interf(i,:) = I3(i)*rx_pow_su(i,:);
end

sinr_su = rx_pow_su(1,:)./(BW_su*N0 + sum(Sat_interf(2:end,:),1));
V_su = 1-(1+sinr_su).^(-2);

et_su = qfunc((BW_su*log2(1+sinr_su)-R_su)*log(2)./sqrt(BW_su*V_su/Dt_su));

dtot_su = Db + Dq_g + Dt_gs./(1-et_gs) + Dp_gs + Dq_s + Dt_su./(1-et_su) + Dp_su;
etot_su = 1-(1-Eb).*(1-Eq_g).*(1-et_gs).*(1-Eq_s).*(1-et_su);

end

