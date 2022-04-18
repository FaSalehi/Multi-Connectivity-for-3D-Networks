function [ dtot_hu, etot_hu ] = HAP_com( xyz_BS, xyz_UAV, xyz_HAP, fc_b, fc_h, BW_bh, BW_hu, R, Ph, Pb, thetaT, ge_max, Ne, Gt_hmax, ... 
    Gr_h, Gr_u, N0, packet_size, hb, hu, hh, nb, r_bh, Pr_intrf, n_smp, Eb, Db, Eq_b, Dq_b, Eq_h, Dq_h, c )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

d_bh = sqrt(r_bh^2 + (hh-hb)^2);  % 3D distance of HAPS and its ground station
theta_bh = atan((hh-hb)./r_bh);

d_hu = sqrt(sum((xyz_HAP - xyz_UAV(1,:))'.^2));  % 3D distance of UAV and HAPS
r_hu = sqrt(sum((xyz_HAP(1,[1 2]) - xyz_UAV(1,[1 2]))'.^2)); 
theta_hu = atan((hh-hu)./r_hu);

Dp_bh = d_bh/c;
Dp_hu = d_hu/c;

R_bh = R; R_hu = R;

% Shadow fading and clutter loss for urban scenario
if theta_bh <= 10*pi/180
    SF = [4*randn(1,0.246*n_smp), 6*randn(1,0.754*n_smp)];
    CL = 0.754*34.3;
    K_bh = 10^(4/10);
elseif theta_bh > 10*pi/180 && theta_bh <= 20*pi/180
    SF = [4*randn(1,0.386*n_smp), 6*randn(1,0.614*n_smp)];
    CL = 0.614*30.9;
    K_bh = 10^(5/10);
elseif theta_bh > 20*pi/180 && theta_bh <= 30*pi/180
    SF = [4*randn(1,0.493*n_smp), 6*randn(1,0.507*n_smp)];
    CL = 0.507*29;
    K_bh = 10^(6/10);
elseif theta_bh > 30*pi/180 && theta_bh <= 40*pi/180
    SF = [4*randn(1,0.613*n_smp), 6*randn(1,0.387*n_smp)];
    CL = 0.387*27.7;
    K_bh = 10^(7/10);
elseif theta_bh > 40*pi/180 && theta_bh <= 50*pi/180
    SF = [4*randn(1,0.726*n_smp), 6*randn(1,0.274*n_smp)];
    CL = 0.274*26.8;
    K_bh = 10^(8/10);
elseif theta_bh > 50*pi/180 && theta_bh <= 60*pi/180
    SF = [4*randn(1,0.805*n_smp), 6*randn(1,0.195*n_smp)];
    CL = 0.195*26.2;
    K_bh = 10^(9/10);
elseif theta_bh > 60*pi/180 && theta_bh <= 70*pi/180
    SF = [4*randn(1,0.919*n_smp), 6*randn(1,0.081*n_smp)];
    CL = 0.081*25.8;
    K_bh = 10^(11/10);
elseif theta_bh > 70*pi/180 && theta_bh <= 80*pi/180
    SF = [4*randn(1,0.968*n_smp), 6*randn(1,0.032*n_smp)];
    CL = 0.032*25.5;
    K_bh = 10^(13/10);
else
    SF = [4*randn(1,0.992*n_smp), 6*randn(1,0.008*n_smp)];
    CL = 0.008*25.5;
    K_bh = 10^(15/10);
end

PL_bh_dB = 32.45 + 20*log10(fc_b) + 20*log10(d_bh) + SF + CL;
PL_hu_dB = 32.45 + 20*log10(fc_h) + 20*log10(d_hu);

PL_bh = 10.^(PL_bh_dB/10);
PL_hu = 10.^(PL_hu_dB/10);

% rician channel generation
if theta_hu < 30*pi/180
    K_hu = 10^(12/10);
elseif theta_hu >= 30*pi/180 && theta_hu < 50*pi/180
    K_hu = 10^(13/10);
elseif theta_hu >= 50*pi/180 && theta_hu < 70*pi/180
    K_hu = 10^(14/10);
else
    K_hu = 10^(15/10);
end

sigma_bh = 1/sqrt(2*(K_bh+1));
sigma_hu = 1/sqrt(2*(K_hu+1));

rho_bh = sqrt(K_bh/(K_bh+1));
rho_hu = sqrt(K_hu/(K_hu+1));

h_bh = rho_bh + sigma_bh*(randn(1,n_smp) + 1i*randn(1,n_smp));
h_hu = rho_hu + sigma_hu*(randn(1,n_smp) + 1i*randn(1,n_smp));

g_bh = h_bh.*conj(h_bh);
g_hu = h_hu.*conj(h_hu);

% BS directional antenna
thetaZ_bh = pi/2 - theta_bh;
ge = 10^(ge_max/10)*sin(thetaZ_bh).^2;
gA = sin(Ne*pi*(cos(thetaZ_bh)-cos(thetaT))/2).^2./(Ne*sin(pi*(cos(thetaZ_bh)-cos(thetaT))/2).^2);
Gt_b = ge.*gA;

% HAP directional antenna
for i = 1:nb
    d_hc = sqrt(sum((xyz_HAP - [xyz_BS(i,[1,2]),0])'.^2));  % 3D distance of HAP and cell i 
    d_uc = sqrt(sum((xyz_UAV(1,:) - [xyz_BS(i,[1,2]),0])'.^2));  % 3D distance of UAV and cell i 
    thetaB_hu = acos((d_hu^2 + d_hc^2 - d_uc^2)/(2*d_hu*d_hc));  % angle between UAV and boresight of cell i
    if thetaB_hu == 0
        Gt_h(i) = 1*10^(Gt_hmax/10);
    else
        Gt_h(i) = 10^(Gt_hmax/10)*4*abs(besselj(1,20*pi*sin(thetaB_hu))./(20*pi*sin(thetaB_hu))).^2;  % linear
    end
end

snr_bh = (Pb*Gt_b*Gr_h*g_bh./PL_bh)./(BW_bh*N0);
V_bh = 1-(1+snr_bh).^(-2);

Dt_bh = packet_size/R_bh;
Dt_hu = packet_size/R_hu;

et_bh = qfunc((BW_bh*log2(1+snr_bh)-R_bh)*log(2)./sqrt(BW_bh*V_bh/Dt_bh));

IBI = rand(1,nb);
I3 = IBI <= Pr_intrf; 
for i = 1:nb
    rx_pow_hu(i,:) = Ph*Gt_h(i)*Gr_u*g_hu./PL_hu;
    HAP_interf(i,:) = I3(i)*rx_pow_hu(i,:);
end

sinr_hu = rx_pow_hu(1,:)./(BW_hu*N0 + sum(HAP_interf(2:end,:),1));
V_hu = 1-(1+sinr_hu).^(-2);

et_hu = qfunc((BW_hu*log2(1+sinr_hu)-R_hu)*log(2)./sqrt(BW_hu*V_hu/Dt_hu));

dtot_hu = Db + Dq_b + Dt_bh./(1-et_bh) + Dp_bh + Dq_h + Dt_hu./(1-et_hu) + Dp_hu;
etot_hu = 1-(1-Eb).*(1-Eq_b).*(1-et_bh).*(1-Eq_h).*(1-et_hu);

end


