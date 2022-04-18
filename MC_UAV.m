clc
clear
close all

% system parameters

fc_b = 2;          % carrier frequency of BS (GHz)
fc_u = 2;          % carrier frequency of UAV (GHz)
fc_h = 2;          % carrier frequency of HAPS (GHz)  
fc_s_S = 2;        % carrier frequency of satellite in S-band (GHz)
fc_s_Ka = 30;      % carrier frequency of satellite in Ka-band (GHz)

BW_RB = 0.2e6;       % bandwidth of each PRB (Hz)
Nc = 2;              % number of PRB allocated to each user (allocated bandwidth): {1,2}
Nf = 2;              % frequency diversity factor for BU links: {1,2}

BW_bu_tot = 10e6;    % total bandwidth of BU (Hz) ---> 50 PRB
BW_bu = Nc*BW_RB;    % bandwidth of each BU link (Hz) < coherence BW

BW_uu_tot = 10e6;    % total bandwidth of UU (Hz)
BW_uu = Nc*BW_RB;    % bandwidth of each UU link (Hz)

BW_hu_tot = 10e6;    % total bandwidth of HU (Hz)
BW_hu = Nc*BW_RB;    % bandwidth of each HU link (Hz)

BW_su_tot = 10e6;    % total bandwidth of SU (Hz)
BW_su = Nc*BW_RB;    % bandwidth of each SU link (Hz)

BW_bh = 0.5e6;       % bandwidth of BH link (Hz): maximum bandwidth of URLLC user = coherence BW
BW_gs = 0.5e6;       % bandwidth of GS link (Hz)

R_max = 10e6;        % data rate (bps)
R_stp = 100e3;

Pb = 10^(46/10);   % BS Tx power (mW)
Pu = 10^(23/10);   % UAV Tx power (mW)
Ph = 10^(46/10);   % HAPS Tx power (mW)
Pg = 10^(50/10);   % gateway Tx power (mW)
Ps = 10^(50/10);   % sat Tx power (mW)

thetaT = 102*pi/180;      % downtilt angle = 102
ge_max = 8;        % maximum element gain (dBi)
Ne = 8;            % number of elements of BS antenna
Gt_hmax = 32;      % HAPS max Tx antenna gain
Gr_h = 32;         % HAPS Rx antenna gain
Gt_u = 0;          % UAV Tx antenna gain
Gr_u = 0;          % UAV Rx antenna gain
Gt_smax = 38;      % Sat max Tx antenna gain
Gr_s = 38;         % Sat Rx antenna gain
Gt_g = 38;         % gateway Tx antenna gain

NF_b = 5;          % BS Rx noise figure (dB)
NF_u = 9;          % UAV Rx noise figure
NF_h = 5;          % HAPS Rx noise figure
NF_s = 5;          % Sat Rx noise figure

N0 = 10^(-174/10);  % single-side noise (mW/Hz)
packet_size = 256;  % bits (32 bytes)
Re = 6371e3;        % radius of Earth (m)
c = 3e8;            % speed of light (m/s)

lambda = 2e3;     % average arrival rate of packets (packet/s)
eff_bw = 7*lambda; % effective bandwidth (packet/s)
Dq_max = 0.3e-3;   % qeueuing delay bound (s)

hb = 25;           % height of BS (meter)
hu = 300;          % altitude of UAV
hh = 20e3;         % altitude of HAPS

hs = 1110e3;          % altitude of satellite (Starlink-I constellation)
ns = 4425;            % number of satellites
theta_0 = 15*pi/180;  % minimum elevation angle of receiver antenna (radian)

nb = 37;              % number of cells (BSs)
Pr_intrf = 0.01;       % Interference probability
ISD = 500;            % inter-site distance in UMa (m)
Rnet = 1.6e3;         % area for UAV distribution = sqrt(2*sqrt(3)/pi)*3*ISD (approximately)
nu = 10;              % number of UAVs
nu_r = 3;             % size of swarm of coordinated UAVs
r_bh = 5e3;           % 2D distance of the HAP from its ground Station
d_2D = 150;           % 2D distance of the UAV from the serving BS (cell radius = sqrt(2*sqrt(3)/pi)*ISD/2 = 260 m)

xyz_HAP = [0,0,hh];
xyz_sat = [0,0,hs];
xyz_gat = [300e3,0,hb];  % suburban

E_th = 1e-5;       % error threshold
D_th1 = 10e-3;     % first delay threshold
D_th2 = 50e-3;     % second delay threshold
n_smp = 1e3;       % number of channel realization
n_itr = 1e1;       % number of iterations
Eb = 1e-6;         % backhual failure probability
Db = 1e-3;         % backhual delay

[Eq_b, Dq_b] = queue(eff_bw, Dq_max, lambda); % queueing delay violation
[Eq_u, Dq_u] = queue(eff_bw, Dq_max, lambda);

Gr_u = Gr_u - NF_u;
Gr_u = 10^(Gr_u/10);
Gt_u = 10^(Gt_u/10);  % omnidirectional
Gr_h = Gr_h - NF_h;
Gr_h = 10^(Gr_h/10);
Gr_s = Gr_s - NF_s;
Gr_s = 10^(Gr_s/10);
Gt_g = 10^(Gt_g/10);

ii = 1;
for R = [100e3:R_stp:R_max]
R_bu = R; R_uu = R;

for jj = 1:n_itr   % average over position
    
% random distribution of UAVs
[ xyz_BS, xyz_UAV ] = position( nb, hb, ISD, Rnet, nu, hu, d_2D );

d_uu = sqrt(sum((xyz_UAV(2:end,:) - repmat(xyz_UAV(1,:),[nu-1,1]))'.^2));    % distance of UAVs
[d_uu,index] = sort(d_uu,'ascend');
xyz_relay = xyz_UAV(2:end,:);
xyz_UAV = [xyz_UAV(1,:); xyz_relay(index,:)];

% BU and UU links
for i = 1:nu_r+1
    d_bu(i,:) = sqrt(sum((xyz_BS - repmat(xyz_UAV(i,:),[nb,1]))'.^2));    % 3D distance of BSs and UAVs
    r_bu(i,:) = sqrt(sum((xyz_BS(:,[1 2]) - repmat(xyz_UAV(i,[1 2]),[nb,1]))'.^2));    % 2D horizontal distance of BSs and UAVs
end

d_bu = sort(d_bu,2,'ascend');
r_bu = sort(r_bu,2,'ascend');
theta_bu = atan((hu-hb)./r_bu);

PL_bu_los = 28 + 20*log10(fc_b) + 22*log10(d_bu);
PL_bu_nlos = -17.5 + 20*log10(40*pi*fc_b/3) + (46 - 7*log10(hu))*log10(d_bu);

a1 = 0.3; a2 = 500; a3 = 20;   % for an urban scenario
K = floor(r_bu*sqrt(a1*a2)/1000-1);
for i = 1:nu_r+1
    for j = 1:nb
        Pr_los = 1;
        for k = 0:K(i,j)
            Pr_elem = 1-exp(-(hb-(k+0.5)*(hb-hu)/(K(i,j)+1))^2/(2*a3^2));
            Pr_los = Pr_elem*Pr_los;
        end
        Pr_los_bu(i,j) = Pr_los;
    end
end

PL_bu_dB = Pr_los_bu.*PL_bu_los + (1-Pr_los_bu).*PL_bu_nlos;
PL_uu_dB = 32.45 + 20*log10(fc_u) + 20*log10(d_uu);

PL_bu = 10.^(PL_bu_dB/10);
PL_uu = 10.^(PL_uu_dB/10);

PL_bu = repmat(PL_bu,[Nf,1]);

% rician channel generation
% rice factor of BU channel
for i = 1:nu_r+1
    for j = 1:nb
        if theta_bu(i,j) < 10*pi/180
            K_bu(i,j) = 10^(4/10);
        elseif theta_bu(i,j) >= 10*pi/180 && theta_bu(i,j) < 20*pi/180
            K_bu(i,j) = 10^(5/10);
        elseif theta_bu(i,j) >= 20*pi/180 && theta_bu(i,j) < 30*pi/180
            K_bu(i,j) = 10^(6/10);
        elseif theta_bu(i,j) >= 30*pi/180 && theta_bu(i,j) < 40*pi/180
            K_bu(i,j) = 10^(7/10);
        elseif theta_bu(i,j) >= 40*pi/180 && theta_bu(i,j) < 50*pi/180
            K_bu(i,j) = 10^(8/10);
        elseif theta_bu(i,j) >= 50*pi/180 && theta_bu(i,j) < 60*pi/180
            K_bu(i,j) = 10^(9/10);
        elseif theta_bu(i,j) >= 60*pi/180 && theta_bu(i,j) < 70*pi/180
            K_bu(i,j) = 10^(10/10);
        elseif theta_bu(i,j) >= 70*pi/180 && theta_bu(i,j) < 80*pi/180
            K_bu(i,j) = 10^(11/10);
        else
            K_bu(i,j) = 10^(12/10);
        end
    end
end
K_uu = 10^(12/10);    % rice factor of UU channel

sigma_bu = 1./sqrt(2*(K_bu+1));  % NLoS
sigma_uu = 1./sqrt(2*(K_uu+1));  % NLoS

rho_bu = sqrt(K_bu./(K_bu+1));   % LoS
rho_uu = sqrt(K_uu./(K_uu+1));   % LoS

sigma_bu = repmat(sigma_bu,[Nf,1]);
rho_bu = repmat(rho_bu,[Nf,1]);

for i = 1:Nf*(nu_r+1)
    for j = 1:nb
        h_bu(i,j,:) = rho_bu(i,j) + sigma_bu(i,j)*(randn(1,n_smp) + 1i*randn(1,n_smp));
    end
end
h_uu = rho_uu + sigma_uu*(randn(nu-1,n_smp) + 1i*randn(nu-1,n_smp));

g_bu = h_bu.*conj(h_bu);   % small-scale channel gain
g_uu = h_uu.*conj(h_uu);

% BS directional antenna
thetaZ_bu = pi/2 - theta_bu;    % zenith angle
ge = 10^(ge_max/10)*sin(thetaZ_bu).^2;     % element radiation pattern
gA = sin(Ne*pi*(cos(thetaZ_bu)-cos(thetaT))/2).^2./(Ne*sin(pi*(cos(thetaZ_bu)-cos(thetaT))/2).^2);     % array radiation pattern
Gt_b = ge.*gA;    % linear
Gt_b = repmat(Gt_b,[Nf,1]);

ICI = rand(1,nb);
I = ICI <= Pr_intrf; 
for i = 1:Nf*(nu_r+1)
    for j = 1:nb
        rx_pow_bu(i,j,:) = Pb*Gt_b(i,j)*Gr_u*g_bu(i,j,:)/PL_bu(i,j);
        BS_interf(i,j,:) = I(j)*rx_pow_bu(i,j,:);
    end
end

for i = 1:Nf*(nu_r+1)
    sinr_bu(i,:) = rx_pow_bu(i,1,:)./(BW_bu*N0 + sum(BS_interf(i,2:end,:),2));
end

IUI = rand(1,nu-1);
II = IUI <= Pr_intrf; 
for j = 1:nu-1
    rx_pow_uu(j,:) = Pu*Gt_u*Gr_u*g_uu(j,:)/PL_uu(j);
    UAV_interf(j,:) = II(j)*rx_pow_uu(j,:);
end

for i = 1:nu_r
    sinr_uu(i,:) = rx_pow_uu(i,:)./(BW_uu*N0 + sum(UAV_interf(nu_r+2:end,:),1));    % nu_r parallel UAVs with JD
end
% sinr_uu2 = sum(rx_pow_uu(1:2,:))./(BW_uu*N0 + sum(rx_pow_uu(3:end,:),1));   % 2 parallel UAVs with MRC

V_bu = 1-(1+sinr_bu).^(-2);
V_uu = 1-(1+sinr_uu).^(-2);
% V_uu2 = 1-(1+sinr_uu2).^(-2);

Dt_bu = packet_size/R_bu;  % transmission delay
Dt_uu = packet_size/R_uu;

et_bu = qfunc((BW_bu*log2(1+sinr_bu)-R_bu)*log(2)./sqrt(BW_bu*V_bu/Dt_bu));  % decoding error
et_uu = qfunc((BW_uu*log2(1+sinr_uu)-R_uu)*log(2)./sqrt(BW_uu*V_uu/Dt_uu));
% et_uu2 = qfunc((BW_uu*log2(1+sinr_uu2)-R_uu)*log(2)./sqrt(BW_uu*V_uu2/Dt_uu));

Etr_bu(:,k) = mean(et_bu,2);
Etr_uu(:,k) = mean(et_uu,2);

end

Et_bu = mean(Etr_bu,2);
Et_uu = mean(Etr_uu,2);

if Nf == 1
    Dtot_bu(ii) = Db + Dq_b + Dt_bu./(1-Et_bu(1));
    Etot_bu(ii) = 1-(1-Eb).*(1-Eq_b).*(1-Et_bu(1));
    
    Dtot_uu(:,ii) = Db + Dq_b + Dt_bu./(1-Et_bu(2:end)) + Dq_u + Dt_uu./(1-Et_uu);
    Etot_uu(:,ii) = 1-(1-Eb).*(1-Eq_b).*(1-Et_bu(2:end)).*(1-Eq_u).*(1-Et_uu);
else
    Dtot_bu(ii) = Db + Dq_b + min(Dt_bu./(1-Et_bu([1,nu_r+2])));
    Etot_bu(ii) = 1-(1-Eb).*(1-Eq_b).*(1-Et_bu(1).*Et_bu(nu_r+2));
    
    Dtot_uu(:,ii) = Db + Dq_b + Dt_bu./(1-[min(Et_bu([2,6]));min(Et_bu([3,7]));min(Et_bu([4,8]))]) + Dq_u + Dt_uu./(1-Et_uu);
    Etot_uu(:,ii) = 1-(1-Eb).*(1-Eq_b).*(1-Et_bu(2:nu_r+1).*Et_bu(nu_r+3:end)).*(1-Eq_u).*(1-Et_uu);
end 

ii = ii+1;
end
    
%% HAPS: BH and HU links

d_bh = sqrt(r_bh^2 + (hh-hb)^2);  % 3D distance of HAPS and its ground station
theta_bh = atan((hh-hb)./r_bh);

d_hu = sqrt(sum((xyz_HAP - xyz_UAV(1,:))'.^2));  % 3D distance of UAV and HAPS
r_hu = sqrt(sum((xyz_HAP(1,[1 2]) - xyz_UAV(1,[1 2]))'.^2)); 
theta_hu = atan((hh-hu)./r_hu);

Dp_bh = d_bh/c;
Dp_hu = d_hu/c;

[Eq_h, Dq_h] = queue(eff_bw, Dq_max, lambda);

ii = 1;
for R = [100e3:R_stp:R_max]
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

% Gateway directional antenna
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
Et_bh(ii) = mean(et_bh);

for q = 1:n_itr
IBI = rand(1,nb);
I3 = IBI <= Pr_intrf; 
for i = 1:nb
    rx_pow_hu(i,:) = Ph*Gt_h(i)*Gr_u*g_hu./PL_hu;
    HAP_interf(i,:) = I3(i)*rx_pow_hu(i,:);
end

sinr_hu = rx_pow_hu(1,:)./(BW_hu*N0 + sum(HAP_interf(2:end,:),1));
V_hu = 1-(1+sinr_hu).^(-2);

et_hu = qfunc((BW_hu*log2(1+sinr_hu)-R_hu)*log(2)./sqrt(BW_hu*V_hu/Dt_hu));
Etr_hu(q) = mean(et_hu);

dtot_hu(q) = Db + Dq_b + Dt_bh./(1-Et_bh(ii)) + Dp_bh + Dq_h + Dt_hu./(1-Etr_hu(q)) + Dp_hu;
etot_hu(q) = 1-(1-Eb).*(1-Eq_b).*(1-Et_bh(ii)).*(1-Eq_h).*(1-Etr_hu(q));

end

Et_hu(ii) = mean(Etr_hu,2);
Dtot_hu(ii) = mean(dtot_hu);
Etot_hu(ii) = mean(etot_hu);

ii = ii+1;
end

%% Satellite: GS and SU links
 
% % d_gs = sqrt((Re+hs)^2-Re^2*cos(theta_gs)^2) - Re*sin(theta_gs);
% % d_su = sqrt((Re+hs)^2-(Re+hu)^2*cos(theta_su)^2) - (Re+hu)*sin(theta_su);
% 
% d_gs = sqrt(sum((xyz_gat(1,:) - xyz_sat)'.^2));  % 3D distance of gateway and sat
% r_gs = sqrt(sum((xyz_gat(1,[1 2]) - xyz_sat(1,[1 2]))'.^2)); 
% theta_gs = atan((hs-hb)./r_gs);
% 
% d_su = sqrt(sum((xyz_sat - xyz_UAV(1,:))'.^2));  % 3D distance of UAV and sat
% r_su = sqrt(sum((xyz_sat(1,[1 2]) - xyz_UAV(1,[1 2]))'.^2)); 
% theta_su = atan((hs-hu)./r_su);
% 
% Dp_gs = d_gs/c;
% Dp_su = d_su/c;
% 
% [Eq_g, Dq_g] = queue(eff_bw, Dq_max, lambda);
% [Eq_s, Dq_s] = queue(eff_bw, Dq_max, lambda);
% 
% % probability of link unavailability 
% dmax_gs = sqrt((Re*sin(theta_0))^2 + hs^2 + 2*Re*hs) - Re*sin(theta_0);
% dmax_su = sqrt(((Re + hu)*sin(theta_0))^2 + (hs - hu)^2 + 2*(Re + hu)*(hs - hu)) - (Re + hu)*sin(theta_0);
% 
% Pvis_gs = 1-(1-(dmax_gs^2 - hs^2)/(4*Re*(Re + hs)))^ns;
% Pvis_su = 1-(1-(dmax_su^2 - (hs - hu)^2)/(4*(Re + hu)*(Re + hs)))^ns;
% 
% El_gs = 1 - Pvis_gs;
% El_su = 1 - Pvis_su;
% 
% ii = 1;
% for R = [100e3:R_stp:R_max]
% R_gs = R; R_su = R;
% 
% % Shadow fading and clutter loss for suburban scenario in S-band
% if theta_gs <= 10*pi/180
%     SF = [1.79*randn(1,0.782*n_smp), 8.93*randn(1,0.218*n_smp)];
%     CL = 0.218*19.52;
%     K_gs = 10^(5/10);
% elseif theta_gs > 10*pi/180 && theta_gs <= 20*pi/180
%     SF = [1.14*randn(1,0.869*n_smp), 9.08*randn(1,0.131*n_smp)];
%     CL = 0.131*18.17;
%     K_gs = 10^(5/10);
% elseif theta_gs > 20*pi/180 && theta_gs <= 30*pi/180
%     SF = [1.14*randn(1,0.919*n_smp), 8.78*randn(1,0.081*n_smp)];
%     CL = 0.081*18.42;
%     K_gs = 10^(6/10);
% elseif theta_gs > 30*pi/180 && theta_gs <= 40*pi/180
%     SF = [0.92*randn(1,0.929*n_smp), 10.25*randn(1,0.071*n_smp)];
%     CL = 0.071*18.28;
%     K_gs = 10^(7/10);
% elseif theta_gs > 40*pi/180 && theta_gs <= 50*pi/180
%     SF = [1.42*randn(1,0.935*n_smp), 10.56*randn(1,0.065*n_smp)];
%     CL = 0.065*18.63;
%     K_gs = 10^(8/10);
% elseif theta_gs > 50*pi/180 && theta_gs <= 60*pi/180
%     SF = [1.56*randn(1,0.94*n_smp), 10.74*randn(1,0.06*n_smp)];
%     CL = 0.06*17.68;
%     K_gs = 10^(9/10);
% elseif theta_gs > 60*pi/180 && theta_gs <= 70*pi/180
%     SF = [0.85*randn(1,0.949*n_smp), 10.17*randn(1,0.051*n_smp)];
%     CL = 0.051*16.5;
%     K_gs = 10^(11/10);
% elseif theta_gs > 70*pi/180 && theta_gs <= 80*pi/180
%     SF = [0.72*randn(1,0.998*n_smp), 11.52*randn(1,0.002*n_smp)];
%     CL = 0.002*16.3;
%     K_gs = 10^(15/10);
% else
%     SF = [0.72*randn(1,0.998*n_smp), 11.52*randn(1,0.002*n_smp)];
%     CL = 0.002*16.3;
%     K_gs = 10^(15/10);
% end
%                 
% PL_gs_dB = 32.45 + 20*log10(fc_s_S) + 20*log10(d_gs) + SF + CL;
% PL_su_dB = 32.45 + 20*log10(fc_s_S) + 20*log10(d_su);
% 
% PL_gs = 10.^(PL_gs_dB/10);
% PL_su = 10.^(PL_su_dB/10);
% 
% % rician channel generation
% if theta_su < 30*pi/180
%     K_su = 10^(12/10);
% elseif theta_su >= 30*pi/180 && theta_su < 50*pi/180
%     K_su = 10^(13/10);
% elseif theta_su >= 50*pi/180 && theta_su < 70*pi/180
%     K_su = 10^(14/10);
% else
%     K_su = 10^(15/10);
% end
% 
% sigma_gs = 1/sqrt(2*(K_gs+1));
% sigma_su = 1/sqrt(2*(K_su+1));
% 
% rho_gs = sqrt(K_gs/(K_gs+1));
% rho_su = sqrt(K_su/(K_su+1));
% 
% h_gs = rho_gs + sigma_gs*(randn(1,n_smp) + 1i*randn(1,n_smp));
% h_su = rho_su + sigma_su*(randn(1,n_smp) + 1i*randn(1,n_smp));
% 
% g_gs = h_gs.*conj(h_gs);
% g_su = h_su.*conj(h_su);
% 
% % Gateway directional antenna
% thetaZ_gs = pi/2 - theta_gs;
% ge = 10^(ge_max/10)*sin(thetaZ_gs).^2;
% gA = sin(Ne*pi*(cos(thetaZ_gs)-cos(thetaT))/2).^2./(Ne*sin(pi*(cos(thetaZ_gs)-cos(thetaT))/2).^2);
% Gt_g = ge.*gA;
% 
% % Satellite directional antenna
% for i = 1:nb
%     d_sc = sqrt(sum((xyz_sat - [xyz_BS(i,[1,2]),0])'.^2));  % 3D distance of Sat and cell i 
%     d_uc = sqrt(sum((xyz_UAV(1,:) - [xyz_BS(i,[1,2]),0])'.^2));  % 3D distance of UAV and cell i 
%     thetaB_su = acos((d_su^2 + d_sc^2 - d_uc^2)/(2*d_su*d_sc));  % angle between UAV and boresight of cell i
%     if thetaB_su == 0
%         Gt_s(i) = 1*10^(Gt_smax/10);
%     else
%         Gt_s(i) = 10^(Gt_smax/10)*4*abs(besselj(1,20*pi*sin(thetaB_su))./(20*pi*sin(thetaB_su))).^2;  % linear
%     end
% end
% 
% snr_gs_S = (Pg*Gt_g*Gr_s*g_gs./PL_gs)/(BW_gs*N0);
% V_gs = 1-(1+snr_gs_S).^(-2);
% 
% Dt_gs = packet_size/R_gs;
% Dt_su = packet_size/R_su;
% 
% et_gs = qfunc((BW_gs*log2(1+snr_gs_S)-R_gs)*log(2)./sqrt(BW_gs*V_gs/Dt_gs));
% Et_gs_S(ii) = mean(et_gs);
% 
% for q = 1:n_itr
% IBI = rand(1,nb);
% I3 = IBI <= Pr_intrf; 
% for i = 1:nb
%     rx_pow_su(i,:) = Ps*Gt_s(i)*Gr_u*g_su./PL_su;
%     Sat_interf(i,:) = I3(i)*rx_pow_su(i,:);
% end
% 
% sinr_su_S = rx_pow_su(1,:)./(BW_su*N0 + sum(Sat_interf(2:end,:),1));
% 
% V_su = 1-(1+sinr_su_S).^(-2);
% et_su = qfunc((BW_su*log2(1+sinr_su_S)-R_su)*log(2)./sqrt(BW_su*V_su/Dt_su));
% 
% Etr_su_S(q) = mean(et_su);
% 
% dtot_su_S(q) = Db + Dq_g + Dt_gs./(1-Et_gs_S(ii)) + Dp_gs + Dq_s + Dt_su./(1-Etr_su_S(q)) + Dp_su;
% etot_su_S(q) = 1-(1-Eb).*(1-Eq_g).*(1-El_gs).*(1-Et_gs_S(ii)).*(1-Eq_s).*(1-El_su).*(1-Etr_su_S(q));
% end
% 
% Et_su_S(ii) = mean(Etr_su_S,2);
% Dtot_su_S(ii) = mean(dtot_su_S);
% Etot_su_S(ii) = mean(etot_su_S);
% 
% ii = ii+1;
% 
% end
% 
% ii = 1;
% for R = [100e3:R_stp:R_max]
% R_gs = R; R_su = R;
% 
% % Shadow fading and clutter loss for suburban scenario in Ka-band
% if theta_gs <= 10*pi/180
%     SF = [1.9*randn(1,0.782*n_smp), 10.7*randn(1,0.218*n_smp)];
%     CL = 0.218*29.5;
%     K_gs = 10^(10/10);
% elseif theta_gs > 10*pi/180 && theta_gs <= 20*pi/180
%     SF = [1.6*randn(1,0.869*n_smp), 10.0*randn(1,0.131*n_smp)];
%     CL = 0.131*24.6;
%     K_gs = 10^(10/10);
% elseif theta_gs > 20*pi/180 && theta_gs <= 30*pi/180
%     SF = [1.9*randn(1,0.919*n_smp), 11.2*randn(1,0.081*n_smp)];
%     CL = 0.081*21.9;
%     K_gs = 10^(12/10);
% elseif theta_gs > 30*pi/180 && theta_gs <= 40*pi/180
%     SF = [2.3*randn(1,0.929*n_smp), 11.6*randn(1,0.071*n_smp)];
%     CL = 0.071*20.0;
%     K_gs = 10^(15/10);
% elseif theta_gs > 40*pi/180 && theta_gs <= 50*pi/180
%     SF = [2.7*randn(1,0.935*n_smp), 11.8*randn(1,0.065*n_smp)];
%     CL = 0.065*18.7;
%     K_gs = 10^(18/10);
% elseif theta_gs > 50*pi/180 && theta_gs <= 60*pi/180
%     SF = [3.1*randn(1,0.94*n_smp), 10.8*randn(1,0.06*n_smp)];
%     CL = 0.06*17.8;
%     K_gs = 10^(21/10);
% elseif theta_gs > 60*pi/180 && theta_gs <= 70*pi/180
%     SF = [3*randn(1,0.949*n_smp), 10.8*randn(1,0.051*n_smp)];
%     CL = 0.051*17.2;
%     K_gs = 10^(24/10);
% elseif theta_gs > 70*pi/180 && theta_gs <= 80*pi/180
%     SF = [3.6*randn(1,0.998*n_smp), 10.8*randn(1,0.002*n_smp)];
%     CL = 0.002*16.9;
%     K_gs = 10^(27/10);
% else
%     SF = [0.4*randn(1,0.998*n_smp), 10.8*randn(1,0.002*n_smp)];
%     CL = 0.002*16.8;
%     K_gs = 10^(30/10);
% end
%                 
% PL_gs_dB = 32.45 + 20*log10(fc_s_Ka) + 20*log10(d_gs) + SF + CL;
% PL_su_dB = 32.45 + 20*log10(fc_s_Ka) + 20*log10(d_su);
% 
% PL_gs = 10.^(PL_gs_dB/10);
% PL_su = 10.^(PL_su_dB/10);
% 
% % rician channel generation
% if theta_su <= 10*pi/180
%     K_su = 10^(20/10);
% elseif theta_su > 10*pi/180 && theta_su <= 20*pi/180
%     K_su = 10^(20/10);
% elseif theta_su > 20*pi/180 && theta_su <= 30*pi/180
%     K_su = 10^(21/10);
% elseif theta_su > 30*pi/180 && theta_su <= 40*pi/180
%     K_su = 10^(22/10);
% elseif theta_su > 40*pi/180 && theta_su <= 50*pi/180
%     K_su = 10^(23/10);
% elseif theta_su > 50*pi/180 && theta_su <= 60*pi/180
%     K_su = 10^(24/10);
% elseif theta_su > 60*pi/180 && theta_su <= 70*pi/180
%     K_su = 10^(26/10);
% elseif theta_su > 70*pi/180 && theta_su <= 80*pi/180
%     K_su = 10^(28/10);
% else
%     K_su = 10^(30/10);
% end
% 
% sigma_gs = 1/sqrt(2*(K_gs+1));
% sigma_su = 1/sqrt(2*(K_su+1));
% 
% rho_gs = sqrt(K_gs/(K_gs+1));
% rho_su = sqrt(K_su/(K_su+1));
% 
% h_gs = rho_gs + sigma_gs*(randn(1,n_smp) + 1i*randn(1,n_smp));
% h_su = rho_su + sigma_su*(randn(1,n_smp) + 1i*randn(1,n_smp));
% 
% g_gs = h_gs.*conj(h_gs);
% g_su = h_su.*conj(h_su);
% 
% snr_gs_Ka = (Pg*Gt_g*Gr_s*g_gs./PL_gs)/(BW_gs*N0);
% V_gs = 1-(1+snr_gs_Ka).^(-2);
% 
% Dt_gs = packet_size/R_gs;
% Dt_su = packet_size/R_su;
% 
% et_gs = qfunc((BW_gs*log2(1+snr_gs_Ka)-R_gs)*log(2)./sqrt(BW_gs*V_gs/Dt_gs));
% Et_gs_Ka(ii) = mean(et_gs);
% 
% for q = 1:n_itr
% IBI = rand(1,nb);
% I3 = IBI <= Pr_intrf; 
% for i = 1:nb
%     rx_pow_su(i,:) = Ps*Gt_s(i)*Gr_u*g_su./PL_su;
%     Sat_interf(i,:) = I3(i)*rx_pow_su(i,:);
% end
% 
% sinr_su_Ka = rx_pow_su(1,:)./(BW_su*N0 + sum(Sat_interf(2:end,:),1));
% V_su = 1-(1+sinr_su_Ka).^(-2);
% 
% et_su = qfunc((BW_su*log2(1+sinr_su_Ka)-R_su)*log(2)./sqrt(BW_su*V_su/Dt_su));
% 
% Etr_su_Ka(q) = mean(et_su);
% 
% dtot_su_Ka(q) = Db + Dq_g + Dt_gs/(1-Et_gs_Ka(ii)) + Dp_gs + Dq_s + Dt_su/(1-Etr_su_Ka(q)) + Dp_su;
% etot_su_Ka(q) = 1-(1-Eb)*(1-Eq_g)*(1-El_gs)*(1-Et_gs_Ka(ii))*(1-Eq_s)*(1-El_su)*(1-Etr_su_Ka(q));
% 
% end
% 
% Et_su_Ka(ii) = mean(Etr_su_Ka,2);
% Dtot_su_Ka(ii) = mean(dtot_su_Ka);
% Etot_su_Ka(ii) = mean(etot_su_Ka);
% 
% ii = ii+1;
% 
% end

%% ------------------Combining schemes---------------------

% combining A2A links
for k = 1:nu_r  
    Dtot_Kuu(k,:) = min(Dtot_uu([1:k],:),[],1);
    Etot_Kuu(k,:) = prod(Etot_uu([1:k],:),1); 
end
% combining DA2G and A2A links
Dtot_BU1 = min([Dtot_bu; Dtot_Kuu(1,:)]);  
Etot_BU1 = Etot_bu .* Etot_Kuu(1,:);

Dtot_BU2 = min([Dtot_bu; Dtot_Kuu(2,:)]);
Etot_BU2 = Etot_bu .* Etot_Kuu(2,:);

Dtot_BU3 = min([Dtot_bu; Dtot_Kuu(3,:)]);
Etot_BU3 = Etot_bu .* Etot_Kuu(3,:);

% combining DA2G and A2A and HAP links
% Dtot_BH = min([Dtot_bu; Dtot_hu]);
% Etot_BH = Etot_bu .* Etot_hu;
% 
% Dtot_BHU1 = min([Dtot_bu; Dtot_hu; Dtot_Kuu(1,:)]);
% Etot_BHU1 = Etot_bu .* Etot_hu .* Etot_Kuu(1,:);
% 
% Dtot_BHU2 = min([Dtot_bu; Dtot_hu; Dtot_Kuu(2,:)]);
% Etot_BHU2 = Etot_bu .* Etot_hu .* Etot_Kuu(2,:);
% 
% Dtot_BHU3 = min([Dtot_bu; Dtot_hu; Dtot_Kuu(3,:)]);
% Etot_BHU3 = Etot_bu .* Etot_hu .* Etot_Kuu(3,:);


figure(1)
R = [100e3:R_stp:R_max]/1e3;
loglog(R,Etot_bu,'-','linewidth',1), hold on, grid on
loglog(R,mean(Etot_uu),'-','linewidth',1)
% loglog(R,Etot_hu,'-','linewidth',1)
% loglog(R,Etot_su_S,'-','linewidth',1)
% loglog(R,Etot_su_Ka,'-','linewidth',1)
loglog(R,Etot_BU1,'--','linewidth',2)
loglog(R,Etot_BU2,'--','linewidth',2)
loglog(R,Etot_BU3,'--','linewidth',2)
% loglog(R,Etot_BH,'--','linewidth',2)
% loglog(R,Etot_BHU1,'--','linewidth',2)
% loglog(R,Etot_BHU2,'--','linewidth',2)
% loglog(R,Etot_BHU3,'--','linewidth',2)
loglog(R,repmat(E_th,size(R)),'k--','LineWidth',1)
xlabel('Data rate (kbps)')
ylabel('Overall error probability')
xlim([100e3 R_max]/1e3)
legend('DA2G','A2A','DA2G + 1-A2A','DA2G + 2-A2A','DA2G + 3-A2A')
% legend('DA2G','HAP','DA2G + 1-A2A','DA2G + 2-A2A','DA2G + 3-A2A','DA2G + HAP','DA2G + 1-A2A + HAP','DA2G + 2-A2A + HAP','DA2G + 3-A2A + HAP')

figure(2)
loglog(R,Dtot_bu,'-','linewidth',1), hold on, grid on
loglog(R,mean(Dtot_uu),'-','linewidth',1)
% loglog(R,Dtot_hu,'-','linewidth',1)
% loglog(R,Dtot_su_S,'-','linewidth',1)
% loglog(R,Dtot_su_Ka,'-','linewidth',1)
loglog(R,Dtot_BU1,'--','linewidth',2)
loglog(R,Dtot_BU2,'--','linewidth',2)
loglog(R,Dtot_BU3,'--','linewidth',2)
% loglog(R,Dtot_BH,'--','linewidth',2)
% loglog(R,Dtot_BHU1,'--','linewidth',2)
% loglog(R,Dtot_BHU2,'--','linewidth',2)
% loglog(R,Dtot_BHU3,'--','linewidth',2)
loglog(R,repmat(D_th1,size(R)),'k--','LineWidth',1)
loglog(R,repmat(D_th2,size(R)),'b--','LineWidth',1)
xlim([100e3 R_max]/1e3)
ylim([1e-3 1])
xlabel('Data rate (kbps)')
ylabel('E2E delay (sec)')
legend('DA2G','A2A','DA2G + 1-A2A','DA2G + 2-A2A','DA2G + 3-A2A')
% legend('DA2G','A2A','HAP','DA2G + 1-A2A','DA2G + 2-A2A','DA2G + 3-A2A','DA2G + HAP','DA2G + 1-A2A + HAP','DA2G + 2-A2A + HAP','DA2G + 3-A2A + HAP')
