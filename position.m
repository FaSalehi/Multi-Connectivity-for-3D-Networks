function [ xyz_BS, xyz_UAV ] = position( nBS,hBS,ISD,R_net,nUAV,hUAV,d_2D )

% %   random position of BS
% R_net = 1000;
% rb = ISD+(R_net-ISD)*sqrt(rand(nBS-1,1));
% theta_b = 2*pi*rand(nBS-1,1);
% xb = rb.*cos(theta_b);
% yb = rb.*sin(theta_b);
% zb = hBS*ones(nBS-1,1);
% xyz_BS = [0,0,hBS; xb,yb,zb];

%   Hexagonal grid of BSs

a = sqrt(3)*ISD/2;
b = ISD/2;
xy_grid = [0,0; a,b; 0,2*b; -a,b; -a,-b; 0,-2*b; a,-b; ...
           2*a,0; 2*a,2*b; a,3*b; 0,4*b; -a,3*b; -2*a,2*b; -2*a,0; -2*a,-2*b; -a,-3*b; 0,-4*b; a,-3*b; 2*a,-2*b; ...
           3*a,b; 3*a,3*b; 2*a,4*b; a,5*b; 0,6*b; -a,5*b; -2*a,4*b; -3*a,3*b; -3*a,b; -3*a,-b; -3*a,-3*b; -2*a,-4*b; -a,-5*b; 0,-6*b; a,-5*b; 2*a,-4*b; 3*a,-3*b; 3*a,-b]; % 37 cells (3 tiers)

xyz_BS = [xy_grid(1:nBS,:),hBS*ones(nBS,1)];

% t=0:0.1:2*pi;
% a=R_net.*sin(t);
% b=R_net.*cos(t);
% plot(a,b,'k');
% hold on;
% plot(xy_grid(:,1),xy_grid(:,2),'s');

%   Distributing the UAVs
% xyz_u = [d_2D,0,hUAV];    % desired UAV
% random positioning
R_cell = sqrt(2*sqrt(3)/pi)*ISD/2;
r_u0 = R_cell*sqrt(rand(1,1));
theta_u0 = 2*pi*rand(1,1);
x_u0 = r_u0.*cos(theta_u0);
y_u0 = r_u0.*sin(theta_u0);
xyz_u = [x_u0,y_u0,hUAV];    % desired UAV

ru = R_net*sqrt(rand(nUAV-1,1));
theta_u = 2*pi*rand(nUAV-1,1);
xu = ru.*cos(theta_u);
yu = ru.*sin(theta_u);
zu = hUAV*ones(nUAV-1,1);
xyz_UAV = [xyz_u; xu,yu,zu];

end
