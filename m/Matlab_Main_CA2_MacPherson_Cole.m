%% ASEN 3111 - Computational Assignment 2 - Main
% Model the flow characteristics over an airfoil using multiple vortices 
% along a given chord length by finding and calculating the stream 
% function, velocity potential, and pressure contour for the airfoil
%
% Author: Cole MacPherson
% Collaborators: R. Block, Z. Lesan, S. Mansfield, A. Uprety
% Date: 21st Feb 2021

%% Housekeeping

clc;
clear;
close all;
tic

%% Constant Declaration

c = 5; % chord length [m]
alpha = deg2rad(10); % angle of attack [degrees]
V_inf = 50; % free stream velocity [m/s]
P_inf = 101.3e3; % free stream pressure [Pa]
rho_inf = 1.225; % free stream density [kg/m^3]
N = 1000; % number of vortices

%% Call function to find S.F., V.P., and P.C. for the desired flow
% last input determines what runs in the function:
%   1 - regular function that calculates S.F., V.P., and P.C. for 50000
%           vortices and then finds the number of vortices to get less 
%           than a 0.5% error (Bullet point 2)
%   2 - regular function that is used to conducts study to determine 
%           changes in 1) chord length, 2) angle of attack, and 3) free 
%           stream flow speed (Bullet points 1 and 3)
%   3 - loads saved data from the running of number 1 (Bullet point 1 & 2 data)
[x,y,Psi,Phi,P,levels,bounds,N2,err] = Plot_Airfoil_Flow(c,alpha,V_inf,P_inf,rho_inf,N,3);

%% Plot stream function at levels
figure
contourf(x,y,Psi,75) % stream function
hold on 
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis equal
axis(bounds)
ylabel('y')
xlabel('x')
title(['Stream Lines for N = ' num2str(N2) ' Vorticies']);

%% Plot velocity potential at levels
figure
contourf(x,y,Phi,100) % velocity potential
hold on 
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis equal
axis(bounds)
ylabel('y')
xlabel('x')
title(['Equipotential Lines for N = ' num2str(N2) ' Vorticies']);

%% Plot pressure contours
figure
contourf(x,y,P,100) % pressure contour
hold on
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis equal
axis(bounds)
k = colorbar;
k.Label.String = 'Pressure [Pa]';
ylabel('y')
xlabel('x')
title(['Pressure Contours for N = ' num2str(N2) ' Vorticies']);

%% Plot error
N_spacing = N:1000:N2;
figure
plot(N_spacing,err(1:length(N_spacing))*100,'linewidth',2);
ylabel('Percent Error')
xlabel('Number of Vortices')
title(['Pressure Contours for N = ' num2str(N2) ' Vorticies']);

fprintf('figure 4 displays how the accuracy of the stream function, \nvelocity potential, and pressure contour changes as the number of\nvortices increases. It was found that %i is the number of\nvortices required to have a accuarcy less than 0.5%%. This\naccuracy was based off of a 50,000 vortice model being taken as the correct\ntheoretical value and error was determined from there. In the end,\nas the number of votices increases the better the model becomes at\nmodelling the system. Theoretically the best model would have an\ninfinite amount of vorticies but that is not fesible.\n\n',N2);

%% Study on how chord length changes S.F., V.P., and P.C.
% last input determines what runs in the function:
%   1 - regular function that calculates S.F., V.P., and P.C. for 50000
%           vortices and then finds the number of vortices to get less 
%           than a 0.5% error (Bullet point 2)
%   2 - regular function that is used to conducts study to determine 
%           changes in 1) chord length, 2) angle of attack, and 3) free 
%           stream flow speed (Bullet points 1 and 3)
%   3 - loads saved data from the running of number 1 (Bullet point 1 & 2 data)
N = 1000;

% chord length of 2
c2 = 2;
[x_c2,y_c2,Psi_c2,Phi_c2,~,~,bounds_c2,~] = Plot_Airfoil_Flow(c2,alpha,V_inf,P_inf,rho_inf,N,2);

% chord length of 5
c5 = 5;
[x_c5,y_c5,Psi_c5,Phi_c5,~,~,bounds_c5,~] = Plot_Airfoil_Flow(c5,alpha,V_inf,P_inf,rho_inf,N,2);

% chord length of 8
c8 = 8;
[x_c8,y_c8,Psi_c8,Phi_c8,~,~,bounds_c8,~] = Plot_Airfoil_Flow(c8,alpha,V_inf,P_inf,rho_inf,N,2);

% plot S.F.
figure
subplot(3,1,1)
contourf(x_c2,y_c2,Psi_c2,40) % stream function
hold on 
plot([0 c2],[0 0],'k','linewidth',2) % airfoil
axis(bounds_c2)
title(['Stream Lines for N = ' num2str(N) ' Vorticies at c = ' num2str(c2) ' [m]']);
subplot(3,1,2)
contourf(x_c5,y_c5,Psi_c5,40) % stream function
hold on
plot([0 c5],[0 0],'k','linewidth',2) % airfoil
axis(bounds_c5)
title(['Stream Lines for N = ' num2str(N) ' Vorticies at c = ' num2str(c5) ' [m]']);
ylabel('y')
subplot(3,1,3)
contourf(x_c8,y_c8,Psi_c8,40) % stream function
hold on
plot([0 c8],[0 0],'k','linewidth',2) % airfoil
axis(bounds_c8)
title(['Stream Lines for N = ' num2str(N) ' Vorticies at c = ' num2str(c8) ' [m]']);
xlabel('x')
hold off

% plot V.P.
figure
subplot(3,1,1)
contourf(x_c2,y_c2,Phi_c2,100) % velocity potential
hold on 
plot([0 c2],[0 0],'k','linewidth',2) % airfoil
axis(bounds_c2)
title(['Equipotential Lines for N = ' num2str(N) ' Vorticies at c = ' num2str(c2) ' [m]']);
subplot(3,1,2)
contourf(x_c5,y_c5,Phi_c5,100) % velocity potential
hold on
plot([0 c5],[0 0],'k','linewidth',2) % airfoil
axis(bounds_c5)
title(['Equipotential Lines for N = ' num2str(N) ' Vorticies at c = ' num2str(c5) ' [m]']);
ylabel('y')
subplot(3,1,3)
contourf(x_c8,y_c8,Phi_c8,100) % velocity potential
hold on
plot([0 c8],[0 0],'k','linewidth',2) % airfoil
axis(bounds_c8)
title(['Equipotential Lines for N = ' num2str(N) ' Vorticies at c = ' num2str(c8) ' [m]']);
xlabel('x')
hold off

fprintf('figures 5 & 6 display how the stream function and the velocity\npotential change as the chord length changes. When looking\nat the figures to view how the chord length affects the two\nit is seen that the chord length does not have a large effect\non the stream function or the velocity potential.\n\n');

%% Study on how angle of attack changes S.F., V.P., and P.C.
% last input determines what runs in the function:
%   1 - regular function that calculates S.F., V.P., and P.C. for 50000
%           vortices and then finds the number of vortices to get less 
%           than a 0.5% error (Bullet point 2)
%   2 - regular function that is used to conducts study to determine 
%           changes in 1) chord length, 2) angle of attack, and 3) free 
%           stream flow speed (Bullet points 1 and 3)
%   3 - loads saved data from the running of number 1 (Bullet point 1 & 2 data)

% angle of attack of 5 deg
alpha_a5 = deg2rad(5);
[x_a5,y_a5,Psi_a5,Phi_a5,~,~,bounds_a5,~] = Plot_Airfoil_Flow(c,alpha_a5,V_inf,P_inf,rho_inf,N,2);

% angle of attack of 10 deg
[x_a10,y_a10,Psi_a10,Phi_a10,~,~,bounds_a10,~] = Plot_Airfoil_Flow(c,alpha,V_inf,P_inf,rho_inf,N,2);

% angle of attack of 15 deg
alpha_a15 = deg2rad(15);
[x_a15,y_a15,Psi_a15,Phi_a15,~,~,bounds_a15,~] = Plot_Airfoil_Flow(c,alpha_a15,V_inf,P_inf,rho_inf,N,2);

% plot S.F.
figure
subplot(3,1,1)
contourf(x_a5,y_a5,Psi_a5,40) % stream function
hold on 
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_a5)
title(['Stream Lines for N = ' num2str(N) ' Vorticies at \alpha = ' num2str(5) '^{o}']);
subplot(3,1,2)
contourf(x_a10,y_a10,Psi_a10,40) % stream function
hold on
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_a10)
title(['Stream Lines for N = ' num2str(N) ' Vorticies at \alpha = ' num2str(10) '^{o}']);
ylabel('y')
subplot(3,1,3)
contourf(x_a15,y_a15,Psi_a15,40) % stream function
hold on
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_a15)
title(['Stream Lines for N = ' num2str(N) ' Vorticies at \alpha = ' num2str(15) '^{o}']);
xlabel('x')
hold off

% plot V.P.
figure
subplot(3,1,1)
contourf(x_a5,y_a5,Phi_a5,100) % velocity potential
hold on 
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_a5)
title(['Equipotential Lines for N = ' num2str(N) ' Vorticies at \alpha = ' num2str(5) '^{o}']);
subplot(3,1,2)
contourf(x_a10,y_a10,Phi_a10,100) % velocity potential
hold on
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_a10)
title(['Equipotential Lines for N = ' num2str(N) ' Vorticies at \alpha = ' num2str(10) '^{o}']);
ylabel('y')
subplot(3,1,3)
contourf(x_a15,y_a15,Phi_a15,100) % velocity potential
hold on
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_a15)
title(['Equipotential Lines for N = ' num2str(N) ' Vorticies at \alpha = ' num2str(15) '^{o}']);
xlabel('x')
hold off

fprintf('figures 7 & 8 display how the stream function and the velocity\npotential change as the angle of attack changes. When looking\nat the figures to view how the angle of attack affects the two\nit is seen that the angle of attack does have an effect\non the stream function or the velocity potential.\n\n');

%% Study on how free-stream flow speed changes S.F., V.P., and P.C.
% last input determines what runs in the function:
%   1 - regular function that calculates S.F., V.P., and P.C. for 50000
%           vortices and then finds the number of vortices to get less 
%           than a 0.5% error (Bullet point 2)
%   2 - regular function that is used to conducts study to determine 
%           changes in 1) chord length, 2) angle of attack, and 3) free 
%           stream flow speed (Bullet points 1 and 3)
%   3 - loads saved data from the running of number 1 (Bullet point 1 & 2 data)

% free-stream flow of 25
V_inf_25 = 10;
[x_v25,y_v25,Psi_v25,Phi_v25,~,~,bounds_v25,~] = Plot_Airfoil_Flow(c,alpha,V_inf_25,P_inf,rho_inf,N,2);

% free-stream flow of 50
[x_v50,y_v50,Psi_v50,Phi_v50,~,~,bounds_v50,~] = Plot_Airfoil_Flow(c,alpha,V_inf,P_inf,rho_inf,N,2);

% free-stream flow of 75
V_inf_75 = 90;
[x_v75,y_v75,Psi_v75,Phi_v75,~,~,bounds_v75,~] = Plot_Airfoil_Flow(c,alpha,V_inf_75,P_inf,rho_inf,N,2);

% plot S.F.
figure
subplot(3,1,1)
contourf(x_v25,y_v25,Psi_v25,40) % stream function
hold on 
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_v25)
title(['Stream Lines for N = ' num2str(N) ' Vorticies at V = ' num2str(V_inf_25) ' [m/s]']);
subplot(3,1,2)
contourf(x_v50,y_v50,Psi_v50,40) % stream function
hold on
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_v50)
title(['Stream Lines for N = ' num2str(N) ' Vorticies at V = ' num2str(V_inf) ' [m/s]']);
ylabel('y')
subplot(3,1,3)
contourf(x_v75,y_v75,Psi_v75,40) % stream function
hold on
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_v75)
title(['Stream Lines for N = ' num2str(N) ' Vorticies at V = ' num2str(V_inf_75) ' [m/s]']);
xlabel('x')
hold off

% plot V.P.
figure
subplot(3,1,1)
contourf(x_v25,y_v25,Phi_v25,100) % velocity potential
hold on 
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_v25)
title(['Equipotential Lines for N = ' num2str(N) ' Vorticies at V = ' num2str(V_inf_25) ' [m/s]']);
subplot(3,1,2)
contourf(x_v50,y_v50,Phi_v50,100) % velocity potential
hold on
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_v50)
title(['Equipotential Lines for N = ' num2str(N) ' Vorticies at V = ' num2str(V_inf) ' [m/s]']);
ylabel('y')
subplot(3,1,3)
contourf(x_v75,y_v75,Phi_v75,100) % velocity potential
hold on
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_v75)
title(['Equipotential Lines for N = ' num2str(N) ' Vorticies at V = ' num2str(V_inf_75) ' [m/s]']);
xlabel('x')
hold off

fprintf('figures 9 & 10 display how the stream function and the velocity\npotential change as the free stream velocity changes. When looking\nat the figures to view how the free stream velocity affects the two\nit is seen that the free stream velocity does have an effect\non the stream function or the velocity potential.\n\n');

toc