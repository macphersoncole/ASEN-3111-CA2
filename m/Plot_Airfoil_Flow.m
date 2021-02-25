function [x,y,Psi,Phi,P,levels,bounds,N_final,err_vec] = Plot_Airfoil_Flow(c,alpha,V_inf,P_inf,rho_inf,N,step)
% Function to calculate and plot the stream function, velocity potential,
% and the pressure contour for the given inputs
%
% Inputs:   c       - chord length of the airfoil desired to model
%           alpha   - angle of attack of the airfoil
%           V_inf   - free stream velocity of the uniform flow
%           P_inf   - free stream pressure of the uniform flow
%           rho_inf - free stream density of the uniform flow
%           N       - number of vortecies to model the airfoil
%
% Outputs:  Psi     - stream function of the uniform flow and vortices
%           Phi     - velocity potential of the uniform flow and vortices
%           P       - pressure contour of the unifrom flow and vortices

    if step == 1
%       1 - regular function that calculates S.F., V.P., and P.C. for 50000
%           vortices and then finds the number of vortices to get less 
%           than a 0.5% error (Bullet point 2)

        %% Domain Declaration
        x_min = -c/2;
        x_max = c+c/2;
        y_min = -c/2;
        y_max = c/2;
        bounds = [x_min x_max y_min y_max];

        %% Define Number of Grid Points
        n_x = 100; % steps in the x direction
        n_y = 100; % steps in the y direction

        %% Define Mesh
        N_2 = 50000;
        [x,y] = meshgrid(linspace(x_min,x_max,n_x),linspace(y_min,y_max,n_y));

        %% Define a function that calculates radius
        r =@(x1) sqrt((x-x1).^2 + (y).^2);

        %% Define the length and the points at which the vortcies exist
        delta_x = c/N_2;
        x_vortex = linspace(delta_x/2,c-delta_x,N_2);

        %% Calculate vortex strength and circulation
        gamma = 2*alpha*V_inf*sqrt((1-(x_vortex/c))./(x_vortex/c)); % vortex strength
        Gamma = gamma*delta_x; % circulation

        %% Calculate psi for uniform flow (Eq. 3.55)
        Psi_uniform = V_inf*(y*cos(alpha)-x*sin(alpha));

        %% Calculate phi for uniform flow (Eq. 3.53)
        Phi_uniform = V_inf*(x*cos(alpha)-y*sin(alpha));

        %% Calculate psi for vortex (Eq. 3.114)
        Psi_vortex = 0;
        for i = 1:N_2
            % calculate the stream function of each vortices and add them
            % together to get the total stream function of the system
            Psi_vortex = Psi_vortex + (Gamma(i)*log(r(x_vortex(i))))/(2*pi);
        end

        %% Calulate phi for vortex (Eq. 3.112)
        Phi_vortex = 0;
        for i = 1:N_2
            % calculate the velocity potential of each vortices and add them
            % together to get the total velocity potential of the system
            theta = atan2(-y,-x+x_vortex(i));
            Phi_vortex = Phi_vortex + -(Gamma(i)*theta)/(2*pi);
        end

        %% Add all the stream functions
        Psi_N = Psi_vortex + Psi_uniform;

        %% Add all the velocity potentials
        Phi_N = Phi_vortex + Phi_uniform;

        %% Calculate P
        q_inf = (1/2)*rho_inf*V_inf^2; % dynamic pressure
        C_p = 1-(gradient(Phi_N,(x_max-x_min)/n_x,(y_max-y_min)/n_y)./V_inf).^2;
        P_N = P_inf + C_p*q_inf; % total pressure

        %% Initialize error
        err = 1; 
        j = 1;
        err_vec = zeros(length(N:1000:24000),1);

        %% While loop to find the number of vortices needed for error less than 0.01
        while err > 0.005

            %% Define Mesh
            [x,y] = meshgrid(linspace(x_min,x_max,n_x),linspace(y_min,y_max,n_y));

            %% Define a function that calculates radius
            r =@(x1) sqrt((x-x1).^2 + (y).^2);

            %% Define the length and the points at which the vortcies exist
            delta_x = c/N;
            x_vortex = linspace(delta_x/2,c-delta_x,N);

            %% Calculate vortex strength and circulation
            gamma = 2*alpha*V_inf*sqrt((1-(x_vortex/c))./(x_vortex/c)); % vortex strength
            Gamma = gamma*delta_x; % circulation

            %% Calculate psi for uniform stream (Eq. 3.55)
            Psi_uniform = V_inf*(y*cos(alpha)-x*sin(alpha));

            %% Calculate phi for uniform stream (Eq. 3.53)
            Phi_uniform = V_inf*(x*cos(alpha)-y*sin(alpha));

            %% Calculate psi for vortex (Eq. 3.114)
            Psi_vortex = 0;
            for i = 1:N
                % calculate the stream function of each vortices and add them
                % together to get the total stream function of the system
                Psi_vortex = Psi_vortex + (Gamma(i)*log(r(x_vortex(i))))/(2*pi);
            end

            %% Calulate phi for vortex (Eq. 3.112)
            Phi_vortex = 0;
            for i = 1:N
                % calculate the velocity potential of each vortices and add them
                % together to get the total velocity potential of the system
                theta = atan2(-y,-x+x_vortex(i));
                Phi_vortex = Phi_vortex + -(Gamma(i)*theta)/(2*pi);
            end

            %% Add all the stream functions
            Psi = Psi_vortex + Psi_uniform;

            %% Add all the velocity potentials
            Phi = Phi_vortex + Phi_uniform;

            %% Calculate P
            q_inf = (1/2)*rho_inf*V_inf^2; % dynamic pressure
            C_p = 1-(gradient(Phi,(x_max-x_min)/n_x,(y_max-y_min)/n_y)./V_inf).^2;
            P = P_inf + C_p*q_inf; % total pressure

            %% Determine color levels for stream function contours
            levmin = Psi(1,n_x); % defines the color levels -> trial and error to find a good representation
            levmax = Psi(n_y,n_x/50);
            levels = linspace(levmin,levmax,100)';

            %% Find the error associated with the 
            Psi_err = mean(abs((Psi_N-Psi)./Psi_N),'all');
            Phi_err = mean(abs((Phi_N-Phi)./Phi_N),'all');
            P_err = mean(abs((P_N-P)./P_N),'all');

            err = max([Psi_err Phi_err P_err]);
            err_vec(j) = err; 
            j = j + 1;

            if err > 0.005
                N = N + 1000;
            end
            % fprintf('%0.5f \n %i \n',err,N);

        end

        N_final = N;
        
    elseif step == 2
%       2 - regular function that is used to conducts study to determine 
%           changes in 1) chord length, 2) angle of attack, and 3) free 
%           stream flow speed (Bullet points 1 and 3)
        
        %% Domain Declaration
        x_min = -c/2;
        x_max = c+c/2;
        y_min = -c/2;
        y_max = c/2;
        bounds = [x_min x_max y_min y_max];
        N_final = N;

        %% Define Number of Grid Points
        n_x = 100; % steps in the x direction
        n_y = 100; % steps in the y direction

        %% Define Mesh
        [x,y] = meshgrid(linspace(x_min,x_max,n_x),linspace(y_min,y_max,n_y));

        %% Define a function that calculates radius
        r =@(x1) sqrt((x-x1).^2 + (y).^2);

        %% Define the length and the points at which the vortcies exist
        delta_x = c/N;
        x_vortex = linspace(delta_x/2,c-delta_x,N);

        %% Calculate vortex strength and circulation
        gamma = 2*alpha*V_inf*sqrt((1-(x_vortex/c))./(x_vortex/c)); % vortex strength
        Gamma = gamma*delta_x; % circulation

        %% Calculate psi for uniform flow (Eq. 3.55)
        Psi_uniform = V_inf*(y*cos(alpha)-x*sin(alpha));

        %% Calculate phi for uniform flow (Eq. 3.53)
        Phi_uniform = V_inf*(x*cos(alpha)-y*sin(alpha));

        %% Calculate psi for vortex (Eq. 3.114)
        Psi_vortex = 0;
        for i = 1:N
            % calculate the stream function of each vortices and add them
            % together to get the total stream function of the system
            Psi_vortex = Psi_vortex + (Gamma(i)*log(r(x_vortex(i))))/(2*pi);
        end

        %% Calulate phi for vortex (Eq. 3.112)
        Phi_vortex = 0;
        for i = 1:N
            % calculate the velocity potential of each vortices and add them
            % together to get the total velocity potential of the system
            theta = atan2(-y,-x+x_vortex(i));
            Phi_vortex = Phi_vortex + -(Gamma(i)*theta)/(2*pi);
        end

        %% Add all the stream functions
        Psi = Psi_vortex + Psi_uniform;

        %% Add all the velocity potentials
        Phi = Phi_vortex + Phi_uniform;

        %% Calculate P
        q_inf = (1/2)*rho_inf*V_inf^2; % dynamic pressure
        C_p = 1-(gradient(Phi,(x_max-x_min)/n_x,(y_max-y_min)/n_y)./V_inf).^2;
        P = P_inf + C_p*q_inf; % total pressure
        
        %% Determine color levels for stream function contours
        levmin = Psi(1,n_x); % defines the color levels -> trial and error to find a good representation
        levmax = Psi(n_y,n_x/2);
        levels = linspace(levmin,levmax,100)';
        err_vec = 1;
        
    elseif step == 3
%       3 - loads saved data from the running of number 1 (Bullet point 1 & 2 data)
        
        %% Load data from bullet 1 and 2
        load('CA2_Data_MacPherson_Cole_1');
        N_final = N2;
        err_vec = err;
        
    end
    
end