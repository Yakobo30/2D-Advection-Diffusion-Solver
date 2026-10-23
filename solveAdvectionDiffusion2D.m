function T = solveAdvectionDiffusion2D(T_initial, vx, vy, Lx, Ly, alpha, dt, timeSteps, limiterType)
    % solveAdvectionDiffusion2D Solves the 2D Advection-Diffusion equation.
    %
    % This solver utilizes a cell-centered Finite Volume Method (FVM) with 
    % Forward Euler time integration and a Second-Order Upwind (SOU) scheme.
    % 
    % Inputs:
    %   T_initial   - [Nx x Ny matrix] Initial concentration/temperature field
    %   vx, vy      - [scalars] Advection velocities in x and y directions
    %   Lx, Ly      - [scalars] Physical domain dimensions
    %   alpha       - [scalar] Diffusion coefficient
    %   dt          - [scalar] Time step size
    %   timeSteps   - [integer] Number of temporal iterations
    %   limiterType - [string] Flux limiter choice: 'MC' or 'minmod'
    %
    % Outputs:
    %   T           - [Nx x Ny matrix] Final concentration field

    % Determine grid size from the initial condition matrix
    [Nx, Ny] = size(T_initial);
    dx = Lx / Nx;
    dy = Ly / Ny;
    
    T = T_initial;
    
    % Time loop
    for k = 1:timeSteps
        % Forward Euler
        dTdt = find_dTdt(T, vx, vy, alpha, dx, dy, Nx, Ny, limiterType);
        T = T + dTdt * dt;
    end
end

function dTdt = find_dTdt(Tcurrent, vx, vy, alpha, dx, dy, Nx, Ny, limiterType)
    flux_x = zeros(Nx+1, Ny);
    flux_y = zeros(Nx, Ny+1);
    dTdt = zeros(Nx, Ny);
    dif_x = zeros(Nx+1, Ny);
    adv_x = zeros(Nx+1, Ny);
    dif_y = zeros(Nx, Ny+1);
    adv_y = zeros(Nx, Ny+1);
    
    % Flux on the Horizontal Plane
    for xSurface = 1:Nx+1
        for y = 1:Ny
            if (xSurface == 1) || (xSurface == Nx+1)
                flux_x(xSurface,y) = 0;
            else
                % Diffusive Flux (using scalar temporary variable for memory efficiency)
                dif_x = -alpha*(Tcurrent(xSurface,y)-Tcurrent(xSurface-1,y))/dx;             
                
                % Advective Flux
                if vx > 0
                    if (xSurface==2)
                        adv_x = vx*Tcurrent(xSurface-1,y);
                    else
                        % Using Limiter (Second-Order Upwind Scheme)
                        fw_dTdx = (Tcurrent(xSurface,y)-Tcurrent(xSurface-1,y))/dx;
                        bw_dTdx = (Tcurrent(xSurface-1,y)-Tcurrent(xSurface-2,y))/dx;
                        
                        switch limiterType
                            case 'MC'
                                dTdx = mclimiter(fw_dTdx, bw_dTdx);
                            case 'minmod'
                                dTdx = minmod(fw_dTdx, bw_dTdx);
                            otherwise
                                error('Unknown limiter type. Choose ''MC'' or ''minmod''.');
                        end
                        
                        adv_x = vx*(Tcurrent(xSurface-1,y)+0.5*dx*dTdx);
                    end
                else
                    if (xSurface == Nx)
                        adv_x = vx*Tcurrent(xSurface,y);
                    else
                        % Using Limiter (Second-Order Upwind Scheme)
                        fw_dTdx = (Tcurrent(xSurface-1,y)-Tcurrent(xSurface,y))/dx;
                        bw_dTdx = (Tcurrent(xSurface,y)-Tcurrent(xSurface+1,y))/dx;
                        
                        switch limiterType
                            case 'MC'
                                dTdx = mclimiter(fw_dTdx, bw_dTdx);
                            case 'minmod'
                                dTdx = minmod(fw_dTdx, bw_dTdx);
                            otherwise
                                error('Unknown limiter type. Choose ''MC'' or ''minmod''.');
                        end
                        
                        adv_x = vx*(Tcurrent(xSurface,y)+0.5*dx*dTdx);
                    end
                end
                
                % Combine fluxes
                flux_x(xSurface,y) = dif_x + adv_x;
            end
        end
    end
    
    % Flux on the Vertical Plane
    for ySurface = 1:Ny+1
        for x = 1:Nx
            if (ySurface == 1) || (ySurface == Ny+1)
                flux_y(x,ySurface) = 0;
            else
                % Diffusive Flux (using scalar temporary variable)
                dif_y = -alpha*(Tcurrent(x,ySurface)-Tcurrent(x,ySurface-1))/dy;
                
                % Advective Flux
                if vy > 0
                    if (ySurface == 2)
                        adv_y = vy*Tcurrent(x,ySurface-1);
                    else
                        % Using Limiter (Second-Order Upwind Scheme)
                        fw_dTdy = (Tcurrent(x,ySurface)-Tcurrent(x,ySurface-1))/dy;
                        bw_dTdy = (Tcurrent(x,ySurface-1)-Tcurrent(x,ySurface-2))/dy;
                        
                        switch limiterType
                            case 'MC'
                                dTdy = mclimiter(fw_dTdy, bw_dTdy);
                            case 'minmod'
                                dTdy = minmod(fw_dTdy, bw_dTdy);
                            otherwise
                                error('Unknown limiter type. Choose ''MC'' or ''minmod''.');
                        end
                        
                        adv_y = vy*(Tcurrent(x,ySurface-1)+0.5*dy*dTdy);
                    end
                else
                    if (ySurface == Ny)
                        adv_y = vy*Tcurrent(x,ySurface);
                    else
                        % Using Limiter (Second-Order Upwind Scheme)
                        fw_dTdy = (Tcurrent(x,ySurface-1)-Tcurrent(x,ySurface))/dy;
                        bw_dTdy = (Tcurrent(x,ySurface)-Tcurrent(x,ySurface+1))/dy;
                        
                        switch limiterType
                            case 'MC'
                                dTdy = mclimiter(fw_dTdy, bw_dTdy);
                            case 'minmod'
                                dTdy = minmod(fw_dTdy, bw_dTdy);
                            otherwise
                                error('Unknown limiter type. Choose ''MC'' or ''minmod''.');
                        end
                        
                        adv_y = vy*(Tcurrent(x,ySurface)+0.5*dy*dTdy);
                    end
                end
                
                % Combine fluxes
                flux_y(x,ySurface) = dif_y + adv_y;
            end
        end
    end
    
    % Calculate dT/dt for each cell
    for x = 1:Nx
        for y = 1:Ny
            dTdt(x,y) = -(flux_x(x+1,y)-flux_x(x,y))/dx - (flux_y(x,y+1)-flux_y(x,y))/dy;
        end
    end
end

function limited_slope = minmod(fw, bw)
    % If slopes are opposite directions, drop to First-Order
    if (fw * bw <= 0)
        limited_slope = 0;
    else
        % Otherwise, pick the absolute minimum
        if abs(fw) < abs(bw)
            limited_slope = fw;
        else
            limited_slope = bw;
        end
    end
end

function limited_slope = mclimiter(fw, bw)
    % If slopes are opposite directions, drop to First-Order
    if (fw * bw <= 0)
        limited_slope = 0;
    else
        % Using MClimiter comparing central difference with minmod 
        % alpha multiplier (In this case let alpha = 2)
        centralDiff = 0.5 * (fw + bw);
        fw2 = 2 * fw;
        bw2 = 2 * bw;
        limited_slope = sign(centralDiff) * min([abs(centralDiff), abs(fw2), abs(bw2)]);
    end
end