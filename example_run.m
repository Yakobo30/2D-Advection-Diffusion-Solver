% example_run.m
% Setup and visualization for the 2D Advection-Diffusion Solver

%% 1. Physical and Grid Parameters
vx = 1;
vy = 0.5;
Lx = 2;
Ly = 2;
Nx = 500;
Ny = 500;
dx = Lx/Nx;
dy = Ly/Ny;
dt = 0.001;
alpha = 0;
timeSteps = 1000;

% CHOOSE YOUR LIMITER HERE: 'MC' or 'minmod'
limiterType = 'MC'; 

%% 2. Initial Conditions Grid Setup
cx = zeros(1, Nx);
cy = zeros(1, Ny);
for i = 1:Nx
    cx(i) = (i - 0.5) * dx;
end
for j = 1:Ny
    cy(j) = (j - 0.5) * dy;
end

T_initial = zeros(Nx, Ny);
for i = 1:Nx
    for j = 1:Ny
        if (abs(cx(i) - 0.5) <= 0.1) && (abs(cy(j) - 0.25) <= 0.1)
            T_initial(i, j) = 1;
        end
    end
end

%% 3. Visualize Initial Conditions
figure;
surf(cx, cy, T_initial', 'EdgeColor', 'none');
title('Initial Pollution Distribution');
xlabel('X-Position'); ylabel('Y-Position'); zlabel('Concentration');

figure;
target_y1 = 0.25; % The physical Y-position you want to inspect
[~, y_idx1] = min(abs(cy - target_y1)); % Finds the closest index
plot(cx, T_initial(:, y_idx1));
title(sprintf('Initial 1D Cross-Section of Pollution at Y = %.2f', target_y1));
xlabel('X-Position');
ylabel('Concentration');

% Calculate the Initial Total Mass (sum of all cells * area of a cell)
initial_mass = sum(T_initial, 'all') * dx * dy;
fprintf('Initial Total Mass: %.6f\n', initial_mass);

%% 4. Run the Solver Tool
disp(['Running simulation using ', limiterType, ' limiter...']);
T_final = solveAdvectionDiffusion2D(T_initial, vx, vy, Lx, Ly, alpha, dt, timeSteps, limiterType);
disp('Simulation complete!');

%% 5. Visualize Final Results
figure;
surf(cx, cy, T_final', 'EdgeColor', 'none');
% Calculate the total physical time elapsed to display in the title
total_time = dt * timeSteps; 
title(sprintf('Final Concentration Distribution (t = %g s)', total_time));
xlabel('X-Position'); ylabel('Y-Position'); zlabel('Concentration');

figure;
target_y2 = 0.75; % The new physical Y-position the scalar moved to
[~, y_idx2] = min(abs(cy - target_y2)); % Finds the closest index
plot(cx, T_final(:, y_idx2));
title(sprintf('Final 1D Cross-Section at Y = %.2f (t = %g s)', target_y2, total_time));
xlabel('X-Position');
ylabel('Concentration');

% Calculate the Final Total Mass (sum of all cells * area of a cell)
final_mass = sum(T_final, 'all') * dx * dy;
fprintf('Final Total Mass: %.6f\n', final_mass);