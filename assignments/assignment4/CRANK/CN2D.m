% Crank-Nicolson Solver for 2D Heat Equation (m=n=32)

clear; close all; 

% Grid Parameters 
m = 32;
n = 32;
a = 0; b = 1;   % Domain in x
c = 0; d = 1;   % Domain in y
h = (b-a)/m;    % Spatial step (hx = hy = h)
k = h;          % Time step (dt = h)
T_final = 0.5;  % Final time
N_steps = round(T_final / k);  % Number of time steps

x = linspace(a, b, m+1);
y = linspace(c, d, n+1);
[X, Y] = meshgrid(x, y);   % Generate 2D mesh

% Initial Condition (t=0)
U = exp(0).*sin(X).*cos(Y);      % u(x,y,0)

tic;
% Build Coefficient Matrix A 
% In 2D, map grid point (i,j) to a 1D index idx = (j-1)*(m+1) + i
N_unknowns = (m+1)*(n+1);
A = sparse(N_unknowns, N_unknowns);

r = k / (2 * h^2);   

for j = 1:n+1
    for i = 1:m+1
        idx = (j-1)*(m+1) + i;
        
        if (i == 1) || (i == m+1) || (j == 1) || (j == n+1)
            % Dirichlet Boundary Point: Matrix row set to 1
            A(idx, idx) = 1; 
        else
            % Interior Point: Apply Crank-Nicolson stencil
            % (1+4r) * U(i,j) - r*U(i-1,j) - r*U(i+1,j) - r*U(i,j-1) - r*U(i,j+1) = RHS
            A(idx, idx) = 1 + 4*r;
            A(idx, idx-1) = -r;     % Left neighbor (i-1, j)
            A(idx, idx+1) = -r;     % Right neighbor (i+1, j)
            A(idx, idx-(m+1)) = -r; % Bottom neighbor (i, j-1)
            A(idx, idx+(m+1)) = -r; % Top neighbor (i, j+1)
        end
    end
end

% Time Iteration 
% Initialize error and solution storage for visualization later
U_final = U; 

for step = 1:N_steps
    t_curr = (step-1) * k;    % Current time
    t_next = step * k;        % Next time
    
    b = zeros(N_unknowns, 1);
    
    % Build Right-Hand Side Vector b
    for j = 1:n+1
        for i = 1:m+1
            idx = (j-1)*(m+1) + i;
            
            if (i == 1)
                % Left Boundary (x=0): u = uexact(t, 0, y)
                b(idx) = exp(-t_next)*sin(x(i))*cos(y(j)); 
            elseif (i == m+1)
                % Right Boundary (x=1): u = uexact(t, 1, y)
                b(idx) = exp(-t_next)*sin(x(i))*cos(y(j));
            elseif (j == 1)
                % Bottom Boundary (y=0): u = uexact(t, x, 0)
                b(idx) = exp(-t_next)*sin(x(i))*cos(y(j));
            elseif (j == n+1)
                % Top Boundary (y=1): u = uexact(t, x, 1)
                b(idx) = exp(-t_next)*sin(x(i))*cos(y(j));
            else
                % Interior Point: Crank-Nicolson RHS
                % RHS = U(i,j) + r*(U(i-1,j)+U(i+1,j)+U(i,j-1)+U(i,j+1)-4*U(i,j)) + k/2*(f_new + f_old)
                % meshgrid gives U(row,column)=U(y_j,x_i), so use U(j,i)
                rhs_old = U(j,i) + r * (U(j,i-1) + U(j,i+1) + U(j-1,i) + U(j+1,i) - 4*U(j,i));
                f_old = exp(-t_curr)*sin(x(i))*cos(y(j));
                f_new = exp(-t_next)*sin(x(i))*cos(y(j));
                b(idx) = rhs_old + (k/2) * (f_new + f_old);
            end
        end
    end
    
    % Solve the linear system A * U_new = b
    U_vec = A \ b;
    
    % Reshape the vector back to a 2D grid
    U = reshape(U_vec, m+1, n+1)';
    U_final = U;   % Store final solution for plotting
end
cpu_time = toc;

% Error Computation and Plotting
% Exact solution at final time
U_exact = exp(-T_final).*sin(X).*cos(Y);

% Calculate errors
error = abs(U_final - U_exact);
max_error = max(error(:));

fprintf('Max error of Crank-Nicolson in 2D (m=n=%d):',m);
fprintf(' = %.15f\n', max_error);
fprintf('CPU time of CN with m=n=%d: %.6f seconds\n',m,cpu_time);

% plot Numerical Solution at t=0.5
figure(1);
hold on;
subplot(1,2,1);
mesh(X, Y, U_final);
xlabel('x'); ylabel('y'); zlabel('u(x,y,t)');
title(sprintf('Numerical solution by CN Method with m=n=%d grid intervals', m));

% Plot Error Distribution at t=0.5
subplot(1,2,2);
mesh(X, Y, abs(U_final - U_exact));
xlabel('x'); ylabel('y'); zlabel('Error');
title(sprintf('Error Distribution of CN Method with m=n=%d grid intervals', m));
