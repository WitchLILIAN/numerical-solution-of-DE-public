% Problem 2: Lax-Wendroff for u_t + u_x = 0, deltat=1
clear; clc;

m = 30;
h = 1/m;
k = h; % deltat=h
N = round(1/k); % N=30 time steps
x = linspace(0,1,m+1);

% Initial condition u(x,0)
U0 = zeros(1,m+1);
U0(m/2+1:m+1) = 1; % for x>=0.5

% Analytical solution at t=1
u_exact = sin(1-x);
u_exact(end) = 0; % u(1,1)=0

% Three NBC methods
methods = {'First-order extrap', 'Second-order extrap', 'Upwind'};
colors = {'b--o', 'r-^', 'g-.s'};

figure('Name', 'Comparison of NBC Methods at t=1');
hold on;
plot(x, u_exact, 'k-', 'LineWidth', 2, 'DisplayName', 'Analytic');

max_errors = zeros(3,1);
err = zeros(3,m+1);

for method = 1:3
    U = U0;
    for n = 1:N
        t_curr = (n-1)*k; 
        
        % Store the old time t=k-1 value of U(m) for the Upwind NBC
        U_old_m = U(m);
        
        % Interior update: Lax-Wendroff is reduced to U(i) = U(i-1)
        % Update the new time t=k U(m) from right to left to avoid overwriting
        for i = m+1:-1:2
            U(i) = U(i-1);
        end
        
        % Left boundary
        U(1) = sin(t_curr + k);
        
        % Right boundary (NBC)
        if method == 1 % First-order extrapolation
            U(m+1) = U(m); % U(m+1,k)=U(m,k)                      
        elseif method == 2 % Second-order extrapolation
            U(m+1) = -U(m-1) + 2*U(m) ; % U(m+1,k)=-U(m-1,k)+2*U(m,k)
        else % Upwind scheme
            U(m+1) = U_old_m; % U(m+1,k)=U(m,k-1)
        end
    end
    
    % Plot this numerical solution 
    plot(x, U, colors{method}, 'LineWidth', 1.5, 'DisplayName', methods{method});
    
    % Compute errors
    err(method,:) = abs(U - u_exact);
    max_errors(method) = max(err(method,:));
end

figure(1);
xlabel('x'); ylabel('u');
title('Comparison of Three NBC Methods (Lax-Wendroff, \Delta t = h)');
legend('Location', 'best');
grid on;
xlim([0 1]);

figure(2);hold on;
for method = 1:3
    plot(x, err(method,:), colors{method}, 'LineWidth', 1.5,...
        'DisplayName', methods{method});
end
xlabel('x'); ylabel('Absolute Error');
title('Error Distribution at t=1');
legend('Location', 'best');
grid on;
xlim([0 1]);

% Print errors
fprintf('Error Comparison at t=1\n');
fprintf('%-25s: Max Error = %.15f\n', methods{1}, max_errors(1));
fprintf('%-25s: Max Error = %.15f\n', methods{2}, max_errors(2));
fprintf('%-25s: Max Error = %.15f\n', methods{3}, max_errors(3));