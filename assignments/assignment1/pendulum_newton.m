function [t,theta,info] = pendulum_newton(K,theta1,theta2,m,tol,kmax)

if nargin < 4
    m = 80;
end
if nargin < 5
    tol = 1.e-10;
end
if nargin < 6
    kmax = 50;
end

h = 2*pi/m;
N = m-1;
t = (0:m)'*h;

theta = linspace(theta1, theta2, m+1)';
u = theta(2:m);

err = inf;
k = 0;

while err > tol && k < kmax
    R = zeros(N,1);
    mainDiag = zeros(N,1);
    offDiag = ones(N-1,1)/h^2;

    for i = 1:N
        if i == 1
            leftValue = theta1;
        else
            leftValue = u(i-1);
        end

        if i == N
            rightValue = theta2;
        else
            rightValue = u(i+1);
        end

        R(i) = (leftValue - 2*u(i) + rightValue)/h^2 + K*sin(u(i));
        mainDiag(i) = -2/h^2 + K*cos(u(i));
    end

    J = spdiags([[offDiag;0], mainDiag, [0;offDiag]], -1:1, N, N);
    delta = -J\R;
    u = u + delta;
    theta = [theta1; u; theta2];

    k = k + 1;
    err = norm(delta, inf);
end

info.iterations = k;
info.correctionNorm = err;
info.residualNorm = final_residual_norm(u,theta1,theta2,K,h);
info.converged = err <= tol;
info.h = h;
end

function residualNorm = final_residual_norm(u,theta1,theta2,K,h)
    N = length(u);
    R = zeros(N,1);

    for i = 1:N
        if i == 1
            leftValue = theta1;
        else
            leftValue = u(i-1);
        end

        if i == N
            rightValue = theta2;
        else
            rightValue = u(i+1);
        end

        R(i) = (leftValue - 2*u(i) + rightValue)/h^2 + K*sin(u(i));
    end

    residualNorm = norm(R, inf);
end
