% Set Lorenz system parameters and define functions
sigma=10; rho=28; beta=8/3;

f=@(u) [sigma*(u(2)-u(1)), u(1)*(rho-u(3))-u(2), u(1)*u(2) - beta*u(3)]';

fA=@(u) [[-sigma, sigma, 0];
         [rho-u(3), -1, -u(1)];
         [u(2), u(1), -beta]];

% Apply Newton methods
x0Goal = [(beta*(rho-1))^0.5, (beta*(rho-1))^0.5, rho-1]';
x0 = x0Goal - 5;
[eps, x] = Newton(f, fA, x0, 10);

% Plot error
semilogy(eps, '-o');
xlabel('$N_{iter}$');
ylabel('Error');
grid on;
disp(x0Goal-x);