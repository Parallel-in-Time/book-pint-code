% Settings for Matlab
close all; clear;
set( 0, 'defaultAxesTickLabelInterpreter', 'latex' );
set( 0, 'defaultLegendInterpreter',        'latex' );
set( 0, 'defaultTextInterpreter',          'latex' );
set( 0, 'defaultAxesFontSize', 12 );
format long 

% Parameters
nF = 5;

n = 0:nF-1;
x = linspace(0, 1, 500);
global nMesh xMesh;
[xMesh, nMesh] = meshgrid(x, n);


% Plot solution
plot(x, u(1/8), 'DisplayName', '$T=1/8$');
hold on;
plot(x, u(1/4), 'DisplayName', '$T=1/4$')
plot(x, u(1/2), 'DisplayName', '$T=1/2$')
legend(); xlabel('$x$'); ylabel('$u(x,T)$');
grid on;

% Function to compute the analytical solution at time t
function sol=u(t)
    global nMesh xMesh; 
    s = (80/pi)./(2*nMesh+1) .* exp(-(2*nMesh+1).^2 * pi^2 * t) ...
        .* sin((2*nMesh+1).*pi .*xMesh);
    sol = sum(s, 1);
end