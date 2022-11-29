% Settings for Matlab
close all; clear;
set( 0, 'defaultAxesTickLabelInterpreter', 'latex' );
set( 0, 'defaultLegendInterpreter',        'latex' );
set( 0, 'defaultTextInterpreter',          'latex' );
set( 0, 'defaultAxesFontSize', 12 );
format long

L = pi;
c = 1;
J = 15;
tSpan = 16*pi;

x = linspace(0, L, J+2);

u0 = sin(2*x);
ut0 = 0;

nStep = 300;
[t1, u1] = WaveEquation(u0, ut0, c, L, tSpan, nStep);

nStep = 600;
[t2, u2] = WaveEquation(u0, ut0, c, L, tSpan, nStep);

nStep = 1200;
[t3, u3] = WaveEquation(u0, ut0, c, L, tSpan, nStep);

figure('Name','u(pi/4,T)','NumberTitle','off')
plot(t1, u1(5, :), '-^', 'DisplayName', 'nStep=300', 'MarkerIndices',1:3:length(t1));
hold on;
plot(t2, u2(5, :), 's-', 'DisplayName', 'nStep=600', 'MarkerIndices',1:6:length(t2));
plot(t3, u3(5, :), 'o-', 'DisplayName', 'nStep=1200', 'MarkerIndices',1:12:length(t3));
uTh = sin(2*pi/4)*cos(2*c*t2);
plot(t2, uTh, '--', 'DisplayName', 'Analytical')
grid on
legend('Location', 'east')
xlim([14*pi, 53]);
xlabel('$t$')
ylabel('$u(\pi/4,t)$')