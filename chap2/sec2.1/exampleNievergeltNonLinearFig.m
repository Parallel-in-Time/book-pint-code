% Settings for Matlab
close all; clear;
format long

f=@(t,u) cos(t)*exp(-u);                   % RHS of the ODE problem
u0=log(2); T=2*pi;                         % initial solution, final time
N=10;                                      % number of subintervals
nSteps=100;                                % fine steps per subintervals
[tPred,uPred]=ForwardEuler(f,[0,T],u0,N);  % approximate prediction
Mn=2; width=0.2;                          % Nievergelt's method parameters
solver=@(t0,t1,u0) UForwardEuler(f,t0,t1,u0,nSteps);
[U,uTraj]=Nievergelt(solver,T,u0,N,uPred,Mn,width);  % Nievergelt method

[tFine,uFine]=ForwardEuler(f,[0,T],u0,N*nSteps);

figure(1)
plot(tFine, uFine, DisplayName='Fine')
hold on; xlabel('Time'); ylabel('Solution');
plot(tPred, uPred, 'o', DisplayName='Prediction')
plot(tPred, U, '^', DisplayName='Nievergelt')
TT=linspace(0,T,N+1);
x = uTraj(1,1,:); x=x(:);
plot(linspace(TT(1),TT(2),nSteps+1), x, 'k--', DisplayName='Trajectories')
for n=2:N
    for m=1:Mn
        x = uTraj(n,m,:); x=x(:);
        plot(linspace(TT(n),TT(n+1),nSteps+1), x, 'k--', ...
            HandleVisibility='off')
    end
end
legend(Location='northeast'), set(gca,'fontsize', 12);
saveas(figure(1), '../Figures/NievergeltExampleNonLinear_1', 'epsc')

Mn=4; width=0.4;                          % Nievergelt's method parameters
solver=@(t0,t1,u0) UForwardEuler(f,t0,t1,u0,nSteps);
[U,uTraj]=Nievergelt(solver,T,u0,N,uPred,Mn,width);  % Nievergelt method

[tFine,uFine]=ForwardEuler(f,[0,T],u0,N*nSteps);

figure(2)
plot(tFine, uFine, DisplayName='Fine')
hold on; xlabel('Time'); ylabel('Solution');
plot(tPred, uPred, 'o', DisplayName='Prediction')
plot(tPred, U, '^', DisplayName='Nievergelt')
TT=linspace(0,T,N+1);
x = uTraj(1,1,:); x=x(:);
plot(linspace(TT(1),TT(2),nSteps+1), x, 'k--', DisplayName='Trajectories')
for n=2:N
    for m=1:Mn
        x = uTraj(n,m,:); x=x(:);
        plot(linspace(TT(n),TT(n+1),nSteps+1), x, 'k--', ...
            HandleVisibility='off')
    end
end
legend(Location='northeast'), set(gca,'fontsize', 12);
saveas(figure(2), '../Figures/NievergeltExampleNonLinear_2', 'epsc')
