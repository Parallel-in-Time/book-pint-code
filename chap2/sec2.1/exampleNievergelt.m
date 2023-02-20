f=@(t,u) cos(t)*u;                         % RHS of the ODE problem
u0=1; T=2*pi;                              % initial value, final time
N=10;                                      % number of subintervals
nSteps=100;                                % fine steps per subinterval
[tPred,uPred]=ForwardEuler(f,[0,T],u0,N);  % approximate prediction
Mn=2; width=0.75;                          % algorithm parameters
solver=@(t0,t1,u0) UForwardEuler(f,t0,t1,u0,nSteps);
[U,uTraj]=Nievergelt(solver,T,u0,N,uPred,Mn,width);  