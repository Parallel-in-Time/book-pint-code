N=1000;                               % number of time steps
f=@(t,u) cos(t)*u;                    % RHS function
u0=1;                                 % initial value
[t,u]=ForwardEuler(f,[0 2*pi],u0,N);  % numerical solution