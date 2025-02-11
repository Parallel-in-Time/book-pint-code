l=6;                                    % can later become number of levels
N=2^l-1;                                % number of gridpoints in time
T=1; dt=T/N; t=0:dt:T;                  % time grid
la=-1;                                  % Dahlquist parameter
e=ones(N,1);
A=spdiags([-e (1-dt*la)*e],[-1 0],N,N); % BE time stepping matrix
u0=0;
al=0.5;                                 % Jacobi relaxation parameter
rng('default'); u=rand(N,1);            % random initial guess
for i=1:20
  plot(t,[u0;u],'-');
  xlabel('t');ylabel(['error iter = ' num2str(i-1)])
  u=u-al/(1-dt*la)*A*u;           % damped Jacobi
  pause
end;
