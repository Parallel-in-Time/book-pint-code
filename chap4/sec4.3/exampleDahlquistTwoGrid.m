l=6;                                    % can later become number of levels
N=2^l-1;                                % number of gridpoints in time
T=1; dt=T/N; t=0:dt:T;                  % time grid
la=-1;                                  % Dahlquist parameter
e=ones(N,1);
A=spdiags([-e (1-dt*la)*e],[-1 0],N,N); % BE time stepping matrix
u0=0;
al=0.5;                                 % Jacobi relaxation parameter
rng('default'); u=rand(N,1);            % random initial guess
Nc=2^(l-1)-1;                           % coarse grid size
P=sparse(N,Nc);                         % prolongation by interpolation
for j=1:Nc
  P(2*j,j)=1; P(2*j-1,j)=0.5; P(2*j+1,j)=0.5;
end;
R=0.5*P';                               % restriction
Ac=R*A*P;                               % coarse matrix by Galerkin
nu=4;                                   % number of presmoothing steps
al=0.5                                  % Jacobi relaxation parameter
rng('default'); u=rand(N,1);
for k=1:10
  err(k)=max(abs(u));
  for i=1:nu                            % presmoothing
    u=u-al/(1-dt*la)*A*u;
    plot(t,[u0;u],'-'); xlabel('t');ylabel('error')
    % axis([0 1 -0.1 1])
    % pause
  end;
  rc=R*(-A*u);                         % compute coarse correction
  u=u+P*(Ac\rc);
  hold on; plot(t,[u0;u],'-r'); hold off
  legend('before coarse','after coarse')
  pause
end;


