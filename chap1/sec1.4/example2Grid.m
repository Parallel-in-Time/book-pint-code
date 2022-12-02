l=4;
J=2^l-1;
e=ones(J,1);
h=1/(J+1);
x=0:h:1;
A=1/h^2*spdiags([e -2*e e],[-1 0 1],J,J);
Jc=2^(l-1)-1;                             % coarse grid size
P=sparse(J,Jc);                           % prolongation by interpolation
for j=1:Jc
  P(2*j,j)=1; P(2*j-1,j)=0.5; P(2*j+1,j)=0.5;
end;
R=0.5*P';                                 % restriction
Ac=R*A*P;                                 % coarse matrix by Galerkin 
nu=2;                                     % number of smoothing steps    
rng('default'); u=rand(J,1);              % use same random initial guess
w=2/3;                                    % Jacobi damping parameter  
for n=1:3
  for i=1:nu
    u=u+w*h^2/2*A*u;                      % Jacobi damping step
  end;
  plot(x,[0;u;0],'-');                    % plot before coarse correction
  xlabel('x');ylabel(['error iter = ' num2str(i)])
  axis([0 1 -0.1 1])
  rc=R*(-A*u);                            % compute residual, f=0  
  u=u+P*(Ac\rc);                          % solve coarse problem
  hold on
  plot(x,[0;u;0],'-r');
  legend('before coarse','after coarse')
  hold off
  pause
end;
