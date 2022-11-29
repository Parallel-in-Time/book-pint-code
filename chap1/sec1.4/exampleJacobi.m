l=4;
J=2^l-1;                                     % number of grid points
e=ones(J,1);
h=1/(J+1);
x=0:h:1;
A=1/h^2*spdiags([e -2*e e],[-1 0 1],J,J);    % discrete Laplacian
w=1;                                         % damping parameter
rng('default'); u=rand(J,1);                 % always same random
for i=1:20
  plot(x,[0;u;0],'-');
  xlabel('x');ylabel(['error iter = ' num2str(i-1)])
  axis([0 1 -0.1 1])
  u=u+w*h^2/2*A*u;                           % u is the error, f=0
  pause
end;
