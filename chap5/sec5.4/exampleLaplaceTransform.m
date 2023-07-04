J=100;                                   % J number of spatial unknowns
e=ones(J,1);                            
dx=1/(J+1);                              % finite difference Laplacian 
L=spdiags([e -2*e e],[-1 0 1],J,J)/dx^2; % Laplace matrix 
u0=20*e;                                 % initial temperature

N=100;                                   % number of quadrature points
t=0.1;                                   % time where we want the solution
ga=pi^2/2;                               % pi^2 smallest eigenvalue of L
tau=t/2;                                 % tau smaller than t
al=tau/(2*t);                            % parameter to reuse nodes
sigma=1/(al*t)*log(N./(1:N));            % quadrature nodes
w=exp(-ga*t)/(al*t*N)*((1:N)/N).^(1/al-1); % trapezoidal rule
w(end)=w(end)/2;                         
ws=2*w/3; ws(1:2:end)=ws(1:2:end)*2;     % transform to Simpson for N even
g=@(t,sigma,uh) 1/pi*imag((1+1i)*exp(-1i*sigma*t)*uh);
u=zeros(J,1);
for n=1:N
  s=ga+(1+1i)*sigma(n);
  uh=-(L+s*eye(J))\u0;
  u=u+w(n)*g(t,sigma(n),uh);
end
plot(dx:dx:1-dx,u);
xlabel('x')
axis([0 1 0 21]);
set(gca,'fontsize',18);
