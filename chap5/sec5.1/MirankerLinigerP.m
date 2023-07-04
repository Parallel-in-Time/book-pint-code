function [t,y]=MirankerLinigerP(f,tspan,y0,N);
% MirankerLinigerP solves ODEs using a predictor corrector method
%   [t,y]=MirankerLinigerP(f,tspan,y0,N) solves dy/dt=f(t,y) with
%   initial value y0 on the time interval tspan doing n steps of
%   a parallel predictor corrector method from Miranker and Liniger

dt=(tspan(2)-tspan(1))/N;
t=(tspan(1):dt:tspan(2))';   
y(:,1)=y0(:);                           % colon to make column vector
y(:,2)=y(:,1)+dt*f(t(1),y(:,1));        % start with an Euler step
yp=y(:,2);
for n=3:N,
  y(:,n)=y(:,n-1)+dt/2*(f(t(n),yp)+f(t(n-1),y(:,n-1))); % parallel predictor
  yp=y(:,n-1)+2*dt*f(t(n),yp);                          % corrector method
end;