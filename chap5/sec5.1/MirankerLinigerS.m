function [t,y]=MirankerLinigerS(f,tspan,y0,N);
% MirankerLiniger1 solves ODEs using a predictor corrector method
%   [t,y]=MirankerLinigerS(f,tspan,y0,N) solves dy/dt=f(t,y) with
%   initial value y0 on the time interval tspan doing n steps of
%   a sequential predictor corrector method from Miranker and Liniger

dt=(tspan(2)-tspan(1))/N;
t=(tspan(1):dt:tspan(2))';   
y(:,1)=y0(:);                               % colon to make column vector
y(:,2)=y(:,1)+dt*f(t(1),y(:,1));            % start with an Euler step
for n=2:N-1,
  yp=y(:,n)+dt/2*(3*f(t(n),y(:,n))-f(t(n-1),y(:,n-1))); % sequential predictor
  y(:,n+1)=y(:,n)+dt/2*(f(t(n+1),yp)+f(t(n),y(:,n)));   % corrector method
end;