function [t,u]=RIDC(lambda,tspan,u0,N,M,K)
% RIDC Solves the Dahlquist problem using RIDC.
%   [t,u]=RIDC(lambda,tspan,u0,N,M,K) solves the Dahlquist equation
%   with lambda on tspan=[tBeg, tEnd], starting from the initial value
%   u0. It uses the equivalent of N IDC time steps for IDC, with M
%   points per time steps (including the left and right boundary
%   point). Then it computes the prediction (K=0) using Backward
%   Euler, and corrects using K Backward Euler sweeps.

nPoints=N*(M-1)+1;
t=linspace(tspan(1),tspan(2),nPoints);      % Time grid for all intervals
dt=t(2)-t(1); u{1}(1)=u0;              
for n=1:nPoints-1                           % Initial guess
  u{1}(n+1)=1/(1-dt*lambda)*u{1}(n);
end
for k=1:K                                   % IDC for first M points
  [~, uI]=IDC(lambda,[t(1),t(M)],u0,1,M,k);
  u{k+1}(1:M)=uI;
end
for n=M:nPoints-1                           % Loop on remainding points
  for k=1:K                                 % Sweep iterations
    integrand=lambda*u{k}(n+2-M:n+1);       % Integral term
    weights=LagrangeWeights(M);
    integral=(M-1)*dt*sum(weights(end,:).*integrand);
    u{k+1}(n+1)=1/(1-dt*lambda)*( ...       % Backward Euler sweep
      u{k+1}(n)-dt*lambda*u{k}(n+1)+integral);
  end
end
u=u{K+1};                                   % Return only last iterate