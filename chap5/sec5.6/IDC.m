function [t,u]=IDC(lambda,tspan,u0,N,M,K)
% IDC Solves the Dahlquist problem using Integral Deferred Correction.
%   [t,u]=IDC(lambda,tspan,u0,N,M,K) solves the Dahlquist equation
%   with lambda on tspan=[tBeg, tEnd], starting from initial value
%   u0. It uses N time steps for IDC, with M points per time step
%   (including the left and right boundary point). Then it computes the
%   prediction (K=0) using Backward Euler, and corrects using K
%   Backward Euler sweeps.

nPoints=N*(M-1)+1;
t=linspace(tspan(1),tspan(2),nPoints);      % Time grid for all intervals
dt=t(2)-t(1); uStep=u0;
for n=1:N
  u{n}{1}(1)=uStep;                         % Initial guess
  for m=1:M-1
    u{n}{1}(m+1)=1/(1-dt*lambda)*u{n}{1}(m);
  end
  for k=1:K                                 % Iteration sweeps 
    u{n}{k+1}(1)=uStep;
    integrand=lambda*u{n}{k};
    for m=1:M-1
      weights=LagrangeWeights(M);           % Integral term
      integral=(M-1)*dt*sum(weights(m, :).*integrand);
      u{n}{k+1}(m+1)=1/(1-dt*lambda)*( ...  % Sweep correction
      u{n}{k+1}(m)-dt*lambda*u{n}{k}(m+1)+integral);
    end
  end
  uStep=u{n}{K+1}(M);                       % Initial value for next step
end
uAll=u0;                                    % Combine each step solutions
for n=1:N
  uAll=[uAll,u{n}{K+1}(2:end)];
end
u=uAll;