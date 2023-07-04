function [t,u]=PIDC(lambda,tspan,u0,N,M,K)
%PIDC Solves the Dahlquist equation with Parallel IDC.
%   [t,u]=PIDC(lambda,tspan,u0,N,M,K) solves the Dahlquist equation
%   with lambda on tspan=[tBeg, tEnd], starting from the initial value
%   u0. It uses N parallel IDC time steps, with M points per time
%   step (including the left and right boundary point). Then it computes
%   the prediction (K=0) using Backward Euler, and corrects it using K
%   Backward Euler sweeps.

nPoints=N*(M-1)+1;
t=linspace(tspan(1),tspan(2),nPoints);      % Time grid for all intervals
dt=t(2)-t(1); uStep=ones(1,K+1)*u0;         % Initial guess
for n=1:N                                   
  u{n}{1}(1)=uStep(1);
  for m=1:M-1
    u{n}{1}(m+1)=1/(1-dt*lambda)*u{n}{1}(m);
  end
  for k=1:K                                 % Sweep iterations
    u{n}{k+1}(1)=uStep(k+1);
    integrand=lambda*u{n}{k};
    for m=1:M-1
      weights=LagrangeWeights(M);           % Integral term
      integral=(M-1)*dt*sum(weights(m, :).*integrand);
      u{n}{k+1}(m+1)=1/(1-dt*lambda)*( ...  % Sweep correction
        u{n}{k+1}(m)-dt*lambda*u{n}{k}(m+1)+integral);
    end
  end
  for k=1:K+1
    uStep(k)=u{n}{k}(M);                    % Initial values for next step
  end
end
uAll=u0;                                    % Combine each step
for n=1:N
  uAll=[uAll,u{n}{K+1}(2:end)];
end
u=uAll;