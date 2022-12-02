function U=MultipleShooting(P,T,u0,N,K,uPred)
% MULTIPLESHOOTING implementation of the multiple shooting method
%   U=MultipleShooting(P,T,u0,N,K,uPred); applies the Multiple
%   shooting algorithm on [0,T] with initial condition u0 at t=0 using 
%   N equidistant coarse time points doing K iterations. It uses a 
%   propagator P(t0,t1,ut0) that returns the solution and the Jacobian 
%   at time t1, and an initial guess uPred containing N+1 starting
%   values for the algorithm.
dT=T/N; TT=0:dT:T;                               % coarse time mesh
U{1}(1,:)=u0;
for n=1:N                                        % initial guess 
    U{1}(n+1,:)=uPred(n+1,:);
end
for k=1:K                                        % iterations
  U{k+1}(1,:)=u0;
  for n=1:N
    % Propagation of the fine solution and the Jacobian
    [u, V] = P(TT(n), TT(n+1), U{k}(n,:));             
    % Newton correction
    U{k+1}(n+1,:) = u + (V*(U{k+1}(n,:) - U{k}(n,:))')';  
  end
end