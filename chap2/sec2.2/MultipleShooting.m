function [t,U] = MultipleShooting(prop,T,u0,N,K,uPred)
%MULTIPLESHOOTING Implementation of MultipleShooting
% ...
t=linspace(0, T, N+1);
nDOF = length(u0); U=zeros(nDOF,N+1,K+1);
U(:,1,:) = repmat(u0,1,1,K+1);              % set u^k_0 = u0
U(:,:,1) = uPred;                           % set u^0_n = uPred
for k=1:K
    for n=1:N
        % Propagate fine solution and Jacobian
        [uF, V] = prop(t(n), t(n+1), U(:,n,k));
        % Newton correction
        U(:,n+1,k+1) = uF + V*(U(:,n,k+1) - U(:,n,k));
    end
end
end

