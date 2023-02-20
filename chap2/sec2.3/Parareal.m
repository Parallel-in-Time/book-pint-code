function U=Parareal(F,G,T,u0,N,K);
% PARAREAL implementation of the parareal algorithm
%   U=Parareal(F,G,T,u0,N,K); applies the parareal algorithm with fine
%   solver F(t0,t1,ut0) and coarse solver G(t0,t1,ut0) on [0,T] with
%   initial condition u0 at t=0 using N equidistant coarse time points
%   doing K iterations. The output U{k} contains the parareal
%   approximations at the coarse time points for each iteration k.

dT=T/N; TT=0:dT:T;                               % coarse time mesh
U{1}(1,:)=u0;                                    
for n=1:N                                        % initial guess with G 
  Go(n+1,:)=G(TT(n),TT(n+1),U{1}(n,:));         
  U{1}(n+1,:)=Go(n+1,:);                         % keep Go for parareal 
end;
for k=1:K                                        % parareal iteration  
  for n=1:N                                  
    Fn(n+1,:)=F(TT(n),TT(n+1),U{k}(n,:));        % parallel with F
  end;
  U{k+1}(1,:)=u0;
  for n=1:N                                  
    Gn(n+1,:)=G(TT(n),TT(n+1),U{k+1}(n,:));      % sequential with G
    U{k+1}(n+1,:)=Fn(n+1,:)+Gn(n+1,:)-Go(n+1,:); % parareal update 
  end;
  Go=Gn;                                         % keep for next iteration
end;
