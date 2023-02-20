la=-1; f=@(t,x) la*x;                               % Dahlquist rhs
MF=20; MG=1;                                        % F and G time steps
F=@(t0,t1,u0) SDahlquistBE(la,t0,t1,u0,MF);         % fine solver F
G=@(t0,t1,u0) SDahlquistBE(la,t0,t1,u0,MG);         % coarse solver G
u0=1; K=10; T=1; N=10;                              % parareal parameters
U=Parareal(F,G,T,u0,N,K);                           % apply parareal     
[t,u]=DahlquistBE(la,[0 T],u0,MF*N);                % fine solution
for k=1:K+1                                         
  err(k)=max(abs(u(1:MF:end)-U{k}));                % compute error
end      
DT=T/N; R0=abs(1/(1-la*DT)); if R0<1, R0=1; end;    % compute error bounds
errsup(1)=err(1); errlin(1)=err(1);
for k=1:K                             
  errsup(k+1)=err(1)*abs(exp(la*DT)-1/(1-la*DT))^k/factorial(k)...
              *R0^(N-k-1)*prod(N-(1:k));
  errlin(k+1)=err(1)*(abs(exp(la*DT)-1/(1-la*DT))/(1-abs(1/(1-la*DT))))^k;
end      
semilogy(0:K,err,'--',0:K,errsup,'-',0:K,errlin,'-')
xlabel('k'); ylabel('error'); 
legend('parareal error','superlinear bound','linear bound')