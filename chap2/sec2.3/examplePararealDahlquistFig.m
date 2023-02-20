la=2*1i;                                              % Dahlquist rhs
f=@(t,x) la*x;
MF=20; MG=1;                                        % F and G time steps
F=@(t0,t1,u0) SDahlquistBE(la,t0,t1,u0,MF);          % fine solver F
G=@(t0,t1,u0) SDahlquistBE(la,t0,t1,u0,MG);         % coarse solver G
K=10; u0=1; 
T=10; 
N=10;                          
U=Parareal(F,G,T,u0,N,K);                                
[t,u]=DahlquistBE(la,[0 T],u0,MF*N);                  % fine solution
TT=0:T/N:T;
for k=1:K                                           % plot fine
  plot(t,u,'-b',TT,U{k},'o');                      % solution and
  xlabel('t');                                     % parareal iterate
  set(gca,'Fontsize',18)
  legend('u','parareal')
  print('-depsc',['PararealDahlquistImT10Iter' num2str(k) '.eps']);
  pause
end

for k=1:K+1                             
  err(k)=max(abs(u(1:MF:end)-U{k}));
end      
DT=T/N; 
R0=abs(1/(1-la*DT));
if R0<1, R0=1; end;
errsup(1)=err(1);
errlin(1)=err(1);
for k=1:K                             
  errsup(k+1)=err(1)*abs(exp(la*DT)-1/(1-la*DT))^k/factorial(k)...
              *R0^(N-k-1)*prod(N-(1:k));
  errlin(k+1)=err(1)*(abs(exp(la*DT)-1/(1-la*DT))/(1-abs(1/(1-la*DT))))^k;
end      


semilogy(0:K,err,'--',0:K,errsup,'-',0:K,errlin,'-')
set(gca,'Fontsize',12)
%axis([0 10 1e-15 1])
axis([0 10 1e-10 1e5])
xlabel('k')          
ylabel('error')
legend('parareal error','superlinear bound','linear bound')
%print -depsc PararealDahlquistErrorImT1.eps
