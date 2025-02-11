m=4; n=2^m;                       % number of grid points power of 2
x0=1; lambda=-1; dt=1/n;          % solve Dahlquist equation
d=ones(n,1); dm=-(1+lambda*dt)*d;
A{1}=spdiags([dm d],[-1 0],n,n);  % original time stepping matrix
f{1}=zeros(n,1); f{1}(1)=x0*(1+lambda*dt); % initial data
for i=1:m                         % cyclic reduction to 2x2
  [A{i+1},f{i+1}]=CyclicReduction(A{i},f{i});
end;
x{m}=A{m}\f{m};                   % solve smallest system
for i=m-1:-1:1                    % cyclic back substitution
  x{i}=CyclicBackSubstitution(A{i},f{i},x{i+1});
end;
t=0:dt:1;                         % compare with exact solution
plot(t,[x0;x{1}],'--',t,exp(lambda*t),'-')
xlabel('t'); legend('cyclic reduction','exact solution');