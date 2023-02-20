function [U,uTraj]=Nievergelt(solver,T,u0,N,uPred,Mn,width)
% NIEVERGELT implementation of the method of Nievergelt for scalar ODEs
%   [U,uTraj]=Nievergelt(solver,T,u0,N,uPred,Mn,width); solve a scalar
%   ODE with Nievergelt's method on N time sub-intervals dividing
%   (0,T]. u0 is the initial solution of the ODE; solver is a function
%   of the form u=solver(t0, t1, u0) that solves numerically the ODE
%   between t0 and t1 with u0 as initial value, and returns the
%   solution over the whole interval [t0,t1]; uPred is the cheap
%   predicted solution, Mn is the number of initial values for each
%   subinterval, and width is the width of the interval, centered
%   around U^0_n, on which the initial values are uniformely
%   distributed.

dT=T/N; TT=0:dT:T;
uTraj(1,1,:)=solver(TT(1),TT(2),u0);            % compute first trajectory
for n=2:N                                       % remaining trajectories
  uStart=uPred(n)+linspace(-width/2,width/2,Mn);% starting around prediction
  for m=1:Mn
    uTraj(n,m,:)=solver(TT(n),TT(n+1),uStart(m));
  end
end
U(1)=u0;                                        % compute Nievergelt solution
U(2)=uTraj(1,1,end);                            % fine solution for U_1^1
for n=2:N                                       % determine p and U_n^1
  m=FindClosestInterval(U(n),uTraj(n,:,1));
  uS1=uTraj(n,m,1); uS2=uTraj(n,m+1,1);
  p=(U(n)-uS2)/(uS1-uS2);
  U(n+1)=p*uTraj(n,m,end)+(1-p)*uTraj(n,m+1,end);
end
