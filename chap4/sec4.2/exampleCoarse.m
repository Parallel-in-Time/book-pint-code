Jc=(J+1)/2-1;                                  % coarse grid points
P=sparse(J,Jc);                                % prolongation by
for j=1:Jc                                     % interpolation
  P(2*j,j)=1; P(2*j-1,j)=0.5; P(2*j+1,j)=0.5;
end;
R=0.5*P';                                      % restriction by transpose
Lc=R*L*P;                                      % coarse matrix by Galerkin
Ac=speye(size(Lc))-dt*Lc;                      % coarsening in space only
xc=x(1:2:end);                                 % coarse spatial mesh