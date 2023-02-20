function m=FindClosestInterval(x,y)
% FINDCLOSESTINTERVAL find closest interval in a vector to a number
%   [m1,m2]=FindClosestInterval(x,y) finds in the vector y the
%   interval which contains the value x, or if there is no such
%   interval, the one which is closest to x, and returns the lower
%   index m of this interval in y.

n=length(y); i=1;
while i<=n & x>y(i), i=i+1; end;
if i==1, m=1; elseif i==n+1, m=n-1; else m=i-1; end;