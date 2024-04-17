function w=LagrangeWeights(M)
% LAGRANGEWEIGHTS Compute Lagrange weights
%   w=LagrangeWeights(M) Returns analytical Lagrange weights for M uniform
%   points.

weights{2}=[1/2, 1/2];
weights{3}=[
    [ 5/24, 1/3, -1/24];
    [-1/24, 1/3,  5/24]];
weights{4}=[
    [ 9/72, 19/72, -5/72,  1/72];
    [-1/72, 13/72, 13/72, -1/72];
    [ 1/72, -5/72, 19/72,  9/72]];
weights{5}=[
    [ 251/2880, 323/1440, -11/120,  53/1440, -19/2880];
    [ -19/2880, 173/1440,  19/120, -37/1440,  11/2880];
    [  11/2880, -37/1440,  19/120, 173/1440, -19/2880];
    [ -19/2880,  53/1440, -11/120, 323/1440, 251/2880]];
weights{6}=[
    [  19/288, 1427/7200, -399/3600,  241/3600, -173/7200,    3/800];
    [  -3/800,  637/7200,  511/3600, -129/3600,   77/7200, -11/7200];
    [ 11/7200,  -31/2400,  401/3600,  401/3600,  -31/2400,  11/7200];
    [-11/7200,   77/7200, -129/3600,  511/3600,  637/7200,   -3/800];
    [   3/800, -173/7200,  241/3600, -399/3600, 1427/7200,   19/288]];
w=weights{M};