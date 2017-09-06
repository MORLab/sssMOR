% Verifying equivalence of rational and cummulative rational Krylov
% subspaces
%
% Description:
%   The function irka.m is tested (5 tests) on:
%    + comparing sysr to a directly calculated solution based on RK.m and
%      the final s0
%    + k>kmax, error>epsilon
%    + Ar, Er purely real
%    + rank(Ar), rank(Er) full
%    + Neither Inf nor NaN in Ar, Er
%    + s0: purely real, purely imaginary, zero, Inf (Markov-parameter)
%    + test systems: building, beam, random, SpiralInductorPeec 
%      (includes E-matrix), LF10 (includes E-matrix).
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Paul Heidenreich
% Last Change:  09 Aug 2017
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%% still to come
% - vergleich mit basis Krylov subspace erstellen

%% compare subspaces
% read in system information
A = [1 2 2; 4 2 6; 7 8 9];      % A-matrix
E = [1 1 1; 1 1 1; 1 1 1];      % E-matrix
%B = [1 2; 4 5; 6 7];            % B-Matrix for MIMO
B = [1; 2; 5];                  % b-vector for SISO
s = [1 2 9];                % shift-vector

% initialize V-matrix
n = size(A,1);      m = size(B,2);        
V_frt = zeros(n,m*size(s,2));
V_sec = V_frt;

%% SISO case
% first method: K = span{A_s1_inv*b,A_s2_inv*b,A_s3_inv*b,...}

% set fixed rhs
rhsB = B;

for ii = 1:1:size(s,2)
    jCol = s(1,ii);
    Anew_V = (A-jCol*E);
    
    % new direction in V
     [vnew] = solveLse(Anew_V,rhsB);
     
    % build basis V
    jbasis = ii*m;
    V_frt(:,jbasis+1-m:ii*m) = vnew;   
end

% second method, recycling old directions: 
% K = span{A_s1_inv*b, A_s2_inv*E*A_s1_inv*b, A_s3_inv*A_s2_inv*E*A_s1_inv*b, ...}

% first direction
jCol = s(1,1);
Anew_V = (A-jCol*E);    [vnew] = solveLse(Anew_V,rhsB);
V_sec(:,1:m) = vnew;

for ii = 2:1:size(s,2)
    jCol = s(1,ii);
    
    jbasis = ii*m;
    rhsB = E*V_sec(:,jbasis-2*m+1:jbasis-m);
    Anew_V = (A-jCol*E);    [vnew] = solveLse(Anew_V,rhsB);
    V_sec(:,jbasis+1-m:ii*m) = vnew;  
end

% rank of both subspaces
rank_frt = rank(V_SISO_frt);
rank_sec = rank(V_SISO_sec);

% angle between both subspaces
alpha = subspace(V_SISO_frt,V_SISO_sec);


%% MIMO case























