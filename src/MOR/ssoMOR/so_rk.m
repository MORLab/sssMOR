function [sysr, V, R] = so_rk(sys, s0, Rt, IP, Opts)
% RK - Model Order Reduction by Rational Krylov
%
% Syntax:
%       sysr        = SO_RK(sys, s0_inp)
%       sysr        = SO_RK(sys, s0_inp, Rt, Opts)
%       [sysr, V, R]= SO_RK(sys, s0_inp,...)
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer, Alessandro Castagnotto 
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  06 Apr 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Parsing

% create the options structure
Def.real = true; %keep the projection matrices real?       
Def.lse = 'sparse'; 

if  exist('Opts','var') && ~isempty(Opts)
    Opts = parseOpts(Opts,Def);
else
    Opts = Def;
end       

%%  Check the inputs

if exist('s0', 'var')
    s0 = shiftVec(s0);
    % sort expansion points & tangential directions
    s0old = s0;
    if Opts.real, 
        s0 = cplxpair(s0); %make sure shifts can be paired 
    else
        s0 = sort(s0);
    end
    if exist('Rt','var') && ~isempty(Rt)
        if size(Rt,2) ~= length(s0),error('Inconsistent size of Rt');end
        
        [~,cplxSorting] = ismember(s0,s0old); 
        Rt = Rt(:,cplxSorting);
    else
        Rt = [];
    end
clear s0old
else
    error('specify s0');
end

if length(s0)> sys.n
    error('sssMOR:arnoldi:reducedOrderExceedsOriginal','The desired reduced order exceeds the original order');
end
%%  Define execution variables
if ~exist('IP', 'var'), 
    IP=@(x,y) (x.'*y);
end
%%  Computation

% input Krylov subspace

% SISO Arnoldi
[V, R] = so_arnoldi(sys, s0, Rt, IP, Opts);
W = V;
sysr = projectiveMor(sys,V,W);
