function [sysr, V, W, sOpt] = rkIcop(sys, s0, q, varargin)
% RKICOP - Rational Krylov with an iteratively calculated optimal point
%
% Syntax:
%       [sysr,V,W,sOpt] = RKICOP(sys,s0,q,Opts)
%
% Description:
%       This function iteratively determines an optimal expansion point 
%       for Krylov-based model order reduction. 
%
%       The computed point is obtained by iterating between the optimal
%       parameter and the reduced system starting from an initial parameter.
%
% Input Arguments:  
%       *Required Input Arguments:*
%       -sys:			full oder model (ss or sss)
%       -s0:            inital expansion point (scalar)
%       -q:             reduction order
%       *Optional Input Arguments:*
%       -Opts:			structure with execution parameters
%			-.rk:       reduction type
%						[{'twoSided'} / 'input' / 'output']
%			-.maxIter:	maximum number of iterations;
%						[{100} / positive integer]
%			-.tol:		convergence tolerance;
%						[{1e-2} / positive float]
%           -.lse:      solve linear system of equations
%                       [{'sparse'} / 'full' / 'gauss' / 'hess' / 'iterative' ]
%
% Output Arguments:
%       -sysr:          reduced system
%       -V, W:          projection matrices spanning Krylov subspaces
%       -sOpt:          optimal expansion point
%
% Examples:
%       This code iteratively approximates the optimal expansion point for 
%       an initial starting point s0=1 and computes a reduced system of 
%       order q=10 for the benchmark model 'building'.
%
%> sys = loadSss('building')
%> [sysr,V,W,sOpt] = rkIcop(sys,1,10);
%> bode(sys,'-',sysr,'--r');
%
% See Also: 
%       rkOp, rk, irka, arnoldi
%
% References:
%       * *[1] Eid, Panzer and Lohmann (2009)*, How to choose a single 
%              expansion point in Krylov-based model reduction? Technical 
%              reports on Automatic Control, vol. TRAC-4(2), Institute of 
%              Automatic Control, Technische Universitaet Muenchen. 
%
%       * *[2] Eid (2009)*, Time domain Model Reduction by Moment Matching, Ph.D
%              thesis, Institute of Automatic Control, Technische 
%              Universitaet Muenchen.
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
% Authors:      Heiko Panzer (heiko@mytum.de), Rudy Eid
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  28 Jun 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

if ~isempty(varargin) && isa(varargin{1},'struct')
    Opts=varargin{1};
end

warning('off','sss:sss:ssRedConversion')
%% Parse the inputs
%   Default execution parameters
Def.rk      = 'twoSided'; % 'twoSided','input','output'
Def.tol     = 1e-2; % stopping tolerance
Def.maxIter = 100; % maximum number of iterations
Def.lse     = 'sparse'; % 'sparse', 'full', 'hess', 'gauss', 'iterative'

% create the options structure
if ~exist('Opts','var') || isempty(fieldnames(Opts))
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

sOpt=s0*ones(sys.p,sys.m);

if q<10
    warning('The results may be unprecise for small q.');
end

%% rkIcop iteration
for i=1:Opts.maxIter
    sOptOld=sOpt;
    
    % calculate reduced system
    if sys.isSiso
        switch(Opts.rk)
            case 'twoSided'
                [sysr, V, W] = rk(sys, [sOpt;q], [sOpt;q]);
            case 'input'
                [sysr, V, W] = rk(sys, [sOpt;q]);
            case 'output'
                [sysr, V, W] = rk(sys, [], [sOpt;q]);
            otherwise
                error('Wrong Opts.');
        end
    else
        switch(Opts.rk)
            case 'twoSided'
                sOpt=sOpt.';
                tempLt=[];
                Rt=[];
                Lt=[];

                for j=1:sys.m
                   Rt=blkdiag(Rt,ones(1,q*sys.p));
                end
                
                for j=1:sys.p
                    tempLt=blkdiag(tempLt,ones(1,q));
                end                    

                for j=1:sys.m
                    Lt=[Lt,tempLt];
                end
                
                [sysr,V,W] = rk(sys,[sOpt(:).';ones(1,sys.m*sys.p)*q],[sOpt(:).';ones(1,sys.m*sys.p)*q],Rt,Lt);
            case 'input'
                sOpt=sOpt.';
                Rt=[];

                for j=1:sys.m
                   Rt=blkdiag(Rt,ones(1,q*sys.p));
                end
                
                [sysr,V,W] = rk(sys,[sOpt(:).';ones(1,sys.m*sys.p)*q],Rt);
            case 'output'
                sOpt=sOpt.';
                tempLt=[];
                Lt=[];
                
                for j=1:sys.p
                    tempLt=blkdiag(tempLt,ones(1,q));
                end                    

                for j=1:sys.m
                    Lt=[Lt,tempLt];
                end
                [sysr,V,W] = rk(sys,[],[sOpt(:).';ones(1,sys.m*sys.p)*q],[],Lt);
            otherwise
                error('Wrong Opts.');
        end
    end
    
    % calculate sOpt
    if isstable(sysr)
        sOpt = rkOp(sysr, Opts);
    else
        warning('sssMOR:rkIcop:sysrUnstable','The ROM in rkIcop was unstable. Using stabsep')
        sOpt = rkOp(stabsep(sysr), Opts);
    end
    
    if abs(sOpt-sOptOld)/sOpt <= Opts.tol
        break
    end
    if i==Opts.maxIter
        warning('sssMor:rkIcop:NotConverged',['rkIcop has not converged after ' num2str(Opts.maxIter) ' steps.']);
    end  
end

%%  Storing additional parameters
%Stroring additional information about thr reduction in the object 
%containing the reduced model:
%   1. Define a new field for the Opts struct and write the information
%      that should be stored to this field
%   2. Adapt the method "parseParamsStruct" of the class "ssRed" in such a
%      way that the new defined field passes the check
Opts.originalOrder = sys.n;
if ~isfield(Opts,'rk') Opts.orth = 'twoSided'; end
if ~isfield(Opts,'lse') Opts.lse = 'sparse'; end
if ~isfield(Opts,'maxIter') Opts.maxIter = 100; end
if ~isfield(Opts,'tol') Opts.tol = 1e-2; end
Opts.s0 = s0;
Opts.sOpt = sOpt;
sysr = ssRed(sysr.A,sysr.B,sysr.C,sysr.D,sysr.E,'rkIcop',Opts,sys);

warning('on','sss:sss:ssRedConversion')


