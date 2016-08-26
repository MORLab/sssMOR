function [sysr, V, W, sOpt] = rkIcop(sys, s0, q, varargin)
% RKICOP - Rational Krylov with an iteratively calculated optimal point
%
% Syntax:
%       [sysr,V,W,sOpt] = RKICOP(sys,s0,q,Opts)
%
% Description:
%       This function iteratively deteremines an optimal expansion point 
%       for Krylov-based model order reduction. 
%
%       The computed point is obtained by iterating between the optimal
%       parameter and the reduced system starting from an inital parameter.
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
%						[{20} / positive integer]
%			-.tol:		convergence tolerance;
%						[{1e-2} / positive float]
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
%       * *[1] R. Eid, H. Panzer and B. Lohmann: How to choose a single 
%               expansion point in Krylov-based model reduction? Technical 
%               reports on Automatic Control, vol. TRAC-4(2), Institute of 
%               Automatic Control, Technische Universitaet Muenchen, Nov. 
%               2009.
%       * *[2] R. Eid: Time domain Model Reduction by Moment Matching, Ph.D
%               thesis, Institute of Automatic Control, Technische 
%               Universitaet Muenchen, 2009.
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

%% Parse the inputs
%   Default execution parameters
Def.rk = 'twoSided'; % 'twoSided','input','output'
Def.tol = 1e-2; % stopping tolerance
Def.maxIter = 100; % maximum number of iterations

% create the options structure
if ~exist('Opts','var') || isempty(fieldnames(Opts))
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

if isa(sys,'ss')
    sys=sss(sys);
end
sOpt=s0*ones(sys.p,sys.m);

if q<10
    warning('The results may be unprecise for small q.');
end

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
    sOpt = rkOp(sysr);
    
    if abs(sOpt-sOptOld)/sOpt <= Opts.tol
        break
    end
    if i==Opts.maxIter
        error(['rkIcop has not converged after ' num2str(Opts.maxIter) ' steps.']);
    end
end

