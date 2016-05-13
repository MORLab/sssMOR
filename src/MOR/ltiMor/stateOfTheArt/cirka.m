function [sysr, s0] = cirka(sys, s0, Opts) 
% cirka - TODO
%
% Syntax:
%       [sysr, s0] 		= cirka(sys, s0) 
%	   [sysr, s0] 		= cirka(sys, s0, Opts) 
%
% Description:
%	   TODO
%
% Input Arguments:
%	   TODO
%
% Output Arguments:
%	   TODO
%
% See Also: 
%       TODO
%
% References:
%	   TODO
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
% Authors:      TODO
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  03 Mai 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------


%% Define auxiliary functions
redFct      = @(sys,s0) irka(sys,s0,struct('stopCrit','s0','suppressverbose',1));
redFctOut   = @(sys,s0) getDesiredOutput(redFct,[1,4],sys,s0);

%% Define execution options
Def.qm0     = length(s0)+2;
Def.s0m     = zeros(1,Def.qm0); 
Def.maxiter = 20; Def.tol = 1e-3;
Def.verbose = 0; Def.plot = 0;
Def.updateModel = 'all';

if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end 

%% Perform reduction
[sysr,s0] = modelFctMor(sys,redFctOut,s0,Opts);