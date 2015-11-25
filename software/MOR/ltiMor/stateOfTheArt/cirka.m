function [sysr, s0] = cirka(sys, s0, Opts) 

%% Define auxiliary functions
redFct      = @(sys,s0) irka(sys,s0,struct('stopCrit','s0','suppressverbose',1));
redFctOut   = @(sys,s0) getDesiredOutput(redFct,[1,4],sys,s0);

%% Define execution options
Def.qm0     = length(s0)+2;
Def.s0m     = zeros(1,Def.qm0); 
Def.maxiter = 20; Def.tol = 1e-3;
Def.verbose = 0; Def.plot = 0;
Def.updateModel = 'new';

if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end 

%% Perform reduction
[sysr,s0] = modelFctMor(sys,redFctOut,s0,Opts);