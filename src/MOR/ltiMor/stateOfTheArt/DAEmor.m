function [sysr,sysSPr,sysIMr] = DAEmor(sys,SpecProjFinC,SpecProjFinB,nu,Opts)
% Model Order Reduction of Structured DAEs
%
% Syntax:
%       sysr = DAEmor(sys,SpecProjFinC,SpecProjFinB,nu)
%       sysr = DAEmor(sys,SpecProjFinC,SpecProjFinB,nu,Opts)
%
% Description:
%       This function implements the reduction of structured DAEs by
%       1.) reduction of the strictly proper subsystem with CUREd SPARK
%       2.) finding a minimal realization of the improper subsystem by BT
%       
% Input Arguments:
%       *Required Input Arguments:*
%       -sys:           a sss-object containing the LTI-DAE-system
%       -SpecProjFinC:  function handle for projecting the output matrix C
%                       onto the deflating subspace of (lambda*E-A) corr.
%                       to the finite eigenvalues
%                       (for explicit knowledge of spectral projectors
%                       use: "C*SpecFinRight")
%       -SpecProjFinB:  function handle for projecting the input matrix B
%                       onto the deflating subspace of (lambda*E-A) corr.
%                       to the finite eigenvalues
%                       (for explicit knowledge of spectral projectors
%                       use: "SpecFinLeft*B")
%       -nu:            Index of the DAE-system
%       *Optional Input Arguments:*
%       -Opts: A structure containing options which are passed through
%
% Output Arguments:     
%       -sysr:   reduced (overall) system
%       -sysSPr: reduced strictly proper subsystem
%       -sysIMr: reduced improper subsystem
%
% Examples:
%
% See Also: 
%
% References:
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
% Authors:      Heiko Panzer, Alessandro Castagnotto, Maria Cruz Varona
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  13 Apr 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% Parse input and load default parameters
% default values
Def.cure.fact = 'V'; %error factorization

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

%% Separation of strictly proper and improper subsystem
% create strictly proper subsystem
sysSP = sys;

% determine method (input-PORK 'V' / output-PORK 'W')
switch Opts.cure.fact
    case 'V'    % input PORK
        %partitioning by projection of output matrix C
        sysSP.C = SpecProjFinC(sysSP.C);
        Bim = sys.B - SpecProjFinB(sys.B);
        Cim = sys.C - sysSP.C;
        
    case 'W'    % output PORK
        %partitioning by projection of input matrix B
        sysSP.B = SpecProjFinB(sysSP.B);
        Bim = sys.B - sysSP.B;
        Cim = sys.C - SpecProjFinC(sys.C);
end

%% Processing of strictly proper subsystem
% run CUREd SPARK
sysSPr = cure(sysSP,Opts);

%% Processing of improper subsystem
% find minimal realization of polynomial part by balanced truncation
sysIMr = tbrDAEimproper(sys,nu,Bim,Cim,Opts);

%% Obtain the reduced model from summing the subsystems
sysr = sysSPr + sysIMr;
