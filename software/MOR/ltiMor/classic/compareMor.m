function RedData = compareMor(sys,redFct,Opts)
% COMPAREMOR - compare different model reduction algorithms
%
% Syntax:
%   COMPAREMOR(sys,redFct)
%   COMPAREMOR(sys,redFct,Opts)
%   RedData = COMPAREMOR(sys,redFct,Opts)
%
% Description:
%   This function executes the reduction scheme(s) defined by the handle 
%   (or a cell array of handles) redFct on the model sys. Optional
%   execution parameters can be defined in the structure Opts.
%
%   sys is a dynamical system given as sss object
%
%   redFct is a function handle (or a cell array thereof) that accepts
%       - a dynamical system sys
%   and returns:
%       - a reduced order dynamical system sysr
%
%   If no output is defined, only a graphical representation of the results
%   will be given. RedData is a structure contaning all information
%   regarding the reduction execution and results.
%
%
%   All other parameters the function may need should be defined when
%   creating the handle.
%
% See also:
%   SSS
%
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Last Change:  02 Sep 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%%  Input parsing
%   Make redFct a cell even if it is a single handle
if ~iscell(redFct), redFct = {redFct}; end

%   Default Opts
Def.verbose = 1;
Def.plotres = 1;

if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end  

%   Initializing output
RedData.timestamp               = datestr(now); 
RedData.pc                      = getComputerName;
RedData.data(length(redFct))    = struct('redFct','','sysr','','tRed','');

%%  Computations

if Opts.verbose, fprintf('Starting compareMor...\n'); end

for iRedFct = 1:length(redFct)
    CurrData = struct();
    CurrData.redFct = func2str(redFct{iRedFct});
    if Opts.verbose, 
        fprintf('Reducing with %s...\n', CurrData.redFct);
    end
    
    % reduction
    tic, CurrData.sysr = redFct{iRedFct}(sys); CurrData.tRed = toc;
    if Opts.verbose
        if CurrData.tRed/3600/24 >1 %more than one day
            tRedString = sprintf('%05.2f %s',CurrData.tRed/3600/24,'days');
        else
            tRedString = datestr(CurrData.tRed/3600/24,'HH:MM:SS.FFF');
        end 
        fprintf('...............done! (time elapsed: %s)\n',tRedString);
    end
    
    RedData.data(iRedFct) = CurrData;
end