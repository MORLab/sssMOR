function sys = loadSss(fname,Opts)
%LOADSSS - create sss object from .mat file data
% ------------------------------------------------------------------
% ssys = loadSss(fname,Opts)
% Inputs:       * fname : string containing either the name of a
%                         benchmark distributed with the sss class or
%                         the full path to the mat file containing 
%                         the matrices
%               * Opts: a structure containing following options
%                   -transf2nd: Type of transfromation from 2nd to
%                       1st order form {'I', 'K', '-K', 'alpha'}
% TODO: add Opts
% Outputs:      * sys: sss object of the benchmark system
% ------------------------------------------------------------------
% USAGE:  loadSss lets you define the path to a .mat file in which 
% the system matrices are stored, either in the form (A,B,C,D,E) - 1st
% order - or (M,Da,K,B,C) - 2nd order.
%
% If not all of the matrices are found in the .mat file, following
% assumptions will be applied:
%   - C = B'
%   - D = zeros(size(C,1),size(B,2))
%   - E = speye(size(A))
%
% The conversion from 2nd order to 1st order gives some additional 
% design freedom that can be parially selected by passing the 
% appropriate options in the Opts structure.
%
% The output is an sss-object containing the 1st order system matrices
% as well as other information about the system that might be stored
% in the .mat file.
%
% See also SSS
%
% ------------------------------------------------------------------
% REFERENCES:
% TODO: add references to the technical report or paper describing sssMOR
% ------------------------------------------------------------------
% This file is part of sssMOR, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Institute 
% of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit www.rt.mw.tum.de
% For any suggestions, submission and/or bug reports, mail us at
%                      -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Last Change:  11 Aug 2015
% Copyright 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%%  Parse input and define default options
if nargin == 0, error('No benchmark system defined'),
else
    if ~exist('Opts','var'), Opts = struct(); end %preallocate
    Def = struct('transf2nd','I');
    Opts = parseOpts(Opts,Def);
end

%%  Load information from the .mat file and transform to 1st order
LoadData = load(fname);

if isfield(LoadData,'A') %1st order form
    if ~isfield(LoadData,'C'), LoadData.C = LoadData.B'; end
    if ~isfield(LoadData,'D'), 
        LoadData.D = zeros(size(LoadData.C,1),size(LoadData.B,2)); 
    end
    if ~isfield(LoadData,'E'), LoadData.E = speye(size(LoadData.A)); end
    
elseif isfield(LoadData,'M') %2nd order form
    warning('The system is in 2nd order form and will be converted to 1st order.')
    if ~isfield(LoadData,'D'), LoadData.D = zeros(size(LoadData.K)); end
    switch Opts.transf2nd %create the matrix E1 multiplying velocities
        case 'I'
            E1 = speye(size(LoadData.K));
        case 'K'
            E1 = LoadData.K;
        case '-K'
            E1 = -LoadData.K;
        case 'alpha'
            % TODO: add the transformation to strictly dissipative form by
            % Panzer
            alpha = 1; 
            E1 = -alpha*LoadData.K;
    end
    % generate 1st order system matrices
    LoadData.E = blkdiag(E1,LoadData.M);
    LoadData.A = [zeros(size(LoadData.K)), E1;
                  -LoadData.K, -LoadData.D]; clear E1
    LoadData.B = [zeros(size(LoadData.B)); LoadData.B];
    LoadData.C = [zeros(size(LoadData.C)), LoadData.C];
    LoadData.D = zeros(size(LoadData.C,1),size(LoadData.B,2)); %overwrite with feedthrough
else
    error('loadSss was not able to determine the form of the given system');
end

%%  Create sss-object
sys = sss(LoadData.A,LoadData.B,LoadData.C,LoadData.D,LoadData.E);
LoadData = rmfield(LoadData,{'A','B','C','D','E'});

%%  Store additional data into the sss object
% TODO: change sss class to allow such attributes

% %   We do not store ANY data given in the .mat file but just the one that
% %   is recognized as being relevant for the sss object. This includes
% %   - w, mag, ph: frequency response data for bode plots
% %   - hsv: Hankel singular values
% %   - S, R: sparse cholesky factors of the Gramian matrices
% %   - lambda: all or some of the eigenvalues of the system
% 
% if ~isempty(fieldnames(LoadData))
%     additionalFields = {'w','mag','ph','hsv','S','R','lambda'};
%     for iFields = 1:length(additionalFields)
%         if isfield(LoadData,additionalFields{iFields})
%             sys.(additionalFields{iFields}) = LoadData.(additionalFields{iFields});
%         end
%     end
% end

    