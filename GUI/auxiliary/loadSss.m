function sys = loadSss(fname,Opts)
% LOADSSS - Creates an sss-object from .mat file data
%
% Syntax:
%       ssy = LOADSSS(fname,Opts)
%
% Description:
%       loadSss lets you define the path to a .mat file in which 
%       the system matrices are stored, either in the form (A,B,C,D,E) - 1st
%       order - or (M,Da,K,B,C) - 2nd order.
%
%       If not all of the matrices are found in the .mat file, following
%       assumptions will be applied:
% 
%       * C = B'
%       * D = zeros(size(C,1),size(B,2))
%       * E = speye(size(A))
% 
%       The conversion from 2nd order to 1st order gives some additional 
%       design freedom that can be parially selected by passing the 
%       appropriate options in the Opts structure.
%
%       The output is an sss-object containing the 1st order system matrices
%       as well as other information about the system that might be stored
%       in the .mat file.
%
% Input Arguments:
%       -fname: string containing either the name of a benchmark
%       distributed with the sss class or the full path to the .mat file
%       containing the matrices
%       -Opts: a structure containing the following options
%           -transf2nd: Type of transformation from 2nd to
%                       1st order form: [{'I'}, 'K', '-K', 'alpha']
%
% Output Arguments:
%       -sys: sss object of the benchmark system
%
% Examples:
%       The following code generates an sss-object from the .mat file data
%       'rail_1357':
%
%> load rail_1357.mat
%> sys = sss(A,B,C,[],E);
%> sysLoadSss = loadSss('rail_1357')
%> sys.D==sysLoadSss.D;
%
%       Another example, this time using the 2nd order model 'gyro' which
%       is converted into 1st order using the default matrix 'I' as 
%       submatrix E1:
%
%> sys = loadSss('gyro')
%
% See Also:
%       sss, load
%
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sss">sss</a>, a Sparse State-Space and System Analysis 
% Toolbox developed at the Chair of Automatic Control in collaboration
% with the Chair of Thermofluid Dynamics, Technische Universitaet Muenchen. 
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  11 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

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
    if ~isfield(LoadData,'B'), 
        if ~isfield(LoadData,'b')
           error('This benchmark does not have a B matrix. Please load it manually'),
        else
           LoadData.B = LoadData.b; 
        end
    end
    if ~isfield(LoadData,'C')
        if ~isfield(LoadData,'c')
            LoadData.C = LoadData.B'; 
        else
            LoadData.C = LoadData.c;
        end
    end
    if ~isfield(LoadData,'D'), 
        LoadData.D = zeros(size(LoadData.C,1),size(LoadData.B,2)); 
    end
    if ~isfield(LoadData,'E'), LoadData.E = speye(size(LoadData.A)); end
    
    sys = sss(LoadData.A,LoadData.B,LoadData.C,LoadData.D,LoadData.E);
    
elseif isfield(LoadData,'M') %2nd order form
%     msgID = 'sssMOR:loadSss:2ndOrder';
    warning('The system is in 2nd order form and will be converted to 1st order.')
    
    if ~isfield(LoadData,'D'), LoadData.D = zeros(size(LoadData.K)); end
    
    if ~isfield(LoadData,'C')
       if isfield(LoadData,'c')
           LoadData.C = LoadData.c;
       end
    end
    
    if ~isfield(LoadData,'B')
       if isfield(LoadData,'b')
          LoadData.B = LoadData.b;
       end
    end
    
    %Use the function second2first to create the system
    
    sys = second2first(LoadData.M,LoadData.D,LoadData.K,LoadData.B,...
                       zeros(size(LoadData.C)),LoadData.C,Opts);
    
else
    error('loadSss was not able to determine the form of the given system');
end

%%  Set name for system
splittedName = strsplit(fname,'\');
    
if length(splittedName) >= 2
    fname = splittedName{1,end};
end

splittedName = strsplit(fname,'.');
fname = splittedName{1,1};

sys.Name = fname;

%%  Store additional data into the sss object
% TODO: change sss class to allow such attributes

% %   We do not store ANY data given in the .mat file but just the one that
% %   is recognized as being relevant for the sss object. This includes
% %   - w, mag, ph: frequency response data for bode plots
% %   - hsv: Hankel singular values
% %   - S, R: sparse cholesky factors of the Gramian matrices
% %   - lambda: all or some of the eigenvalues of the system
% 
% LoadData = rmfield(LoadData,{'A','B','C','D','E'});
%
% if ~isempty(fieldnames(LoadData))
%     additionalFields = {'w','mag','ph','hsv','S','R','lambda'};
%     for iFields = 1:length(additionalFields)
%         if isfield(LoadData,additionalFields{iFields})
%             sys.(additionalFields{iFields}) = LoadData.(additionalFields{iFields});
%         end
%     end
% end

    