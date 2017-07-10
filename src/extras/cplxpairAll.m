function [s0out,varargout] = cplxpairAll(s0in,varargin)
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
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  20 Jan 2017
% Copyright (c) 2016,2017 Chair of Automatic Control, TU Muenchen


s0out = cplxpair(s0in);

if nargin > 1
    [~,cplxSorting] = ismember(s0out,s0in); %determine sorting order
    varargout = varargin; %preallocate
    for iArg = 1:length(varargin)
        varargout{iArg} = varargin{iArg}(:,cplxSorting);
    end
end