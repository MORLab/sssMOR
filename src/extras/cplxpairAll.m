function [s0out,varargout] = cplxpairAll(s0in,varargin)
% CPLXPAIRALL - Sort arrays of complex numbers into complex conjugate pairs
% 
% Syntax:
%		s0out				= CPLXPAIRALL(s0out)
%		s0out,A1,A2,... 	= CPLXPAIRALL(s0out,A1,A2,...)
% 
% Description:
%       This function sorts the first input |s0out| by using the function
%       |cplxpair|. If several inputs arrays A1,A2,...,An are passed, then,
%       in addition, also A1,A2,...,An will be sorted according to the new
%       sorting of s0out.
%
%       This function is useful when sorting different arrays according to
%       one, as |cplxpair| does not return the respective indices.
%
% Input Arguments:
%		*Required Input Arguments:*
%		-s0in:          Array of complex conjugate pairs
%
%		*Optional Input Arguments:*
%		-A1,A2,...:   	Arrays of complex conjugate pairs
%
% Output Arguments:
%       -s0out:         Array sorted according to cplxpair
%       -A1,A2,...:     Arrays sorted according to s0out
%
% See Also: 
%		cplxpair, ismember
%
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
% Last Change:  02 Aug 2017
% Copyright (c) 2016,2017 Chair of Automatic Control, TU Muenchen


s0out = cplxpair(s0in);

if nargin > 1
    [~,cplxSorting] = ismember(s0out,s0in); %determine sorting order
    varargout = varargin; %preallocate
    for iArg = 1:length(varargin)
        varargout{iArg} = varargin{iArg}(:,cplxSorting);
    end
end