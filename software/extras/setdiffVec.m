function setdiffVec = setdiffVec(v1,v2)
% SETDIFFVEC - compute difference between two unsorted vectors
%
% Syntax:
% 		setdiffVec = SETDIFFVEC(v1,v2)
%
% Description:
% 		This function computes the difference between the elements in two
% 		(generally complex valued) unsorted vectors v1 and v2. If the elements 
% 		in v1 and v2 are complex, then they must be closed under conjugation.
% 
% 		Sometimes we wish to verify if two sets of (complex) values are equal 
% 		and if not, compare some distance measure.
% 		MATLAB's minus function operates elementwise. However, if the sorting 
% 		is not known a priori, computing the right difference might be a 
% 		challengin task.
% 
% 		This function takes the vectors v1 and v2 as set of numbers, eliminates
% 		all entries that are equal (using the builtin setdiff function), and
% 		computes a difference vector of the same size by substracting the
% 		remaining entries after appropriate sorting.
% 
% 		// Warning: since the entries in the vectors are sorted, it is not
% 		possible to assign the entries of setdiffVec to any entries of v1 or
% 		v2. This function is hence useful only if the order of the result can
% 		be arbitrary.
%
% Example:
%		TODO description of what this code does
%
%> v1 = [38, 0.43-13.44i, 0.43+13.44i, 0.38-5.21i,0.38+5.21i]';
%> v2 = [0.38-5.21i, 0.38+5.21i, 0.43-13.44i, 0.43+13.44i, 38]';
%> diff1 = v1 - v2 %standard difference in MATLAB
%> diff2 = setdiffVec(v1,v2)
%> norm(diff1),norm(diff2)
%
% See also:
%		setdiff, cplxpair
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
% Last Change:  01 Sep 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Parse input and initialize
if length(v1)~=length(v2),error('Both vectors must be of the same length'),end

setdiffVec = zeros(size(v1));
%% Computations
%   Compute the entries that differ (all others have to match exactly)
v1diff = setdiff(v1,v2,'stable');
v2diff = setdiff(v2,v1,'stable');

%   Repeated entries are returned only once by setdiff. Make sure you get
%   the right multiplicity
idx = [];
for iElements = 1:length(v1diff),
    idx = [idx, find(v1==v1diff(iElements))];
end
v1diff = v1(idx);
idx = [];
for iElements = 1:length(v2diff),
    idx = [idx, find(v2==v2diff(iElements))];
end
v2diff = v2(idx);

%   Make some appropriate sorting
tol = 1e-6;
v1diff = cplxpair(v1diff,tol); %sort v1diff and v2diff by the same principle
v2diff = cplxpair(v2diff,tol);

%   Compute the difference
setdiffVec(1:length(v1diff)) = v1diff-v2diff;