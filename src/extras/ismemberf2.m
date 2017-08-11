function idx = ismemberf2(sNew,sOld,tol)
% ISMEMBERF2 - determine which elements are within a radius
%
% Syntax:
%       idx = ISMEMBERF2(sNew,sOld,tol)
%
% Description:
%       Determine which elements of a vector sNew are already included (within
%       a relative tolerance tol) in a comparison set sOld.
% 
% Input Arguments:  
%       *Required Input Arguments:*
%       -sNew:			new set
%       -sOld:			old set as comparison basis
%       -tol:           tolerance radius (relative to values)
%
%
% See also:
%   ismember, modelFct, cirka
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/?sssMOR">www.rt.mw.tum.de/?sssMOR</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sssMOR">www.rt.mw.tum.de/?sssMOR</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  22 Nov 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%   doing it with loops, there might be a faster version

idx = zeros(1,length(sNew));
if true
    %implementation 1: for loop
    for iS = 1:length(sNew)
        dist = abs(sOld-sNew(iS));
        idx(iS) = any(dist<abs(sOld)*tol);
    end
else
%     %implementation 2: vectorized
%     [Snew,Sold] = meshgrid(sNew,sOld);
%     dist = abs(Sold-Snew);
%     idx = any(dist<abs(sOld)*tol);
end