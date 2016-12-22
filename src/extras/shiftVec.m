function s0=shiftVec(s0)
% SHIFTVEC - convert shift definition to single row notation
% 
% Syntax:
%		s0				= SHIFTVEC(s0)
% 
% Description:
%       Converts a two-raw matrix with [s1, s2, s3, ...; q1, q2, q3, ...]
%       where si denotes the value of the complex shift and qi its
%       multiplicity to a one row notation where the shift value si is
%       repeated qi times, e.g. [s1,s1,s1,s2,s3,s3,...]
%
% Input Arguments:
%		*Required Input Arguments:*
%		-s0: 		shift vector or matrix
%
% Output Arguments:
%       -s0: shift vector 
%
%
% See Also: 
%		rk, irka, arnoldi
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/?sssMOR">www.rt.mw.tum.de/?sssMOR</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sssMOR">www.rt.mw.tum.de/?sssMOR</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  16 Dec 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

    if size(s0,1)==2
        temp=zeros(1,sum(s0(2,:)));
        for j=1:size(s0,2)
            k=sum(s0(2,1:(j-1))); k=(k+1):(k+s0(2,j));
            temp(k)=s0(1,j)*ones(1,s0(2,j));
        end
        s0=temp;
    end

    if size(s0,1)>size(s0,2)
        s0=transpose(s0);
    end