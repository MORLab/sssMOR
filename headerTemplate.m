function [output1,output2] = headerTemplate(input1,input2,input3)
% HEADERTEMPLATE - sssMOR toolbox header template
% 
% Syntax:
%       HEADERTEMPLATE
%       HEADERTEMPLATE(input1)
%       HEADERTEMPLATE(input1,input2)
%       [output1,output2] = HEADERTEMPLATE(input1,input2)
% 
% Description:
%       This header file is to be copied and pasted at the beginning
%       of every script and function within the MORLab toolbox for reference.
% 
%       The description should be kept as short as possible but as long as
%       required to give a good impression of what the function/script is 
%       about and how it should be used.
%
% Input Arguments:
%		-input1: Description
%       -input2: Description
%       -input3: Description
%
% Output Arguments:
%       -output1: Description 
%       -output2: Description
%
% Examples:
%		This code loads a benchmark model included in the toolbox
%		and plots its bode diagram using the sparse state space class:
%
%> load build; 
%> sys = sss(A,B,C);
%> bode(sys);
%
%// 	Note: the build.mat file for the example can be found in the benchmarks folder
%
% See Also: 
%		arnoldi
%
% References:
%		* *[1] Grimme (1997)*, Krylov projection methods for model reduction
%		* *[2] Antoulas (2005)*, Approximation of large-scale dynamical systems
%		* *[3] Giraud (2005)*, The loss of orthogonality in the Gram-Schmidt...
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Author 1, Author 2
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  30 Oct 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------