function sys = second2first(M, D, K, B, Cx, Cv, Opts)
%SECOND2FIRST - Convert a 2nd order system to state space representation
%
% Syntax:
%       sys = SECOND2FIRST(M,D,K,B,Cx,Cv)
%       sys = SECOND2FIRST(M,D,K,B,Cx,Cv,Opts)
%
% Description:
%       This function takes the system matrices of a 2nd order system and
%       converts this system to state space representation. 
%       
%       The representation of a 2nd order system thereby is defined as
%       follows:
%
%       $M  ~ \ddot x + D ~ \dot x + K ~ x = B ~ u$
%
%       $y = C_x ~ x + C_v ~ \dot x$
%
%       State space representation uses the following form:
%
%       ${ F ~~ 0 \choose 0 ~~ M}$  ${\dot x \choose \ddot x}$
%       $= {0 ~~ F \choose -K ~~ -D}$ ${x \choose \dot x }$ + $B ~ u$
%
%       $y = \left( C_x ~~ C_v \right) ~ { x \choose \dot x }$
%
%       The conversion from second order to first order allows some freedom
%       in the choise of the matrix F. Which value from the possible
%       options should be used to create the matrix F can be specified by
%       passing a option-structure with the desired choise as an additional
%       input argument to the function. 
%
%       If the representation of the system to convert contains some
%       matrices which are all zero (i.e. no damping), then the value [] 
%       can be specified for these matrices.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -M:     Mass matrix of the second order system
%       -D:     Damping matrix of the second order system
%       -K:     Stiffness matrix of the second order system
%       -B:     Input matrix of the second order system
%       -Cx:    Output matrix of the second order system 
%       -Cv:    Output matrix of the second order system
%       *Optional Input Arguments:*
%       -Opts:	a structure containing following options
%			-transf2nd:     Type of transformation from 2nd to
%                           1st order form: 
%                           [{'I'}, 'K', '-K', 'alpha']
%
% Output Arguments:
%       -sys: Sparse State Space Model (1st order form)
%
% Examples:
%       The following code converts a randomly created 2nd order system
%       without damping to a 1st order system:
%> clear all;
%> sys = second2first(rand(20),[],rand(20),rand(20,1),rand(1,20),[]);
%	
% See Also:
%       loadSss
%
% References:
%       * *[1] B. Salimbahrami*: "Structure Preserving Order Reduction of Large
%       Scale Second Order Models", PhD Thesis at TUM, p. 36f
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sss">sss</a>, a Sparse State-Space and System Analysis 
% Toolbox developed at the Chair of Automatic Control in collaboration
% with the Professur fuer Thermofluiddynamik, Technische Universitaet Muenchen. 
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sss">sss</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer, Rudy Eid 
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  21 Mar 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

% check wheather options are specified or not 
if ~exist('Opts','var')
    Opts = struct(); 
end

Def = struct('transf2nd','I');
Opts = parseOpts(Opts,Def);


% check the matrix-dimensions of M and D
if size(M,1)~=size(M,2)
    error('M must be symmetric.')
elseif any(size(M)-size(K))
    error('M and K must have same size.');
end

if isempty(D) || ~any(any(D))
    D = sparse(size(M,1), size(M,1));
elseif any(size(D)-size(K))
    error('D must have same size as M and K.');
end

% evaluate the specified option
switch Opts.transf2nd
        case 'I'
            F = speye(size(K));
        case 'K'
            F = K;
        case '-K'
            F = -K;
        case 'alpha'
            % TODO: add the transformation to strictly dissipative form by
            % Panzer
            warning('alpha option not implemented yet, using K instead');
            alpha = 1; 
            F = alpha*LoadData.K;
end

% create the matrices A and E for the first-order-system
n = size(M,1);
O = sparse(n,n);

E = [F O; O M];
A = [O F; -K -D];

% check the matrix-dimensions of Cx and Cv and B
if isempty(Cx) && isempty(Cv)
    % both output vectors are empty
    error('All outputs are zero.');
elseif isempty(Cx)
    % only Cv is given
    Cx = sparse(size(Cv,1),n);
elseif ~exist('Cv', 'var') || isempty(Cv)
    % only Cx is given
    Cv = sparse(size(Cx,1),n);
elseif any(size(Cx)-size(Cv))
    error('Cx and Cv must have same size.');
elseif size(Cx,2)~=size(M,2)
    error('Cx, Cv must have same column dimension as M, D, K')
end

if isempty(B)
    error('All inputs zero.');
elseif size(B,1)~=size(M,1)
    error('B must have same row dimension as M, D, K')
end

% create the matrices C, D and B for the first-order-system 
B = [sparse(n,size(B,2)); B];
C = [Cx, Cv];
D = zeros(size(C,1),size(B,2));

% create the first-order-system
sys = sss(A, B, C, D, E);
