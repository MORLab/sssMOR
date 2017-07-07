function [sysr, V, Rt, s0, kIter] = so_irka(sys,q,Opts)
% Iterative-Rational-Krylov-Algorithmus f¨¹r
% Systeme zweiter Ordnung (SO-IRKA)
%
% Input    -      tol:     Tolernz der Hauptschleife   
%                  nr:     Reduktionsordnung
%                nmax:     Maximale Anzahl der Iterationen
%                   M:     Massenmatrix des Originalsystems
%                   K:     Steifigkeitsmatrix des Originalsystems
%                   D:     D?mpfungsmatrix des Originalsystems
%                   B:     Eingangsmatrix des Originalsystems
%
% 
% Output    -       V:     Rechte Projektionsmatrix
%                sysr:     Reduziertes Modell
%                iter:     Anzahl der Iterationen
%
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
% Authors:      Xuwei Wu, Alessandro Castagnotto
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  07 Jul 2017
% Copyright (c) 2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Define optional execution parameters
Def.tol     = 1e-6;
Def.nmax    = 50;

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end  

%% Initialize

kIter       = 0;
s0          = zeros(1,q);
Rt          = ones(sys.m,q);

es          = Inf; 

%% Run SO-IRKA iteration

while kIter <=Opts.nmax && es > Opts.tol
    kIter = kIter+1;
    
    % reduction
    [sysr, V]        = so_rk(sys,s0,Rt);

    % new shifts
    [X,e] = polyeig(full(sysr.K),full(sysr.D),full(sysr.M));   
    %make a selection of q out of 2q eigenvalues -> smallest magnitude
    [~,idx] = sort(e,'ascend');
    
    s0new = -e(idx(1:q)).'; 
    Rt = (X(:,idx(1:q)).'*sysr.B).';
        
    es = norm(s0new-s0)/norm(s0); 
    fprintf(1,'so_irka iteration %03i convergence:\t %4.3e\n',kIter,es)
    s0 = s0new;
end

end