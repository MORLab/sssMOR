classdef testL2norm < sssTest
% testL2norm - testing of ssRed.l2norm
%
% Description:
%   The functionality of the l2norm is tested in the following way:
%    + Generation of random system with stable and unstable eigenvalues
%    + computation of the l2 norm and verification through numerical
%    quadrature
%    + Reduction with IRKA to a stable model and comparison between l2norm
%    and norm
%
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Last Change:  02 Aug 2017
% Copyright (c) 2017 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------     
    
    methods(Test)
        function unstableSISOmodel(testCase)
            
            % Create random ssRed-object          
            nS  = ceil(10*rand);
            nAs = ceil(20*rand);
            iS  = 1;
            oS  = 1;
            A   = diag([-100*rand(1,nS), 100*rand(1,nAs)]);
            B   = ones(nS+nAs,iS);
            C   = ones(oS,nS+nAs);
            
            sys = ssRed(A,B,C);
            
            % Define verfication tolerance
            tol = 1e-3;
            
            % l2norm
            n   = l2norm(sys);

            % Verify through numerical integration
            w = logspace(-5,5,1e5);
            m = squeeze(bode(sys,w));

            nInt = sqrt(1/pi*trapz(w,m.^2));

            verifyEqual(testCase,n,nInt,'AbsTol',tol,'RelTol',tol)
        end
        function stableModel(testCase)
            
            sys     = sss('building');
            sysr    = ssRed(1,1,1);
            while ~isstable(sysr)
                sysr    = irka(sys,ceil(20*rand));
            end
            
            % Define verfication tolerance
            tol = 1e-3;
            
            % l2norm
            L2n   = l2norm(sysr);
            H2n   = norm(sysr);
            verifyEqual(testCase,L2n,H2n,'AbsTol',tol,'RelTol',tol)
        end
    end
end

