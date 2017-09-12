classdef testStabsep < sssTest
% testStabsep - testing of ssRed.stabsep
%
% Description:
%   The functionality of the stabsep is tested in the following way:
%    + Generation of random system with stable and unstable eigenvalues
%    + Decomposition using stabsep and verification of the resulting model
%
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> morlab@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Last Change:  02 Aug 2017
% Copyright (c) 2017 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------     
    
    methods(Test)
        function mainFunctionality(testCase)
            
            % Create random ssRed-object          
            nS  = ceil(10*rand);
            nAs = ceil(20*rand);
            nMs = ceil(2*rand);
            iS  = ceil(5*rand);
            oS  = ceil(5*rand);
            A   = diag([-100*rand(1,nS), zeros(1,nMs), 100*rand(1,nAs)]);
            B   = ones(nS+nAs+nMs,iS);
            C   = ones(oS,nS+nAs+nMs);
            
            sys = ssRed(A,B,C);
            
            % Stabsep
            [sysS, sysAs, Vs,Ws,Vas,Was] = stabsep(sys);
            
            % Define verfication tolerance
            tol = 1e-10;
            
            % Verify sys = sysS+sysAs
            [m, p,w]    = bode(sys);
            [m2,p2]     = bode(sysS+sysAs,w);
            verifyEqual(testCase,m,m2,'AbsTol',tol)
            verifyEqual(testCase,p,p2,'AbsTol',tol)
            
            % Verify real part of eigenvalues
            verifyTrue(testCase,all(real(eig(sysS))<0))
            verifyTrue(testCase,all(real(eig(sysAs))>=0))
            verifyEqual(testCase,length(find(real(eig(sysAs))==0)),nMs)
            
            % Verify dimensions
            verifyEqual(testCase,sysS.n,nS)
            verifyEqual(testCase,sysAs.n,nAs+nMs)
            
            verifySize(testCase,Vs,[sys.n,nS])
            verifySize(testCase,Ws,[sys.n,nS])
            
            verifySize(testCase,Vas,[sys.n,nAs+nMs])
            verifySize(testCase,Was,[sys.n,nAs+nMs])      
            
            % real orthogonal bases
            verifyTrue(testCase,isreal(Vs))
            verifyTrue(testCase,isreal(Ws))
            verifyTrue(testCase,isreal(Vas))
            verifyTrue(testCase,isreal(Was))
            
            verifyLessThanOrEqual(testCase,norm(Vs'*Vs-eye(nS)),tol)
            verifyLessThanOrEqual(testCase,norm(Ws'*Ws-eye(nS)),tol)
            verifyLessThanOrEqual(testCase,norm(Vas'*Vas-eye(nAs+nMs)),tol)
            verifyLessThanOrEqual(testCase,norm(Was'*Was-eye(nAs+nMs)),tol)
        end
    end
end

