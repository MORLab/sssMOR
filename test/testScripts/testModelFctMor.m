classdef testModelFctMor < sssTest
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
% Last Change:  22 Nov 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------ 
 
    methods(Test)
        function generalFunctionality(testCase) 
            warning('off','sssMOR:irka:maxiter')
            n = 2;
            s0 = zeros(1,n);
            
            redFct   = @(sys,s) irka(sys,s);
            redFctOut= @(sys,s) getDesiredOutput(redFct,[1,4],sys,s);
            
            for i=1:length(testCase.sysCell)
                %  test system
                sys=testCase.sysCell{i};  
                sys= sys(1,1);
                
                [sysr, ~, sysm, relH2err] = modelFctMor(sys,redFctOut,s0);  
                
                verifyClass(testCase,sysr,'ssRed')
                verifyTrue(testCase,isstable(sysm))
                verifyGreaterThanOrEqual(testCase,size(sysm.A,1),n);
                verifyEqual(testCase,norm(sysm-sysr)/norm(sysm),relH2err)           
            end
            warning('on','sssMOR:irka:maxiter')
        end
    end
end