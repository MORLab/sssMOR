classdef testProjectiveMor < sssTest
% testProjectiveMor - testing of projectiveMor.m
%
% Description:
%   The function projectiveMor is tested on:
%    + input/output usage
%    + Galerkin projection using real shift
%    + Galerkin projection using complex shift
%    + Petrov-Galerkin projection using real shift
%    + Petrov-Galerkin projection using complex shift
%
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> morlab@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto, Lisa Jeschek
% Last Change:  06 Apr 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------
    
    methods(Test)
        
         function IOusage(testCase) 
            %Correct parsing of input/outputs
            sys = testCase.sysCell{1};
            s0 = 1; r = ones(sys.m,length(s0));

            V = (sys.A-s0*sys.E)\(sys.B*r);
            sysr    = projectiveMor(sys,V);
            sysr2   = projectiveMor(sys,V,[]);
            sysr3   = projectiveMor(sys,V,struct('trans','T'));
            sysr4   = projectiveMor(sys,V,[],struct('trans','T')); 

            verifyEqual(testCase,{dssdata(sysr)},{dssdata(sysr2)},'AbsTol', 1e-10,'Wrong parsing')
            verifyEqual(testCase,{dssdata(sysr)},{dssdata(sysr3)},'AbsTol', 1e-10,'Wrong parsing')
            verifyEqual(testCase,{dssdata(sysr)},{dssdata(sysr4)},'AbsTol', 1e-10,'Wrong parsing')
         end 
         function GalerkinRealSingleT (testCase) 
              %Galerkin projection using real shifts
                for i=1:length(testCase.sysCell)

                    sys = testCase.sysCell{i};
                    s0 = 1; r = ones(sys.m,length(s0));
                
                    V = (sys.A-s0*sys.E)\(sys.B*r);
                    sysr = projectiveMor(sys,V);
                    
                    Me = moments(sys-sysr,s0,1)*r;
                    
                    verifyEqual(testCase,sysr.n,1,'Reduced order not correct')
                    verifyLessThanOrEqual(testCase,Me,1e-10,'Moments do not match')  
                    
                    verification(testCase,sysr)
                end
         end 
         function GalerkinRealSingleH (testCase) 
              %Galerkin projection using real shifts
                for i=1:length(testCase.sysCell)

                    sys = testCase.sysCell{i};
                    s0 = 1; l = ones(sys.p,length(s0));
                
                    W = (sys.A-s0*sys.E)'\(sys.C'*l);
                    sysr = projectiveMor(sys,W,struct('trans','H'));
                    
                    Me = l'*moments(sys-sysr,s0,1);
                    
                    verifyEqual(testCase,sysr.n,1,'Reduced order not correct')
                    verifyLessThanOrEqual(testCase,Me,1e-10,'Moments do not match')  
                    
                    verification(testCase,sysr)
                end
         end 
         function GalerkinComplexSingleT (testCase) 
                for i=1:length(testCase.sysCell)

                    sys = testCase.sysCell{i};
                    s0 = 1+1i; r = ones(sys.m,length(s0));
                
                    V = (sys.A-s0*sys.E)\(sys.B*r);
                    sysr = projectiveMor(sys,V);
                    
                    Me = moments(sys-sysr,s0,1)*r;
                    
                    verifyEqual(testCase,sysr.n,1,'Reduced order not correct')
                    verifyLessThanOrEqual(testCase,Me,1e-10,'Moments do not match')  
                    
                    verification(testCase,sysr)
                end
         end 
         function GalerkinComplexSingleH (testCase) 
              %Galerkin projection using real shifts
                for i=1:length(testCase.sysCell)

                    sys = testCase.sysCell{i};
                    s0 = 1+1i; l = ones(sys.p,length(s0));
                
                    W = (sys.A-s0*sys.E)'\(sys.C'*l);
                    sysr = projectiveMor(sys,W,struct('trans','H'));
                    
                    Me = l'*moments(sys-sysr,s0,1);
                    
                    verifyEqual(testCase,sysr.n,1,'Reduced order not correct')
                    verifyLessThanOrEqual(testCase,Me,1e-10,'Moments do not match')  
                    
                    verification(testCase,sysr)
                end
         end 
         function PetrovGalerkinComplexSingleT (testCase) 
              %Galerkin projection using real shifts
                for i=1:length(testCase.sysCell)

                    sys = testCase.sysCell{i};
                    s0 = 1+1i; 
                    r = ones(sys.m,length(s0));
                    l = ones(sys.p,length(s0));
                    
                    V = (sys.A-s0*sys.E)\(sys.B*r);
                    W = (sys.A-s0*sys.E).'\(sys.C.'*l);
                    sysr = projectiveMor(sys,V,W);
                    
                    MeL = l.'*moments(sys-sysr,s0,1);
                    MeR = moments(sys-sysr,s0,1)*r;
                    MeH = mmat(...
                            mmat(repmat(l.',[1,1,2*length(s0)]),...
                                    moments(sys-sysr,s0,2)),...
                               repmat(r,[1,1,2*length(s0)]));

                    verifyEqual(testCase,sysr.n,1,'Reduced order not correct')
                    verifyLessThanOrEqual(testCase,MeL,1e-10,'Moments do not match')  
                    verifyLessThanOrEqual(testCase,MeR,1e-10,'Moments do not match')  
                    verifyLessThanOrEqual(testCase,MeH,1e-10,'Moments do not match')  

                    
                    verification(testCase,sysr)
                end
         end 
         function PetrovGalerkinComplexSingleH (testCase) 
              %Galerkin projection using real shifts
                for i=1:length(testCase.sysCell)

                    sys = testCase.sysCell{i};
                    s0 = 1+1i; 
                    r = ones(sys.m,length(s0));
                    l = ones(sys.p,length(s0));
                    
                    V = (sys.A-s0*sys.E)\(sys.B*r);
                    W = (sys.A-s0*sys.E)'\(sys.C'*l);
                    sysr = projectiveMor(sys,V,W,struct('trans','H'));
                    
                    MeL = l'*moments(sys-sysr,s0,1);
                    MeR = moments(sys-sysr,s0,1)*r;
                    MeH = mmat(...
                            mmat(repmat(l',[1,1,2*length(s0)]),...
                                    moments(sys-sysr,s0,2)),...
                               repmat(r,[1,1,2*length(s0)]));

                    verifyEqual(testCase,sysr.n,1,'Reduced order not correct')
                    verifyLessThanOrEqual(testCase,MeL,1e-10,'Moments do not match')  
                    verifyLessThanOrEqual(testCase,MeR,1e-10,'Moments do not match')  
                    verifyLessThanOrEqual(testCase,MeH,1e-10,'Moments do not match')  

                    
                    verification(testCase,sysr)
                end
         end 
    end  
end

function verification(testCase, sysr)
    %additonal verifications for each ROM
     
    verifyEqual(testCase, rank(full(sysr.A)), size(sysr.A,1),...
            'Rank(Ar) is not full');
       verifyEqual(testCase, rank(full(sysr.E)), size(sysr.A,1),...
            'Rank(Er) is not full');
       verifyEqual(testCase, nnz(isinf(sysr.A)), 0, ...
            'Ar contains Inf');
       verifyEqual(testCase, nnz(isinf(sysr.E)), 0, ...
            'Er contains Inf');
       verifyEqual(testCase, nnz(isnan(sysr.A)), 0, ...
            'Ar contains Nan');
       verifyEqual(testCase, nnz(isnan(sysr.E)), 0, ...
            'Er contains Nan');
end
         
         