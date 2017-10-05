classdef testCrksm < sssTest
% testCrksm - testing of crksm.m
%
% Description:
%   The function crksm.m is tested on....:
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> morlab@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Maria Cruz Varona, Paul Heidenreich
% Last Change:  05 Oct 2017
% Copyright (c) 2017 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------ 
clear all
clc
clearvars -global

sys = sss('CDPlayer');
Opts.purpose = 'MOR'; Opts.rctol = 1e-3;
[sysr, data] = crksm(sys, zeros(1,10), Opts);
% [sysr, data] = crksm(sys, [], zeros(1,10), Opts);
bode(sys,'-',sysr,'--r');

 methods(Test)
        function SISO_Vsided(testCase) 
             
        end
        
        function SISO_Wsided(testCase) 
             
        end
        
        function SISO_TwoSided(testCase) 
             
        end
        
        function SISO_TwoSidedHermite(testCase) 
            
        end
        
        function MIMO_Block_Vsided(testCase) 
            
        end
        
        function MIMO_Block_Wsided(testCase)
            
        end
        
        function MIMO_Block_TwoSided(testCase)
            
        end
        
        function MIMO_Block_TwoSidedHermite(testCase)
            
        end
        
        function MIMO_Tangential_Vsided(testCase) 
            
        end
        
        function MIMO_Tangential_Wsided(testCase)
            
        end
        
        function MIMO_Tangential_TwoSided(testCase)
            
        end
        
        function MIMO_Tangential_TwoSidedHermite(testCase)
            
        end
        
        % MISO_SIMO_Block not supported in crksm!!
        
        function MISO_SIMO_Tangential_Vsided(testCase)
            
        end
        
        function MISO_SIMO_Tangential_Wsided(testCase)
            
        end
        
        function MISO_SIMO_Tangential_TwoSided(testCase)
            
        end
        
        function MISO_SIMO_Tangential_TwoSidedHermite(testCase)
            
        end
    end
end