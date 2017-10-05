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
% clear all
% clc
% clearvars -global
% 
% sys = sss('CDPlayer');
% Opts.purpose = 'MOR'; Opts.rctol = 1e-3;
% [sysr, data] = crksm(sys, zeros(1,10), Opts);
% % [sysr, data] = crksm(sys, [], zeros(1,10), Opts);
% bode(sys,'-',sysr,'--r');

%TODO: komplexe Shifts nutzen
methods(Test)
%% SISO
    function SISO_Vsided(testCase) 
        SISObenchmarksSysCell = getSISObenchmarks();
        for i=1:length(SISObenchmarksSysCell)
             %  get system
             sys = SISObenchmarksSysCell{i};
             %  get interpolation data
             r = 10; %s0_inp = 100*randn(1,r); 
             s0_inp = -eigs(sys,r).';
             %  test system with crksm
             sysrCrksm = crksm(sys, s0_inp);
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end
        
    function SISO_Wsided(testCase) 
        SISObenchmarksSysCell = getSISObenchmarks();
        for i=1:length(SISObenchmarksSysCell)
             %  get system
             sys = SISObenchmarksSysCell{i};
             %  get interpolation data
             r = 10; s0_out = 100*randn(1,r); %s0_out = -eigs(sys,r).';
             %  test system with crksm
             sysrCrksm = crksm(sys, [], s0_out);
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    function SISO_TwoSided(testCase) 
        SISObenchmarksSysCell = getSISObenchmarks();
        for i=1:length(SISObenchmarksSysCell)
             %  get system
             sys = SISObenchmarksSysCell{i};
             %  get interpolation data
             r = 10; s0_inp = 100*randn(1,r); s0_out = 100*randn(1,r); %s0_out = -eigs(sys,r).';
             %  test system with crksm
             sysrCrksm = crksm(sys, s0_inp, s0_out);
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    function SISO_TwoSidedHermite(testCase) 
        SISObenchmarksSysCell = getSISObenchmarks();
        for i=1:length(SISObenchmarksSysCell)
             %  get system
             sys = SISObenchmarksSysCell{i};
             %  get interpolation data
             r = 10; s0_inp = zeros(1,r); s0_out = zeros(1,r); %s0_out = -eigs(sys,r).';
             %  test system with crksm
             sysrCrksm = crksm(sys, s0_inp, s0_out);
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

%% MIMO
%     function MIMO_Block_Vsided(testCase) 
%         MIMObenchmarksSysCell = getMIMObenchmarks();
%     end
% 
%     function MIMO_Block_Wsided(testCase)
%         MIMObenchmarksSysCell = getMIMObenchmarks();
%     end
% 
%     function MIMO_Block_TwoSided(testCase)
%         MIMObenchmarksSysCell = getMIMObenchmarks();
%     end
% 
%     function MIMO_Block_TwoSidedHermite(testCase)
%         MIMObenchmarksSysCell = getMIMObenchmarks();
%     end
% 
%     function MIMO_Tangential_Vsided(testCase) 
%         MIMObenchmarksSysCell = getMIMObenchmarks();
%     end
% 
%     function MIMO_Tangential_Wsided(testCase)
%         MIMObenchmarksSysCell = getMIMObenchmarks();
%     end
% 
%     function MIMO_Tangential_TwoSided(testCase)
%         MIMObenchmarksSysCell = getMIMObenchmarks();
%     end
% 
%     function MIMO_Tangential_TwoSidedHermite(testCase)
%         MIMObenchmarksSysCell = getMIMObenchmarks();
%     end

%% MISO_SIMO
%     function MISO_SIMO_Block_Vsided(testCase)
%         MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks();
%     end
% 
%     function MISO_SIMO_Block_Wsided(testCase)
%         MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks();
%     end
% 
%     function MISO_SIMO_Block_TwoSided(testCase)
%         MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks();
%     end
% 
%     %MISO_SIMO_Block_TwoSidedHermite not supported in crksm!!
% 
%     function MISO_SIMO_Tangential_Vsided(testCase)
%         MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks();
%     end
% 
%     function MISO_SIMO_Tangential_Wsided(testCase)
%         MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks();
%     end
% 
%     function MISO_SIMO_Tangential_TwoSided(testCase)
%         MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks();
%     end
% 
%     function MISO_SIMO_Tangential_TwoSidedHermite(testCase)
%         MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks();
%     end
end  
 
end

%% ***************************** AUXILIARY ********************************
function SISObenchmarksSysCell = getSISObenchmarks()

%     SISOBenchmarks = {'beam.mat','building.mat','eady.mat','fom.mat','heat-cont.mat',...
%         'random','SpiralInductorPeec.mat'};
    
    SISOBenchmarks = {'building.mat','fom.mat','heat-cont.mat'};
 
    nFiles = length(SISOBenchmarks);
    
    SISObenchmarksSysCell=cell(1,nFiles);
    for i=1:nFiles
        sys = sss(SISOBenchmarks{i});
        SISObenchmarksSysCell{i}=sys;
    end
 
end
 
function MIMObenchmarksSysCell = getMIMObenchmarks()

    MIMOBenchmarks = {'CDplayer','iss'};
 
    nFiles = length(MIMOBenchmarks);
    
    MIMObenchmarksSysCell=cell(1,nFiles);
    for i=1:nFiles
        sys = sss(MIMOBenchmarks{i});
        MIMObenchmarksSysCell{i}=sys;
    end
 
end

function MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks()

    MISO_SIMOBenchmarks = {'rail_1357','rail_5177','gyro'}; %'rail_79841'
 
    nFiles = length(MISO_SIMOBenchmarks);
    
    MISO_SIMObenchmarksSysCell=cell(1,nFiles);
    for i=1:nFiles
        sys = sss(MISO_SIMOBenchmarks{i});
        MISO_SIMObenchmarksSysCell{i}=sys;
    end
 
end
 