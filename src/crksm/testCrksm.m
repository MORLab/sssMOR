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
% Last Change:  10 Oct 2017
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
             Opts.strategy = 'ADI';
             Opts.method = 'heur';
             Opts.nShifts = 10;
             s0_inp = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.restolLyap = 1e-6;
             Opts.strategy = 'adaptive';
             Opts.stopCrit = 'residualLyap';
             Opts.shifts = 'dynamical';
             sysrCrksm = crksm(sys, s0_inp, Opts);
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end
        
    function SISO_Wsided(testCase) 
        Opts.restolLyap = 1e-7;
        Opts.shifts = 'dynamical';
        
        SISObenchmarksSysCell = getSISObenchmarks();
        for i=1:length(SISObenchmarksSysCell)
             %  get system
             sys = SISObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'ADI';
             Opts.method = 'heur';
             Opts.nShifts = 10;
             s0_out = initializeShifts(sys,10,1,Opts);
             
             %  test system with crksm
             Opts.restolLyap = 1e-6;
             Opts.strategy = 'adaptive';
             Opts.stopCrit = 'residualLyap';
             Opts.shifts = 'dynamical';
             sysrCrksm = crksm(sys, [], s0_out, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

%     function SISO_TwoSided(testCase) 
%         Opts.rctol = 1e-3;
% %         Opts.shifts = 'dynamical';
%         
%         SISObenchmarksSysCell = getSISObenchmarks();
%         for i=1:length(SISObenchmarksSysCell)
%              %  get system
%              sys = SISObenchmarksSysCell{i};
%              %  get interpolation data
%              r = 10; s0_inp = 100*randn(1,r); s0_out = 100*randn(1,r); %s0_out = -eigs(sys,r).';
%              %  test system with crksm
%              sysrCrksm = crksm(sys, s0_inp, s0_out, Opts);
%              
%              % evaluation
% %              NORM{i} = norm(sys - sysrCrksm);
%         end
%     end

    function SISO_TwoSidedHermite(testCase) 
        Opts.rctol = 1e-3;
        Opts.shifts = 'dynamical';
        
        SISObenchmarksSysCell = getSISObenchmarks();
        for i=1:length(SISObenchmarksSysCell)
             %  get system
             sys = SISObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'ADI';
             Opts.method = 'heur';
             Opts.nShifts = 10;
             [s0_inp,~,s0_out] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.restolLyap = 1e-6;
             Opts.strategy = 'adaptive';
             Opts.stopCrit = 'residualLyap';
             Opts.shifts = 'dynamical';
             sysrCrksm = crksm(sys, s0_inp, s0_out, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

%% MIMO Block
    function MIMO_Block_Vsided(testCase) 
        
        MIMObenchmarksSysCell = getMIMObenchmarks();
        for i=1:length(MIMObenchmarksSysCell)
             %  get system
             sys = MIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'ADI';
             Opts.method = 'heur';
             Opts.nShifts = 10;
             [s0_inp] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.restolLyap = 1e-6;
             Opts.strategy = 'adaptive';
             Opts.stopCrit = 'residualLyap';
             Opts.shifts = 'dynamical';
             sysrCrksm = crksm(sys, s0_inp, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    function MIMO_Block_Wsided(testCase)        
        MIMObenchmarksSysCell = getMIMObenchmarks();
        for i=1:length(MIMObenchmarksSysCell)
             %  get system
             sys = MIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'ADI';
             Opts.method = 'heur';
             Opts.nShifts = 10;
             [~,~,s0_out] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.restolLyap = 1e-6;
             Opts.strategy = 'adaptive';
             Opts.stopCrit = 'residualLyap';
             Opts.shifts = 'dynamical';
             sysrCrksm = crksm(sys, [], s0_out, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

%     function MIMO_Block_TwoSided(testCase)       
%         Opts.rctol = 1e-3;
% %         Opts.shifts = 'dynamical';
%         
%         MIMObenchmarksSysCell = getMIMObenchmarks();
%         for i=1:length(MIMObenchmarksSysCell)
%              %  get system
%              sys = MIMObenchmarksSysCell{i};
%              %  get interpolation data
%              r = 10; s0_inp = 100*randn(1,r); s0_out = 100*randn(1,r); %s0_out = -eigs(sys,r).';
%              %  test system with crksm
%              sysrCrksm = crksm(sys, s0_inp, s0_out, Opts);
%              
%              % evaluation
% %              NORM{i} = norm(sys - sysrCrksm);
%         end
%     end

%     function MIMO_Block_TwoSidedHermite(testCase)        
%         Opts.rctol = 1e-3;
% %         Opts.shifts = 'dynamical';
%         
%         MIMObenchmarksSysCell = getMIMObenchmarks();
%         for i=1:length(MIMObenchmarksSysCell)
%              %  get system
%              sys = MIMObenchmarksSysCell{i};
%              %  get interpolation data
%              r = 10; s0_inp = zeros(1,r); s0_out = zeros(1,r); %s0_out = -eigs(sys,r).';
%              %  test system with crksm
%              sysrCrksm = crksm(sys, s0_inp, s0_out, Opts);
%              
%              % evaluation
% %              NORM{i} = norm(sys - sysrCrksm);
%         end
%     end
    
%% MIMO Tangential
% 
%     function MIMO_Tangential_Vsided(testCase) 
%         MIMObenchmarksSysCell = getMIMObenchmarks();
% [s0,Rt,s0,Lt] = initiali,....
%     
% norm(V_rk - V_crksm); subspace(V_rk,V_crksm); norm(sysr_Rk - sysr_Crksm); 

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

%% MISO_SIMO Block
    function MISO_SIMO_Block_Vsided(testCase)
        
        MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks();
        for i=1:length(MISO_SIMObenchmarksSysCell)
             %  get system
             sys = MISO_SIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'ADI';
             Opts.method = 'heur';
             Opts.nShifts = 10;
             s0_inp = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.restolLyap = 1e-6;
             Opts.strategy = 'adaptive';
             Opts.stopCrit = 'residualLyap';
             Opts.shifts = 'dynamical';
             sysrCrksm = crksm(sys, s0_inp, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    function MISO_SIMO_Block_Wsided(testCase)
        
        MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks();
        for i=1:length(MISO_SIMObenchmarksSysCell)
            
             %  get system
             sys = MISO_SIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'ADI';
             Opts.method = 'heur';
             Opts.nShifts = 10;
             [~,~,s0_out] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.restolLyap = 1e-6;
             Opts.strategy = 'adaptive';
             Opts.stopCrit = 'residualLyap';
             Opts.shifts = 'dynamical';
             Opts.maxiter = 20;
             sysrCrksm = crksm(sys, [], s0_out, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    function MISO_SIMO_Block_TwoSided(testCase)
        MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks();
        for i=1:length(MISO_SIMObenchmarksSysCell)
             %  get system
             sys = MISO_SIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'ADI';
             Opts.method = 'heur';
             Opts.nShifts = 10;
             [s0_inp,~,s0_out] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.restolLyap = 1e-6;
             Opts.strategy = 'adaptive';
             Opts.stopCrit = 'residualLyap';
             Opts.shifts = 'dynamical';
             Opts.maxiter = 20;
             sysrCrksm = crksm(sys, s0_inp, s0_out, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    %MISO_SIMO_Block_TwoSidedHermite not supported in crksm!!


%% MISO_SIMO Tangential
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

    MISO_SIMOBenchmarks = {'rail_1357','rail_5177'}; %'rail_79841'  'gyro'
 
    nFiles = length(MISO_SIMOBenchmarks);
    
    MISO_SIMObenchmarksSysCell=cell(1,nFiles);
    for i=1:nFiles
        sys = sss(MISO_SIMOBenchmarks{i});
        MISO_SIMObenchmarksSysCell{i}=sys;
    end
 
end
 