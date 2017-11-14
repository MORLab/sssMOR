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
    % test ok
        SISObenchmarksSysCell = getSISObenchmarks();
        % building, fom, heat-cont
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
             %sysrCrksm = crksm(sys, s0_inp, Opts);
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end
        
    function SISO_Wsided(testCase) 
    % test ok
    
        SISObenchmarksSysCell = getSISObenchmarks();
        % building, fom, heat-cont
        
        for i=1:length(SISObenchmarksSysCell)
             %  get system
             sys = SISObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'ADI';
             Opts.method = 'heur';
             Opts.nShifts = 10;
             s0_out = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.restolLyap = 1e-6;
             Opts.strategy = 'adaptive';
             Opts.stopCrit = 'residualLyap';
             Opts.shifts = 'dynamical';
             %sysrCrksm = crksm(sys, [], s0_out, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    function SISO_TwoSided(testCase)
    % test ok   
        SISObenchmarksSysCell = getSISObenchmarks();
        % building, fom, heat-cont
        
        for i=1:length(SISObenchmarksSysCell)
             %  get system
             sys = SISObenchmarksSysCell{i};
             %  get interpolation data
             Opts.strategy = 'ADI';
             Opts.method = 'heur';
             Opts.nShifts = 10;
             s0_inp = initializeShifts(sys,Opts.nShifts,1,Opts);
             Opts.strategy = 'eigs';
             [~,~,s0_out] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.purpose = 'MOR';
             Opts.strategy = 'adaptive';
             Opts.restolMOR = 1e-3;
             %sysrCrksm = crksm(sys, s0_inp, s0_out, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    function SISO_TwoSidedHermite(testCase) 
    % Test ok
    
        SISObenchmarksSysCell = getSISObenchmarks();
         % building, fom, heat-cont
         
        for i=1:length(SISObenchmarksSysCell)
             %  get system
             sys = SISObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'ADI';
             Opts.method = 'heur';
             Opts.nShifts = 10;
             s0_inp = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm 
             Opts.purpose = 'MOR';
             Opts.strategy = 'adaptive';
             Opts.restolMOR = 1e-3;
             %sysrCrksm = crksm(sys, s0_inp, s0_inp, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

%% MIMO Block
    function MIMO_Block_Vsided(testCase) 
    % Test ok
        
        MIMObenchmarksSysCell = getMIMObenchmarks();
        % 'CDplayer','iss'
        
        for i=1:length(MIMObenchmarksSysCell)
             % get system
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
             %sysrCrksm = crksm(sys, s0_inp, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    function MIMO_Block_Wsided(testCase)
    % Test ok
        MIMObenchmarksSysCell = getMIMObenchmarks();
        % 'CDplayer','iss'
         
        for i=1:length(MIMObenchmarksSysCell)
             % get system
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
             %sysrCrksm = crksm(sys, [], s0_out, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    function MIMO_Block_TwoSided(testCase)       
    % Test ok, endet in warning, instabil 
        
        MIMObenchmarksSysCell = getMIMObenchmarks();
        % 'CDplayer','iss'
        for i=1:length(MIMObenchmarksSysCell)
             % get system
             sys = MIMObenchmarksSysCell{i};
             %  get interpolation data
             Opts.strategy = 'ADI';
             Opts.method = 'heur';
             Opts.nShifts = 10;
             s0_inp = initializeShifts(sys,Opts.nShifts,1,Opts);
             Opts.strategy = 'eigs';
             [~,~,s0_out] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.strategy = 'adaptive';
             Opts.purpose = 'MOR';
             Opts.restolMOR = 1e-3;
             %sysrCrksm = crksm(sys, s0_inp, s0_out, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    function MIMO_Block_TwoSidedHermite(testCase) 
    % Test ok, iss endet in warning, ist instabil
        
        MIMObenchmarksSysCell = getMIMObenchmarks();
        % 'CDplayer','iss'
        for i=1:length(MIMObenchmarksSysCell)
             %  get system
             sys = MIMObenchmarksSysCell{i};
             %  get interpolation data
             Opts.strategy = 'ROM';
             %Opts.method = 'heur';
             Opts.nShifts = 10;
             s0_inp = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.strategy = 'adaptive';
             Opts.purpose = 'MOR';
             Opts.restolMOR = 1e-3;
             %sysrCrksm = crksm(sys, s0_inp, s0_inp, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end
    
%% MIMO Tangential
% diese fälle noch analysieren
% 
    function MIMO_Tangential_Vsided(testCase) 
        % test ok
        MIMObenchmarksSysCell = getMIMObenchmarks();
        
        for i=1:length(MIMObenchmarksSysCell)
             %  get system
             sys = MIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'eigs';
             Opts.nShifts = 10;
             [s0_inp,Rt] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.restolLyap = 1e-6;
             Opts.strategy = 'adaptive';
             Opts.stopCrit = 'residualLyap';
             Opts.shifts = 'dynamical';
             %sysrCrksm = crksm(sys, s0_inp,Rt,Opts);
             %norm(V_rk - V_crksm); subspace(V_rk,V_crksm); norm(sysr_Rk - sysr_Crksm); 
        end
    end

    function MIMO_Tangential_Wsided(testCase)
    % test ok
        MIMObenchmarksSysCell = getMIMObenchmarks();
        
        for i=1:length(MIMObenchmarksSysCell)
             %  get system: 'CDplayer','iss'
             sys = MIMObenchmarksSysCell{i};
             
             %  get interpolation data
%              Opts.strategy = 'eigs';
%              Opts.nShifts = 10;
%              [s0_inp,Rt,s0_out,Lt] = initializeShifts(sys,Opts.nShifts,1,Opts);
             Opts.strategy = 'ROM';
             Opts.nShifts = 10;
             [s0_inp,Rt,s0_out,Lt] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.restolLyap = 1e-6;
             Opts.strategy = 'adaptive';
             Opts.stopCrit = 'residualLyap';
             Opts.shifts = 'dynamical';
             %sysrCrksm = crksm(sys,[],s0_out,[],Lt,Opts);
             %norm(V_rk - V_crksm); subspace(V_rk,V_crksm); norm(sysr_Rk - sysr_Crksm); 
        end
    end

    function MIMO_Tangential_TwoSided(testCase)
        % test ok
        MIMObenchmarksSysCell = getMIMObenchmarks();
        
        for i=1:length(MIMObenchmarksSysCell)
             %  get system: 'CDplayer','iss'
             sys = MIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'ROM';
             Opts.nShifts = 10;
             [s0_inp,Rt,~,~] = initializeShifts(sys,Opts.nShifts,1,Opts);
             Opts.strategy = 'eigs';
             Opts.nShifts = 10;
             [~,~,s0_out,Lt] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.purpose = 'MOR';
             Opts.strategy = 'eigs';
             Opts.restolMOR = 1e-3;
             %sysrCrksm = crksm(sys,s0_inp,s0_out,Rt,Lt,Opts);
             %norm(V_rk - V_crksm); subspace(V_rk,V_crksm); norm(sysr_Rk - sysr_Crksm); 
             
        end
    end

    function MIMO_Tangential_TwoSidedHermite(testCase)
        % test ok
        MIMObenchmarksSysCell = getMIMObenchmarks();
        
         for i=1:length(MIMObenchmarksSysCell)
             %  get system
             sys = MIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'eigs';
             Opts.nShifts = 10;
             [s0_inp,Rt,~,Lt] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.purpose = 'MOR';
             Opts.restolMOR = 1e-3;
              Opts.strategy = 'eigs';
             %sysrCrksm = crksm(sys,s0_inp,s0_inp,Rt,Lt,Opts);
             %norm(V_rk - V_crksm); subspace(V_rk,V_crksm); norm(sysr_Rk - sysr_Crksm); 
         end
    end

%% MISO_SIMO Block
    function MISO_SIMO_Block_Vsided(testCase)
        
        MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks();
        % Test ok
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
             %sysrCrksm = crksm(sys, s0_inp, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    function MISO_SIMO_Block_Wsided(testCase)
        % Test ok
        
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
             Opts.maxiter = 100;
             %sysrCrksm = crksm(sys, [], s0_out, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    function MISO_SIMO_Block_TwoSided(testCase)
        MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks();
        % fall geht nicht, weil m~=p und auffuellen Dimensionsprobleme
        % macht
        for i=1:length(MISO_SIMObenchmarksSysCell)
             %  get system
             sys = MISO_SIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'ADI';
             Opts.method = 'heur';
             Opts.nShifts = 10;
             s0_inp = initializeShifts(sys,Opts.nShifts,1,Opts);
             Opts.strategy = 'eigs';
             [~,~,s0_out] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.restolLyap = 1e-6;
             Opts.strategy = 'adaptive';
             Opts.stopCrit = 'residualLyap';
             Opts.shifts = 'dynamical';
             Opts.maxiter = 20;
             %sysrCrksm = crksm(sys, s0_inp, s0_out, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    %MISO_SIMO_Block_TwoSidedHermite not supported in crksm!!


%% MISO_SIMO Tangential
    function MISO_SIMO_Tangential_Vsided(testCase)
    % test ok
        MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks();
        
        for i=1:length(MISO_SIMObenchmarksSysCell)
             %  get system
             sys = MISO_SIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'eigs';
             Opts.nShifts = 10;
             [s0_inp,Rt] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.restolLyap = 1e-6;
             Opts.strategy = 'adaptive';
             Opts.stopCrit = 'residualLyap';
             Opts.shifts = 'dynamical';
             %sysrCrksm = crksm(sys, s0_inp,Rt,Opts);
             %norm(V_rk - V_crksm); subspace(V_rk,V_crksm); norm(sysr_Rk - sysr_Crksm); 
        end
    end
 
    function MISO_SIMO_Tangential_Wsided(testCase)
    % test ok
        MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks();
        
        for i=1:length(MISO_SIMObenchmarksSysCell)
             %  get system
             sys = MISO_SIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'eigs';
             Opts.nShifts = 10;
             [~,~,s0_out,Lt] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.restolLyap = 1e-6;
             Opts.strategy = 'adaptive';
             Opts.stopCrit = 'residualLyap';
             Opts.shifts = 'dynamical';
             %sysrCrksm = crksm(sys,[],s0_out,[],Lt,Opts);
             %norm(V_rk - V_crksm); subspace(V_rk,V_crksm); norm(sysr_Rk - sysr_Crksm); 
        end
    end
 
    function MISO_SIMO_Tangential_TwoSided(testCase)
    % test ok
        MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks();
        
        for i=1:length(MISO_SIMObenchmarksSysCell)
             %  get system
             sys = MISO_SIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'ROM';
             Opts.nShifts = 10;
             [s0_inp,Rt,s0_out,Lt] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.purpose = 'MOR';
             Opts.maxiter = 300;
             Opts.strategy = 'eigs';
             Opts.restolMOR = 1e-3;
             %sysrCrksm = crksm(sys,s0_inp,s0_out,Rt,Lt,Opts);
             %norm(V_rk - V_crksm); subspace(V_rk,V_crksm); norm(sysr_Rk - sysr_Crksm); 
             
        end
    end
 
    function MISO_SIMO_Tangential_TwoSidedHermite(testCase)
    % test ok
        MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks();
        
        for i=1:length(MISO_SIMObenchmarksSysCell)
             %  get system
             sys = MISO_SIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'eigs';
             Opts.nShifts = 10;
             [s0_inp,Rt,~,Lt] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.maxiter = 300;
             Opts.purpose = 'MOR';
             Opts.restolMOR = 1e-3;
             Opts.strategy = 'eigs';
             %sysrCrksm = crksm(sys,s0_inp,s0_inp,Rt,Lt,Opts);
             %norm(V_rk - V_crksm); subspace(V_rk,V_crksm); norm(sysr_Rk - sysr_Crksm); 
             
        end
    end
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
 