classdef testCrksm < sssTest
% testCrksm - testing of crksm.m
%
% Description:
%   The function crksm.m is tested for Lyapunov and MOR purposes for both
%   SISO and MIMO models
%
%   For SISO: use preferably 'ADI'+'heur' for initializeShifts 
%   (other options: 'ROM' or 'eigs')
%
%   For MIMO-Block: use preferably 'ADI'+'heur' for initializeShifts 
%   (other options: 'ROM' or 'eigs')
%
%   For MIMO-Tangential: only use 'ROM' or 'eigs' for initializeShifts 
%   (other options DO NOT compute right and left tangential directions Rt and Lt)
%
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> morlab@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Maria Cruz Varona, Paul Heidenreich
% Last Change:  17 Nov 2017
% Copyright (c) 2017 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------ 

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
             Opts.adiShiftsMethod = 'heur';
             Opts.nShifts = 10;
             s0_inp = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.stopCrit = 'residualLyap'; Opts.restolLyap = 1e-6;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'adaptive';
             [sysrCrksm,V,W,S,data] = crksm(sys, s0_inp, Opts);
             
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
             Opts.adiShiftsMethod = 'heur';
             Opts.nShifts = 10;
             s0_out = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.stopCrit = 'residualLyap'; Opts.restolLyap = 1e-6;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'adaptive';
             [sysrCrksm,V,W,R,data] = crksm(sys, [], s0_out, Opts);
             
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
             Opts.adiShiftsMethod = 'heur';
             Opts.nShifts = 10;
             s0_inp = initializeShifts(sys,Opts.nShifts,1,Opts);
             Opts.strategy = 'eigs';
             [~,~,s0_out] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.purpose = 'MOR';
             Opts.restolMOR = 1e-3;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'adaptive'; %[{'adaptive'} / 'eigs']
             [sysrCrksm,V,W,Z,data] = crksm(sys, s0_inp, s0_out, Opts); % error with building, eady, fom, SpiralInductorPeec
             
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
             Opts.adiShiftsMethod = 'heur';
             Opts.nShifts = 10;
             s0_inp = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm 
             Opts.purpose = 'MOR';
             Opts.restolMOR = 1e-3;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'adaptive'; %[{'adaptive'} / 'eigs']
             [sysrCrksm,V,W,Z,data] = crksm(sys, s0_inp, s0_inp, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

%% MIMO Block (m=p)
    function MIMO_Block_Vsided(testCase) 
    % Test ok
        MIMObenchmarksSysCell = getMIMObenchmarks(); % 'CDplayer', 'iss'
        
        for i=1:length(MIMObenchmarksSysCell)
             % get system
             sys = MIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'ADI';
             Opts.adiShiftsMethod = 'heur';
             Opts.nShifts = 10;
             [s0_inp] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.stopCrit = 'residualLyap'; Opts.restolLyap = 1e-6;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'adaptive';
             [sysrCrksm,V,W,S,data] = crksm(sys, s0_inp, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    function MIMO_Block_Wsided(testCase)
    % Test ok
        MIMObenchmarksSysCell = getMIMObenchmarks(); % 'CDplayer', 'iss'
         
        for i=1:length(MIMObenchmarksSysCell)
             % get system
             sys = MIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'ADI';
             Opts.adiShiftsMethod = 'heur';
             Opts.nShifts = 10;
             [~,~,s0_out] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.stopCrit = 'residualLyap'; Opts.restolLyap = 1e-6;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'adaptive';
             [sysrCrksm,V,W,R,data] = crksm(sys, [], s0_out, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    function MIMO_Block_TwoSided(testCase)       
    % Test ok, ends in warning, unstable 
        MIMObenchmarksSysCell = getMIMObenchmarks(); % 'CDplayer', 'iss'
        
        for i=1:length(MIMObenchmarksSysCell)
             % get system
             sys = MIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'ADI';
             Opts.adiShiftsMethod = 'heur';
             Opts.nShifts = 10;
             s0_inp = initializeShifts(sys,Opts.nShifts,1,Opts);
             Opts.strategy = 'eigs';
             [~,~,s0_out] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.purpose = 'MOR';
             Opts.restolMOR = 1e-3;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'adaptive'; %[{'adaptive'} / 'eigs']
             [sysrCrksm,V,W,Z,data] = crksm(sys, s0_inp, s0_out, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    function MIMO_Block_TwoSidedHermite(testCase) 
    % Test ok, iss ends in warning, unstable
        
        MIMObenchmarksSysCell = getMIMObenchmarks(); % 'CDplayer', 'iss'
        
        for i=1:length(MIMObenchmarksSysCell)
             %  get system
             sys = MIMObenchmarksSysCell{i};
             %  get interpolation data
             Opts.strategy = 'ROM';
             %Opts.adiShiftsMethod = 'heur';
             Opts.nShifts = 10;
             s0_inp = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.purpose = 'MOR';
             Opts.restolMOR = 1e-3;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'eigs'; %[{'adaptive'} / 'eigs']
             [sysrCrksm,V,W,Z,data] = crksm(sys, s0_inp, s0_inp, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end
    
%% MIMO Tangential (m=p)
% diese fälle noch analysieren
% 
    function MIMO_Tangential_Vsided(testCase) 
        % test ok
        MIMObenchmarksSysCell = getMIMObenchmarks(); % 'CDplayer', 'iss'
        
        for i=1:length(MIMObenchmarksSysCell)
             %  get system
             sys = MIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'eigs';
             Opts.nShifts = 10;
             [s0_inp,Rt] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.stopCrit = 'residualLyap'; Opts.restolLyap = 1e-6;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'adaptive';
             [sysrCrksm,V,W,S,data] = crksm(sys, s0_inp,Rt,Opts);
             %norm(V_rk - V_crksm); subspace(V_rk,V_crksm); norm(sysr_Rk - sysr_Crksm); 
        end
    end

    function MIMO_Tangential_Wsided(testCase)
    % test ok
        MIMObenchmarksSysCell = getMIMObenchmarks(); % 'CDplayer', 'iss'
        
        for i=1:length(MIMObenchmarksSysCell)
             %  get system
             sys = MIMObenchmarksSysCell{i};
             
             %  get interpolation data
%              Opts.strategy = 'eigs';
%              Opts.nShifts = 10;
%              [s0_inp,Rt,s0_out,Lt] = initializeShifts(sys,Opts.nShifts,1,Opts);
             Opts.strategy = 'ROM';
             Opts.nShifts = 10;
             [s0_inp,Rt,s0_out,Lt] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.stopCrit = 'residualLyap'; Opts.restolLyap = 1e-6;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'adaptive';
             [sysrCrksm,V,W,R,data] = crksm(sys,[],s0_out,[],Lt,Opts);
             %norm(V_rk - V_crksm); subspace(V_rk,V_crksm); norm(sysr_Rk - sysr_Crksm); 
        end
    end

    function MIMO_Tangential_TwoSided(testCase)
        % test ok
        MIMObenchmarksSysCell = getMIMObenchmarks(); % 'CDplayer', 'iss'
        
        for i=1:length(MIMObenchmarksSysCell)
             %  get system
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
             Opts.restolMOR = 1e-3;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'eigs'; %[{'adaptive'} / 'eigs']
             [sysrCrksm,V,W,Z,data] = crksm(sys,s0_inp,s0_out,Rt,Lt,Opts);
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
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'eigs'; %[{'adaptive'} / 'eigs']
             [sysrCrksm,V,W,Z,data] = crksm(sys,s0_inp,s0_inp,Rt,Lt,Opts);
             %norm(V_rk - V_crksm); subspace(V_rk,V_crksm); norm(sysr_Rk - sysr_Crksm); 
         end
    end

%% MISO_SIMO Block (m~=p)
    function MISO_SIMO_Block_Vsided(testCase)
        % Test ok
        MISO_SIMObenchmarksSysCell = getMISO_SIMObenchmarks();

        for i=1:length(MISO_SIMObenchmarksSysCell)
             %  get system
             sys = MISO_SIMObenchmarksSysCell{i};
             
             %  get interpolation data
             Opts.strategy = 'ADI';
             Opts.adiShiftsMethod = 'heur';
             Opts.nShifts = 10;
             s0_inp = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.stopCrit = 'residualLyap'; Opts.restolLyap = 1e-6;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'adaptive';
             [sysrCrksm,V,W,S,data] = crksm(sys, s0_inp, Opts);
             
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
             Opts.adiShiftsMethod = 'heur';
             Opts.nShifts = 10;
             [~,~,s0_out] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.maxiter = 100;
             Opts.stopCrit = 'residualLyap'; Opts.restolLyap = 1e-6;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'adaptive';
             [sysrCrksm,V,W,R,data] = crksm(sys, [], s0_out, Opts);
             
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
             Opts.adiShiftsMethod = 'heur';
             Opts.nShifts = 10;
             s0_inp = initializeShifts(sys,Opts.nShifts,1,Opts);
             Opts.strategy = 'eigs';
             [~,~,s0_out] = initializeShifts(sys,Opts.nShifts,1,Opts);
             
             %  test system with crksm
             Opts.maxiter = 20;
             Opts.stopCrit = 'residualLyap'; Opts.restolLyap = 1e-6;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'adaptive'; %[{'adaptive'} / 'eigs']
             [sysrCrksm,V,W,Z,data] = crksm(sys, s0_inp, s0_out, Opts);
             
             % evaluation
%              NORM{i} = norm(sys - sysrCrksm);
        end
    end

    %MISO_SIMO_Block_TwoSidedHermite not supported in crksm!!


%% MISO_SIMO Tangential (m~=p)
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
             Opts.stopCrit = 'residualLyap'; Opts.restolLyap = 1e-6;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'adaptive';
             [sysrCrksm,V,W,S,data] = crksm(sys, s0_inp,Rt,Opts);
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
             Opts.stopCrit = 'residualLyap'; Opts.restolLyap = 1e-6;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'adaptive';
             [sysrCrksm,V,W,R,data] = crksm(sys,[],s0_out,[],Lt,Opts);
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
             Opts.restolMOR = 1e-3;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'eigs'; %[{'adaptive'} / 'eigs']
             [sysrCrksm,V,W,Z,data] = crksm(sys,s0_inp,s0_out,Rt,Lt,Opts);
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
             Opts.purpose = 'MOR';
             Opts.maxiter = 300;
             Opts.restolMOR = 1e-3;
             Opts.shifts = 'dynamical'; %[{'dynamical'} / 'fixedCyclic']
             Opts.strategy = 'eigs'; %[{'adaptive'} / 'eigs']
             [sysrCrksm,V,W,Z,data] = crksm(sys,s0_inp,s0_inp,Rt,Lt,Opts);
             %norm(V_rk - V_crksm); subspace(V_rk,V_crksm); norm(sysr_Rk - sysr_Crksm); 
             
        end
    end
end  
 
end

%% ***************************** AUXILIARY ********************************
function SISObenchmarksSysCell = getSISObenchmarks()

%     SISOBenchmarks = {'building.mat','random','heat-cont.mat','beam.mat','eady.mat',...
%         'fom.mat','SpiralInductorPeec.mat'};
    
    SISOBenchmarks = {'building.mat','fom.mat','heat-cont.mat'};
%     SISOBenchmarks = {'heat-cont.mat','beam.mat','eady.mat','fom.mat'};
 
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
 