classdef testSsRed < sssTest
% testSsRed - testing of ssRed.m
%
% Description:
%   The functionality of the ssRed-class is tested in the following way:
%    + Reduction of a SISO and of MIMO system to create a ssRed-object
%    + For the resulting ssRed-object some functions of the Control System
%      Toolbox are called. Currently this functions are
%       - tf
%       - bode
%       - impulse
%       - zero
%       - pole
%       - zpk
%       - tzero (only for MIMO-models)
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Niklas Kochdumper
% Last Change:  21 Aug 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------     
    
    methods(Test)
        function testSsRedTbr(testCase)
            
            % Create ssRed-object througth reduction with tbr (SISO)          
            sys = loadSss('building.mat');
            sysr = tbr(sys,10);
            
            % Test some functions of the Control System Toolbox on the
            % ssRed-object (SISO)
            tf(sysr);
            bode(sysr);
            impulse(sysr);
            pole(sysr);
            zero(sysr);
            zpk(sysr);
            
            % Try to reduce the ssRed-object again (SISO)
            sysr = tbr(sysr,5);
            
            % Create ssRed-object througth reduction with tbr (MIMO)          
            sys = loadSss('CDplayer.mat');
            sysr = tbr(sys,10);
            
            % Test some functions of the Control System Toolbox on the
            % ssRed-object (MIMO)
            tf(sysr);
            bode(sysr);
            impulse(sysr);
            pole(sysr);
            zero(sysr);
            zpk(sysr);
            tzero(sysr);
            
            % Try to reduce the ssRed-object again (MIMO)
            sysr = tbr(sysr,5);
        end
        
        function testSsRedIrka(testCase)
            
            warning('off','all')
            
            % Create ssRed-object througth reduction with irka (SISO)          
            sys = loadSss('building.mat');
            s0 = [0 1+i 1-i];
            sysr = irka(sys,s0);
            
            % Test some functions of the Control System Toolbox on the
            % ssRed-object (SISO)
            tf(sysr);
            bode(sysr);
            impulse(sysr);
            pole(sysr);
            zero(sysr);
            zpk(sysr);
            
            % Try to reduce the ssRed-object again (SISO)
            s0 = [0 1];
            sysr = irka(sysr,s0);
            
            % Create ssRed-object througth reduction with irka (MIMO)          
            sys = loadSss('CDplayer.mat');
            s0 = [0 1+i 1-i];
            L = [1 0 0; 0 1 1];
            sysr = irka(sys,s0,L,L);
            
            % Test some functions of the Control System Toolbox on the
            % ssRed-object (MIMO)
            tf(sysr);
            bode(sysr);
            impulse(sysr);
            pole(sysr);
            zero(sysr);
            zpk(sysr);
            tzero(sysr);
            
            % Try to reduce the ssRed-object again (MIMO)
            s0 = [0 1];
            L = [1 0; 0 1];
            sysr = irka(sysr,s0,L,L);
        end
        
        function testSsRedRk(testCase)
            
            warning('off','all')
            
            % Create ssRed-object througth reduction with rk (SISO)          
            sys = loadSss('building.mat');
            s0 = [0 1+i 1-i];
            sysr = rk(sys,s0);
            
            % Test some functions of the Control System Toolbox on the
            % ssRed-object (SISO)
            tf(sysr);
            bode(sysr);
            impulse(sysr);
            pole(sysr);
            zero(sysr);
            zpk(sysr);
            
            % Try to reduce the ssRed-object again (SISO)
            s0 = [0 1];
            sysr = rk(sysr,s0);
            
            % Create ssRed-object througth reduction with rk (MIMO)          
            sys = loadSss('CDplayer.mat');
            s0 = [0 1+i 1-i];
            L = [1 0 0; 0 1 1];
            sysr = rk(sys,s0,L,L);
            
            % Test some functions of the Control System Toolbox on the
            % ssRed-object (MIMO)
            tf(sysr);
            bode(sysr);
            impulse(sysr);
            pole(sysr);
            zero(sysr);
            zpk(sysr);
            tzero(sysr);
            
            % Try to reduce the ssRed-object again (MIMO)
            s0 = [0 1];
            L = [1 0; 0 1];
            sysr = rk(sysr,s0,L,L);
        end
        
        function testSsRedModalMor(testCase)
            
            warning('off','all')
            
            % Create ssRed-object througth reduction with modalMor (SISO)          
            sys = loadSss('building.mat');
            Opts.type = 'SR';
            Opts.real = 'real';
            sysr = modalMor(sys,10,Opts);
            
            % Test some functions of the Control System Toolbox on the
            % ssRed-object (SISO)
            tf(sysr);
            bode(sysr);
            impulse(sysr);
            pole(sysr);
            zero(sysr);
            zpk(sysr);
            
            % Try to reduce the ssRed-object again (SISO)
            sysr = modalMor(sysr,5,Opts);
            
            % Create ssRed-object througth reduction with modalMor (MIMO)          
            sys = loadSss('CDplayer.mat');
            Opts.type = 'SR';
            Opts.real = 'real';
            sysr = modalMor(sys,10,Opts);
            
            % Test some functions of the Control System Toolbox on the
            % ssRed-object (MIMO)
            tf(sysr);
            bode(sysr);
            impulse(sysr);
            pole(sysr);
            zero(sysr);
            zpk(sysr);
            tzero(sysr);
            
            % Try to reduce the ssRed-object again (MIMO)
            sysr = modalMor(sysr,5,Opts);
        end
        
        function testSsRedCure(testCase)
            
            warning('off','all')
            
            % Create ssRed-object througth reduction with cure_spark (SISO)          
            sys = loadSss('building.mat');
            Opts.cure = struct('nk',4, 'redfun', 'spark', 'verbose', 1, 'stopval',12);
            sysr = cure(sys,Opts);
            
            % Test some functions of the Control System Toolbox on the
            % ssRed-object (SISO)
            tf(sysr);
            bode(sysr);
            impulse(sysr);
            pole(sysr);
            zero(sysr);
            zpk(sysr);
            tzero(sysr);
            
            % Try to reduce the ssRed-object again (SISO)
            sysr = cure(sysr,Opts);
            
            % Create ssRed-object througth reduction with cure_irka (SISO)          
            Opts.cure = struct('nk',4, 'redfun', 'irka', 'verbose', 1, 'stopval',12);
            sysr = cure(sys,Opts);
            
            % Test some functions of the Control System Toolbox on the
            % ssRed-object (SISO)
            tf(sysr);
            bode(sysr);
            impulse(sysr);
            pole(sysr);
            zero(sysr);
            zpk(sysr);
            tzero(sysr);
            
            % Try to reduce the ssRed-object again (SISO)
            sysr = cure(sysr,Opts);
            
            % Create ssRed-object througth reduction with cure_rk+pork (SISO)          
            Opts.cure = struct('nk',4, 'redfun', 'rk+pork', 'verbose', 1, 'stopval',12);
            sysr = cure(sys,Opts);
            
            % Test some functions of the Control System Toolbox on the
            % ssRed-object (SISO)
            tf(sysr);
            bode(sysr);
            impulse(sysr);
            pole(sysr);
            zero(sysr);
            zpk(sysr);
            tzero(sysr);
            
            % Try to reduce the ssRed-object again (SISO)
            sysr = cure(sysr,Opts);
            
            % Create ssRed-object througth reduction with cure_rk+pork (MIMO)          
            sys = loadSss('CDplayer.mat');
            Opts.cure = struct('nk',4, 'redfun', 'rk+pork', 'verbose', 1, 'stopval',12);
            sysr = cure(sys,Opts);
            
            % Test some functions of the Control System Toolbox on the
            % ssRed-object (MIMO)
            tf(sysr);
            bode(sysr);
            impulse(sysr);
            pole(sysr);
            zero(sysr);
            zpk(sysr);
            tzero(sysr);
        end
        
    end
    
end

