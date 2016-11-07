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
            sys = loadSss('build.mat');
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
            sys = loadSss('build.mat');
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
            sys = loadSss('build.mat');
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
            sys = loadSss('build.mat');
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
            sys = loadSss('build.mat');
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
        
        function testSssFunctionality(testCase)
            % create a ssRed-model
            sys1 = loadSss('building');
            sys2 = loadSss('CDplayer');
            sys3 = loadSss('beam');
            sysr1 = ssRed('userDefined',[],sys1.A,sys1.B,sys1.C,sys1.D,sys1.E);
            sysr2 = ssRed('userDefined',[],sys2.A,sys2.B,sys2.C,sys2.D,sys2.E);
            sysr3 = ssRed('userDefined',[],sys3.A,sys3.B,sys3.C,sys3.D,sys3.E);
            
            % test append
            sys_appended = append(sys1,sys2);
            sysr_appended = append(sysr1,sysr2);
            verification(testCase,sysr_appended,sys_appended, 1e-8, ...
                        'Test for the append-function failed!');
            
            % test bode
            [mag, phase, omega] = bode(sys1);
            [magr, phaser] = bode(sysr1,omega);
            verifyEqual(testCase,mag,magr,'AbsTol',1e-5, ... 
                        'Test for the bode-function (magnitude) failed!');
            verifyEqual(testCase,phase,phaser,'AbsTol',1e-5, ... 
                        'Test for the bode-function (phase) failed!');
                    
            % test bodemag
            bodemag(sysr1,omega);
            
            % test bodeplot
            bodeplot(sysr1,omega);
            
            % test c2d
            % (can not be tested against the sss-function, because the 
            % selectable methods do not match)
            sysrd = c2d(sysr1,0.2);
            
            % test clear (does not exist for ss-objects)
            
            % test connect
            sysr1.u = 'w'; sysr1.y = 'u';
            sysr3.u = 'u'; sysr3.y = 'y';
            sys1.u = {'w'}; sys1.y = {'u'};
            sys3.u = {'u'}; sys3.y = {'y'};
            sys_connected = connect(sys1,sys3,{'w'},{'y'});
            sysr_connected = connect(sysr1,sysr3,'w','y');
            verification(testCase,sysr_connected,sys_connected, 1e-8, ...
                        'Test for the append-function failed!');
                    
            % test connectSss (does not exist for ss-objects)
            
            % test dcgain
            gain = dcgain(sys3);
            gainr = dcgain(sysr3);
            verifyEqual(testCase,gain,gainr,'AbsTol',1e-5, ... 
                        'Test for the dcgain-function failed!');
                    
            % test decayTime
            decTime = decayTime(sys1);
            decTimer = decayTime(sysr1);
            verifyEqual(testCase,decTime,decTimer,'AbsTol',1e-5, ... 
                        'Test for the dcgain-function failed!');
                    
            % test diag (does not exist for ss-objects)
            
            % test disp (does something different for ss-objects)
            
            % test eig
            [V,D,W] = eig(sys1);
            [Vr,Dr,Wr] = eig(sysr1);
            verifyEqual(testCase,V,Vr,'AbsTol',1e-5, ... 
                        'Test for the eig-function failed (matrix V)!');
            verifyEqual(testCase,D,Dr,'AbsTol',1e-5, ... 
                        'Test for the eig-function failed (matrix D)!');
            verifyEqual(testCase,W,Wr,'AbsTol',1e-5, ... 
                        'Test for the eig-function failed (matrix W)!');
                    
            % test eigs
            [V,D,flag] = eigs(sys1);
            [Vr,Dr,flagr] = eigs(sysr1);
            verifyEqual(testCase,V,Vr,'AbsTol',1, ... 
                        'Test for the eigs-function failed (matrix V)!');
            verifyEqual(testCase,D,Dr,'AbsTol',1, ... 
                        'Test for the eigs-function failed (matrix D)!');
                    
            % test freqresp (frd-object generation does not work for ssRed)
            [G, w] = freqresp(sys1);
            Gr = freqresp(sysr1,w);
            verifyEqual(testCase,G,Gr,'AbsTol',1e-5, ... 
                        'Test for the freqresp-function failed!');
                    
            % test impulse
            Opts.tf = 1;
            [tf,h,t] = impulse(sys1,Opts);
            [tfr,hr,tr] = impulse(sysr1,t,Opts);
            verifyEqual(testCase,h,hr,'AbsTol',1e-5, ... 
                        'Test for the impulse-function failed!');
                    
            % test issd
            [i,nA] = issd(sys1);
            [ir,nAr] = issd(sysr1);
            verifyEqual(testCase,nA,nAr,'AbsTol',1e-5, ... 
                        'Test for the issd-function failed!');
                    
            % test isstable
            stab = isstable(sys1);
            stabr = isstable(sysr1);
            verifyEqual(testCase,+stab,+stabr,'Test for the isstable-function failed!');
            
            % test lsim
            u = rand(100,1);
            [Y,~,t] = lsim(sys1,u,0.01,'RK4',0.01);
            [Yr,tr,~] = lsim(sysr1,u,t);
            verifyEqual(testCase,Y,Yr','AbsTol',1e-5, ... 
                        'Test for the lsim-function failed!');
            
                   
                    
                    
            
        end
      
        
    end
    
end

function [] = verification (testCase, sys, sys_sss, absTol,message)
          % test if a ssRed-model and a sss-Model are equal
          verifyEqual(testCase,sys.A,full(sys_sss.A),'AbsTol', absTol, message);
          verifyEqual(testCase,sys.B,full(sys_sss.B),'AbsTol', absTol, message);
          verifyEqual(testCase,sys.C,full(sys_sss.C),'AbsTol', absTol, message);
          verifyEqual(testCase,sys.D,full(sys_sss.D),'AbsTol', absTol, message);
          %verifyEqual(testCase,sys.E,full(sys_sss.E),'AbsTol', 1e-8, message);
end

