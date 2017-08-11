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
%    + Conversion from sss, ss and dss to ssRed
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> morlab@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Niklas Kochdumper, Alessandro Castagnotto
% Last Change:  02 Aug 2017
% Copyright (c) 2016,2017 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------     
    
    methods(Test)
        function testSsRedTbr(testCase)
            
            % Create ssRed-object througth reduction with tbr (SISO)          
            sys = sss('building.mat');
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
            sys = sss('CDplayer.mat');
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
            sys = sss('building.mat');
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
            sys = sss('CDplayer.mat');
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
            sys = sss('building.mat');
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
            sys = sss('CDplayer.mat');
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
            sys = sss('building.mat');
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
            sys = sss('CDplayer.mat');
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
            sys = sss('building.mat');
            sysr = cure(sys);
            
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
            sysr = cure(sysr);
            
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
            Opts.cure = struct('nk',4, 'redfun', 'rk+pork', 'verbose', 1,...
                'stop','nmax','stopval',12,'init','sm');
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
            sys = sss('CDplayer.mat');
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
            sys1 = sss('building');
            sys2 = sss('CDplayer');
            sys3 = sss('beam');
            sysr1 = ssRed(sys1.A,sys1.B,sys1.C,sys1.D,sys1.E);
            sysr2 = ssRed(sys2.A,sys2.B,sys2.C,sys2.D,sys2.E);
            sysr3 = ssRed(sys3.A,sys3.B,sys3.C,sys3.D,sys3.E);
            
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
            
            % test clear
            emptySys = clear(sysr1);
            verifyEqual(testCase,double(isempty(emptySys)),1,'AbsTol',1e-5, ... 
                        'Test for the clear-function failed!');
            
            % test connect
            sysr1.u = 'w'; sysr1.y = 'u';
            sysr3.u = 'u'; sysr3.y = 'y';
            sys1.u = {'w'}; sys1.y = {'u'};
            sys3.u = {'u'}; sys3.y = {'y'};
            sys_connected = connect(sys1,sys3,{'w'},{'y'});
            sysr_connected = connect(sysr1,sysr3,'w','y');
            verification(testCase,sysr_connected,sys_connected, 1e-8, ...
                        'Test for the append-function failed!');
                    
            % test connectSss (not implemented for ssRed-objects)
            
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
            sysd = diag(sys1);
            sysdr = diag(sysr1);
            verification(testCase,sysdr,sysd, 1e-8, ...
                        'Test for the diag-function failed!');
            
            % test disp
            infostr = disp(sysr1);
            
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
                    
            % test freqresp
            [G, w] = freqresp(sys1);
            Gr = freqresp(sysr1,w);
            verifyEqual(testCase,G,Gr,'AbsTol',1e-5, ... 
                        'Test for the freqresp-function failed!');
                 
            % test frd
            frdObj = frd(sys1,w);
            frdObjr = frd(sysr1,w);
            verifyEqual(testCase,frdObj.ResponseData,frdObjr.ResponseData,'AbsTol',1e-5, ... 
                        'Test for the freqresp-function failed (frd-object)!');
                    
            % test impulse
            Opts.tf = 1;
            [tf,h,t] = impulse(sys1,Opts);
            [tfr,hr,tr] = impulse(sysr1,t,Opts);
            verifyEqual(testCase,h,hr,'AbsTol',1e-3, ... 
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
            Opts.method = 'RK4';
            Opts.Ts_sample = 0.01;
            u = rand(100,1);
            t = 0:0.01:(length(u)-1)*0.01;
            [Y,~,~] = lsim(sys1,u,0.01,[],Opts);
            [Yr,~,~] = lsim(sysr1,u',t);
            verifyEqual(testCase,Y,Yr,'AbsTol',1e-5, ... 
                        'Test for the lsim-function failed!');
             
            % test lyapchol
            [R,L] = lyapchol(sys1);
            [Rr,Lr] = lyapchol(sysr1);
            verifyEqual(testCase,R,Rr,'AbsTol',1e-5, ... 
                        'Test for the lyapchol-function failed!');
            verifyEqual(testCase,L,Lr,'AbsTol',1e-5, ... 
                        'Test for the lyapchol-function failed!');
                    
            % test minus
            sys_diff    = sys1-sys3;
            sysr_diff   = sysr1-sysr3;
            verifyLessThan(testCase,norm(sys_diff-sysr_diff), 1e-12, ...
                        'Test for the minus-function failed!')
                    
            % test mtimes
            sys_prod    = mtimes(sys1,sys3);
            sysr_prod   = mtimes(sysr1,sysr3);
            verifyLessThan(testCase,norm(sys_prod-sysr_prod),1e-12, ...
                        'Test for the mtimes-function failed!')
                    
            % test norm
            nrm     = norm(sys1, 'inf');
            nrmr    = norm(sysr1,Inf);
            verifyEqual(testCase,nrm,nrmr,'AbsTol',1e-3, ... 
                        'Test for the norm-function failed (inf-norm)!');                    
            nrm     = norm(sys1, 2);
            nrmr    = norm(sysr1,2);
            verifyEqual(testCase,nrm,nrmr,'AbsTol',1e-3, ... 
                        'Test for the norm-function failed (2-norm)!');
                    
            % test plus
            sys_add     = sys1+sys3;
            sysr_add    = sysr1+sysr3;
            verifyLessThan(testCase,norm(sys_add-sysr_add), 1e-12, ...
                        'Test for the plus-function failed!');
                    
            % test pole
            k = 10;
            p = pole(sys1,k);
            p = sort(p,'descend');
            p = cplxpair(p);
            pr = pole(sysr1);
            pr = sort(pr,'descend');
            pr = cplxpair(pr(1:k));
            verifyEqual(testCase,p,pr,'AbsTol',1e-5, ... 
                        'Test for the poles-function failed!');
                    
            % test pzmap
            [p,z] = pzmap(sys1);
            [pr,zr] = pzmap(sysr1);
            pr = pr(1:size(p));
            zr = zr(1:size(z));
            verifyEqual(testCase,real(p),real(pr),'AbsTol',1e-5, ... 
                        'Test for the pzmap-function failed (poles)!');
            verifyEqual(testCase,real(z),real(zr),'AbsTol',1e-5, ... 
                        'Test for the pzmap-function failed (zeros)!');
                    
            % test residue
            [r,p,d] = residue(sys1);
            [rr,pr,dr] = residue(sysr1);
            verifyEqual(testCase,r,rr,'AbsTol',1, ... 
                        'Test for the residue-function failed (residuals)!');
            verifyEqual(testCase,p,pr,'AbsTol',1, ... 
                        'Test for the residue-function failed (eigenvalues)!');
            verifyEqual(testCase,full(d),dr,'AbsTol',1, ... 
                        'Test for the residue-function failed (feedthrough)!');
                    
            % test sigma
            [s, omega] = sigma(sys1);
            sr = sigma(sysr1,omega);
            verifyEqual(testCase,s,sr,'AbsTol',1e-8, ... 
                        'Test for the sigma-function failed!');
                   
            % test sim (not implemented for ssRed-objects)
            
            % test size
            p = size(sys1,1);
            pr = size(sysr1,1);
            m = size(sys1,2);
            mr = size(sysr1,2);
            verifyEqual(testCase,p,pr, ... 
                        'Test for the size-function failed (output dimension)!');
            verifyEqual(testCase,m,mr, ... 
                        'Test for the size-function failed (input dimension)!'); 
                    
            % test spy
            spy(sysr1,'test plot'); drawnow; close(gcf)
            
            % test ss (not implemented for ssRed-objects)
            
            % test sss (not implemented for ssRed-objects)
            
            % test step
            Opts.tf = 1;
            [tf,h,t] = step(sys1,Opts);
            [tfr,hr,tr] = step(sysr1,t,Opts);
            verifyEqual(testCase,h,hr,'AbsTol',1e-5, ... 
                        'Test for the step-function failed!');
                    
            % test truncate
            sys_sub = truncate(sys2,1,[1 2]);
            sysr_sub = truncate(sysr2,1,[1 2]);
            verification(testCase,sysr_sub,sys_sub, 1e-8, ...
                        'Test for the truncate-function failed!');
                    
            % test zero
            k = 10;
            z = zero(sys1,k);
            z = sort(z,'descend');
            z = cplxpair(z);
            zr = zero(sysr1);
            zr = sort(zr,'descend');
            zr = cplxpair(zr(1:k));
            verifyEqual(testCase,z,zr,'AbsTol',1e-5, ... 
                        'Test for the zeros-function failed!');
                    
            % test zpk
            k = 10;
            zpkData = zpk(sys1,k,'lm');
            zpkDatar = zpk(sysr1);
            z = cell2mat(zpkData.z);
            zr = cell2mat(zpkDatar.z);
            z = sort(z,'descend');
            zr = sort(zr,'descend');
            z = cplxpair(z);
            zr = cplxpair(zr(1:k));
            verifyEqual(testCase,z,zr,'AbsTol',1e-5, ... 
                        'Test for the zpk-function failed (invariant zeros)!');
            p = cell2mat(zpkData.p);
            pr = cell2mat(zpkDatar.p);
            p = sort(p,'descend');
            pr = sort(pr,'descend');
            p = cplxpair(p);
            pr = cplxpair(pr(1:k));
            verifyEqual(testCase,p,pr,'AbsTol',1e-5, ... 
                        'Test for the zpk-function failed (invariant zeros)!');     
            
            
        end
        
        function testConversionToSsRed(testCase)
            tol = 1e-6;
            
            % sss
            sys     = sss('building');
            sys2    = ssRed(sys); 
            
            verifyClass(testCase,sys2,'ssRed')
            verifyLessThanOrEqual(testCase,norm(sys-sys2),tol);
            
            % ss
            sys     = stabsep(rss(ceil(100*rand)));
            sys2    = ssRed(sys); 
            
            verifyClass(testCase,sys2,'ssRed')
            verifyLessThanOrEqual(testCase,norm(sys-sys2),tol);
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

