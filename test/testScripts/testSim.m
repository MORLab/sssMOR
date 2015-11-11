classdef testSim < matlab.unittest.TestCase
% testSim - Testing of sim.m
%
% Description:
%   The function sim.m is tested on:
%    + Simulation of a SSS, SISO benchmark system (build).
%    + Simulation of a SSS, SISO random system.
%    + Simulation of a DSSS, SISO benchmark system (SpiralInductorPeec).
%    + Simulation of a DSS, MISO random system.
%    + Simulation of a DSS, SIMO random system.
%    + Simulation of a SSS, MIMO benchmark system (iss).
%    + Simulation of a DSSS, MIMO benchmark system (rail_1357).
%    + Simulation of a SSS, MIMO random system.
%    + Verification of the obtained results.
%
%------------------------------------------------------------------
% This file is part of sssMOR, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Institute
% of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit www.rt.mw.tum.de
% For any suggestions, submission and/or bug reports, mail us at
%                   -> sssMOR@rt.mw.tum.de <-
%------------------------------------------------------------------
% Authors:      Maria Cruz Varona
% Last Change:  7 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------
    methods(Test)
        function testSISObench(testCase)
            load('build.mat');
            sysSparse=sss(A,B,C);
            sys=ss(A,B,C,zeros(1,1));
            nOutputs = sysSparse.p;
            
%             isstable = isstable(sysSparse)
            
            Ts = 1e-4;
            t = 0:Ts:10;
            u = idinput(length(t),'rgs',[0 0.5/(1/2/Ts)])';
            datau = iddata([],u',Ts);
            
            dataSparseRK4 = sim(sysSparse,datau,'RK4');
            dataSparseFoEuler = sim(sysSparse,datau,'forwardEuler');
            dataSparseBaEuler = sim(sysSparse,datau,'backwardEuler');
            dataSparseDiscrete = sim(c2d(sysSparse,Ts),datau,'discrete');
            
            actSim = [dataSparseRK4.y,dataSparseFoEuler.y,dataSparseBaEuler.y,dataSparseDiscrete.y];
            expSim = lsim(sys,u,t);
            
%             figure; plot(t,actSim); hold on; plot(t,expSim); 
%             xlabel('Time [s]'); ylabel('Amplitude');
            
            for iCase=1:4
                verification(testCase, actSim(:,(iCase*nOutputs-nOutputs+1):iCase*nOutputs), expSim);
            end
        end
        function testSISOrandom(testCase)
            sys=rss(35);
            sysSparse=sss(sys);
            nOutputs = sysSparse.p;
            
%             isstable = isstable(sysSparse)
            
            Ts = 1e-4;
            t = 0:Ts:10;
            u = sin(2*pi*1*t);
            datau = iddata([],u',Ts);
            
            dataSparseRK4 = sim(sysSparse,datau,'RK4');
            dataSparseFoEuler = sim(sysSparse,datau,'forwardEuler');
            dataSparseBaEuler = sim(sysSparse,datau,'backwardEuler');
            dataSparseDiscrete = sim(c2d(sysSparse,Ts),datau,'discrete');
            
            actSim = [dataSparseRK4.y,dataSparseFoEuler.y,dataSparseBaEuler.y,dataSparseDiscrete.y];
            expSim = lsim(sys,u,t);
            
            figure; plot(t,actSim); hold on; plot(t,expSim); 
            xlabel('Time [s]'); ylabel('Amplitude');
            
            for iCase=1:4
                verification(testCase, actSim(:,(iCase*nOutputs-nOutputs+1):iCase*nOutputs), expSim);
            end
        end
        function testDSSSISObench(testCase)
            load('SpiralInductorPeec.mat');
            sysSparse=sss(A,B,C,[],E);
            sys=ss(sysSparse);
            nOutputs = sysSparse.p;
            
%             isstable = isstable(sysSparse)
            
            Ts = 1e-10;
            t = 0:Ts:3e-7;
%             u = idinput(length(t),'rgs',[0 0.5/(1/2/Ts)])';
            u = double(t>=Ts);
            datau = iddata([],u',Ts);
            
            dataSparseRK4 = sim(sysSparse,datau,'RK4');
            dataSparseFoEuler = sim(sysSparse,datau,'forwardEuler');
            dataSparseBaEuler = sim(sysSparse,datau,'backwardEuler');
            dataSparseDiscrete = sim(c2d(sysSparse,Ts),datau,'discrete');
            
            actSim = [dataSparseRK4.y,dataSparseFoEuler.y,dataSparseBaEuler.y,dataSparseDiscrete.y];
            expSim = lsim(sys,u,t);
            
            figure; plot(t,actSim); hold on; plot(t,expSim); 
            xlabel('Time [s]'); ylabel('Amplitude');
            
            for iCase=1:4
                verification(testCase, actSim(:,(iCase*nOutputs-nOutputs+1):iCase*nOutputs), expSim);
            end
        end
        function testDSSSMISOrandom(testCase)
            n=35;
            nInputs=5;
            sys=rss(n);
            sys=dss(sys.A,rand(n,nInputs),sys.C,rand(1,nInputs),rand(size(sys.A)));
            sysSparse=sss(sys);
            nOutputs = sysSparse.p;
            
%             isstable = isstable(sysSparse)

            Ts = 1e-4;
            t = 0:Ts:10;
            u = idinput(length(t),'rgs',[0 0.5/(1/2/Ts)])';
                        
            U = repmat(u,nInputs,1);
            dataU = iddata([],U',Ts);
            
            dataSparseRK4 = sim(sysSparse,dataU,'RK4');
            dataSparseFoEuler = sim(sysSparse,dataU,'forwardEuler');
            dataSparseBaEuler = sim(sysSparse,dataU,'backwardEuler');
            dataSparseDiscrete = sim(c2d(sysSparse,Ts),dataU,'discrete');
            
            actSim = [dataSparseRK4.y,dataSparseFoEuler.y,dataSparseBaEuler.y,dataSparseDiscrete.y];
            expSim = lsim(sys,U,t);
            
%             figure; plot(t,actSim); hold on; plot(t,expSim); 
%             xlabel('Time [s]'); ylabel('Amplitude');

            for iCase=1:4
                verification(testCase, actSim(:,(iCase*nOutputs-nOutputs+1):iCase*nOutputs), expSim);
            end
        end
        function testDSSSSIMOrandom(testCase)
            n=35;
            nOutputs=5;
            sys=rss(n);
            sys=dss(sys.A,sys.B,rand(nOutputs,n),rand(nOutputs,1),rand(size(sys.A)));
            sysSparse=sss(sys);
            nOutputs = sysSparse.p;
            
%             isstable = isstable(sysSparse)
            
            Ts = 1e-4;
            t = 0:Ts:10;
            u = idinput(length(t),'rgs',[0 0.5/(1/2/Ts)])';
            datau = iddata([],u',Ts);
            
            dataSparseRK4 = sim(sysSparse,datau,'RK4');
            dataSparseFoEuler = sim(sysSparse,datau,'forwardEuler');
            dataSparseBaEuler = sim(sysSparse,datau,'backwardEuler');
            dataSparseDiscrete = sim(c2d(sysSparse,Ts),datau,'discrete');
            
            actSim = [dataSparseRK4.y,dataSparseFoEuler.y,dataSparseBaEuler.y,dataSparseDiscrete.y];
            expSim = lsim(sys,u,t);
            
%             figure; plot(t,actSim); hold on; plot(t,expSim); 
%             xlabel('Time [s]'); ylabel('Amplitude');

            for iCase=1:4
                verification(testCase, actSim(:,(iCase*nOutputs-nOutputs+1):iCase*nOutputs), expSim);
            end
        end
%         function testMIMObench(testCase)
%             load('cdplayer.mat');
%             sysSparse=sss(A,B,C);
%             sys=ss(full(A),full(B),full(C),zeros(2,2));
%             nOutputs = sysSparse.p;
%             
%             isstable = isstable(sysSparse)
% 
%             Ts = 1e-4;
%             t = 0:Ts:1;
%             u = idinput(length(t),'rgs',[0 0.5/(1/2/Ts)])';
% %             u = double(t>=Ts);
%                         
%             U = repmat(u,size(B,2),1);
%             dataU = iddata([],U',Ts);
%             
%             dataSparseRK4 = sim(sysSparse,dataU,'RK4');
%             dataSparseFoEuler = sim(sysSparse,dataU,'forwardEuler');
%             dataSparseBaEuler = sim(sysSparse,dataU,'backwardEuler');
%             dataSparseDiscrete = sim(c2d(sysSparse,Ts),dataU,'discrete');
%             
%             actSim = [dataSparseRK4.y,dataSparseFoEuler.y,dataSparseBaEuler.y,dataSparseDiscrete.y];
%             expSim = lsim(sys,U,t);
%             
%             figure; plot(t,actSim); hold on; plot(t,expSim); 
%             xlabel('Time [s]'); ylabel('Amplitude');
% 
%             for iCase=1:4
%                 verification(testCase, actSim(:,(iCase*nOutputs-nOutputs+1):iCase*nOutputs), expSim);
%             end
%         end

        function testMIMObench(testCase)
            load('rail_1357.mat');
            sysSparse=sss(A,B,C,[],E);
            sys=dss(full(A),full(B),full(C),zeros(size(C,1),size(B,2)),full(E));
            nOutputs = sysSparse.p;
            
%             isstable = isstable(sysSparse)

            Ts = 1e-4;
            t = 0:Ts:2;
%             u = idinput(length(t),'rgs',[0 0.5/(1/2/Ts)])';
            u = double(t>=Ts);
                        
            U = repmat(u,size(B,2),1);
            dataU = iddata([],U',Ts);
            
            dataSparseRK4 = sim(sysSparse,dataU,'RK4');
            dataSparseFoEuler = sim(sysSparse,dataU,'forwardEuler');
            dataSparseBaEuler = sim(sysSparse,dataU,'backwardEuler');
            dataSparseDiscrete = sim(c2d(sysSparse,Ts),dataU,'discrete');
            
            actSim = [dataSparseRK4.y,dataSparseFoEuler.y,dataSparseBaEuler.y,dataSparseDiscrete.y];
            expSim = lsim(sys,U,t);
            
            figure; plot(t,actSim); hold on; plot(t,expSim); 
            xlabel('Time [s]'); ylabel('Amplitude');

            for iCase=1:4
                verification(testCase, actSim(:,(iCase*nOutputs-nOutputs+1):iCase*nOutputs), expSim);
            end
        end
        function testDSSSMIMOrandom(testCase)
            n=35;
            nInputs=7;
            nOutputs=5;
            sys=rss(n);
            sys=dss(sys.A,rand(n,nInputs),rand(nOutputs,n),zeros(nOutputs,nInputs),rand(size(sys.A)));
            sysSparse=sss(sys);
            
%             isstable = isstable(sysSparse)

            Ts = 1e-4;
            t = 0:Ts:3;
            u = double(t>=Ts);
                        
            U = repmat(u,nInputs,1);
            dataU = iddata([],U',Ts);
            
            dataSparseRK4 = sim(sysSparse,dataU,'RK4');
            dataSparseFoEuler = sim(sysSparse,dataU,'forwardEuler');
            dataSparseBaEuler = sim(sysSparse,dataU,'backwardEuler');
            dataSparseDiscrete = sim(c2d(sysSparse,Ts),dataU,'discrete');
            
            actSim = [dataSparseRK4.y,dataSparseFoEuler.y,dataSparseBaEuler.y,dataSparseDiscrete.y];
            expSim = lsim(sys,U,t);
            
            figure; plot(t,actSim); hold on; plot(t,expSim); 
            xlabel('Time [s]'); ylabel('Amplitude');
            
            for iCase=1:4
                verification(testCase, actSim(:,(iCase*nOutputs-nOutputs+1):iCase*nOutputs), expSim);
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1e-3,...
    'Difference between actual and expected exceeds relative tolerance');
end
