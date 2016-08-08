classdef testTbr < sssTest
% testtbr - testing of tbr.m
%
% Description:
%   The function tbr.m is tested (5 tests) on:
%    + q: 5, 10, 15, 20, 25
%    + test systems: building, beam, fom, random, LF10 (with E-matrix).
%    + comparing sysr with directly calculated solution 
%    + W'*V identity matrix
%    + Ar, Er purely real
%    + rank(Ar), rank(Er) full
%    + Neither Inf nor NaN in Ar, Er
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
%               Lisa Jeschek
% Last Change:  11 Feb 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------ 
    
    methods(Test)
        function testTbr1(testCase) %q=5
            load('building.mat');      
            q=5;
            
            [sysr, V, W, calchsv] = tbr(sss(A,B,C,0),q);
            load building hsv;
            actSolution={full(sysr.A),full(sysr.B),full(sysr.C),V,W, calchsv};
           
            S=lyapchol(full(A),full(B));
            R=lyapchol(full((A)'),full(C)');
            [Usvd,sigma,Vsvd] = svd(S*R');
            
            expV=S'*Usvd(:,1:q)*(sigma(1:q,1:q))^(-1/2);
            expW=((sigma(1:q,1:q))^(-1/2)*Vsvd(:,1:q)'*R)'; 
            
            expSolution={expW'*full(A)*expV, expW'*full(B),full(C)*expV,  expV, expW, hsv};
            verification(testCase, actSolution, expSolution, sysr);
            verifyLessThan(testCase, norm((W'*V)-eye(size(V,2))),1.01,...
               'W*V not identity matrix');
        end
        
        function testTbr2(testCase) %q=15
            load('beam.mat');
            q=15;
            
            [sysr, V, W] = tbr(sss(A,B,C,0),q);
            actSolution={full(sysr.A),full(sysr.B),full(sysr.C),V,W};

            S=lyapchol(full(A),full(B));
            R=lyapchol(full((A)'),full(C)');
            [Usvd,sigma,Vsvd] = svd(S*R');
            
            expV=S'*Usvd(:,1:q)*(sigma(1:q,1:q))^(-1/2);
            expW=((sigma(1:q,1:q))^(-1/2)*Vsvd(:,1:q)'*R)'; 
            
            expSolution={expW'*full(A)*expV, expW'*full(B),full(C)*expV,  expV, expW};
            verification(testCase, actSolution, expSolution, sysr);
            verifyLessThan(testCase, norm((W'*V)-eye(size(V,2))),1.01,...
               'W*V not identity matrix');
        end        
        function testTbr3(testCase) %q=25
            load('fom.mat');
            q=25;
            Opts.type='tbr';
            
            [sysr, V, W] = tbr(sss(A,B,C,0),q, Opts);
            actSolution={full(sysr.A),full(sysr.B),full(sysr.C),V,W};
            
            S=lyapchol(full(A),full(B));
            R=lyapchol(full((A)'),full(C)');
            [Usvd,sigma,Vsvd] = svd(S*R');
            
            expV=S'*Usvd(:,1:q)*(sigma(1:q,1:q))^(-1/2);
            expW=((sigma(1:q,1:q))^(-1/2)*Vsvd(:,1:q)'*R)'; 
            
            expSolution={expW'*full(A)*expV, expW'*full(B),full(C)*expV,  expV, expW};
            verification(testCase, actSolution, expSolution, sysr);
            verifyLessThan(testCase, norm((W'*V)-eye(size(V,2))),1.01,...
                'W*V not identity matrix');
        end
        function testTbr4(testCase) %q=20
            load('random.mat');
            q=20;
            
            [sysr, V, W] = tbr(sss(A,B,C,0),q);
            actSolution={full(sysr.A),full(sysr.B),full(sysr.C),V,W};
            
            S=lyapchol(full(A),full(B));
            R=lyapchol(full((A)'),full(C)');
            [Usvd,sigma,Vsvd] = svd(S*R');
            
            expV=S'*Usvd(:,1:q)*(sigma(1:q,1:q))^(-1/2);
            expW=((sigma(1:q,1:q))^(-1/2)*Vsvd(:,1:q)'*R)'; 

            expSolution={expW'*full(A)*expV, expW'*full(B),full(C)*expV,  expV, expW};
            verification(testCase, actSolution, expSolution, sysr);
            verifyLessThan(testCase, norm((W'*V)-eye(size(V,2))),1.01,...
               'W*V not identity matrix');
        end
        function testTbr5(testCase) %q=10 (with E-matrix)
            load('LF10.mat');
            q=10;
            
            E=blkdiag(speye(size(M)), M);
            A=[zeros(size(M)),speye(size(M)); -K, -D];
            B=[zeros(size(M,1),1); B];
            C=[C, zeros(1,size(M,1))];
            
            [sysr, V, W] = tbr(sss(E\A,E\B,C,0),q);
            actSolution={full(sysr.A),full(sysr.B),full(sysr.C),V,(W'/E)'};
            
            S=lyapchol(full(E\A),full(E\B));
            R=lyapchol(full((E\A)'),full(C)');
            [Usvd,sigma,Vsvd] = svd(S*R');
            
            expV=S'*Usvd(:,1:q)*(sigma(1:q,1:q))^(-1/2);
            expW=((sigma(1:q,1:q))^(-1/2)*Vsvd(:,1:q)'*R/E)'; %warum W/E?
            
            expSolution={expW'*full(A)*expV, expW'*full(B),full(C)*expV, expV, expW};
            
            verification(testCase, actSolution, expSolution, sysr);
            verifyEqual(testCase, full(sysr.E),  expW'*E*expV,'RelTol',0.2,'AbsTol',0.0000000001,...
                    'sysr.E');
        end
%         function testTbr6(testCase) %adi
%             for i=1:length(testCase.sysCell)
%                 sys=testCase.sysCell{i};
%                 if ~sys.isDae && sys.n>100
%                     q=50;
%                     opts.type='adi';
%                     [~,~,~,actHsv]=tbr(sys,q,opts);
%                     opts.type='tbr';
%                     [~,~,~,expHsv]=tbr(sys,q,opts);
%                     
%                     actSolution={actHsv(1:5)};
%                     expSolution={expHsv(1:5)};
% 
%                     verifyEqual(testCase, actSolution, expSolution,'RelTol',0.3,...
%                         'Difference between actual and expected exceeds relative tolerance');
%                 end
%             end
%         end
        function testTbr7(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae && sys.isSiso && sys.n>100
                    Opts.type='tbr';
                    Opts.redErr=1e-10;
                    [sysr,~,~,hsv]=tbr(sys,Opts);
                    [impResSysr,t]=step(ss(sysr));
                    impResSys=step(ss(sys),t);
                    hsvError=(sum(hsv(sysr.n+1:end))+hsv(end)*(sys.n-length(hsv)))/hsv(1)*2;
                    verifyLessThanOrEqual(testCase, norm(impResSys-impResSysr)/length(t), hsvError);
                end 
            end
        end
        function testMatchDcGain(testCase)
            warning('on','tbr:rcond');
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    Opts.type='matchDcGain';
                    sysr=tbr(sys,10,Opts);
                    w = warning('query','last');
                    if isempty(w) || ~strcmp(w.identifier,'tbr:rcond') %A22 not close to singular
                        actSolution= freqresp(sysr,0);
                        expSolution= freqresp(sys,0);
                        verifyLessThanOrEqual(testCase, norm(abs(actSolution)-abs(expSolution)), 1e-6);
                    end
                end
            end
        end
        function testLyapchol(testCase)
            Opts.lse='gauss';
            Opts.hsvTol=1e-15;
            Opts.type='adi';
            Opts.warnOrError='warn';
            q=10;
            
            for k=1:length(testCase.sysCell)
                sys=testCase.sysCell{k};
                if ~sys.isDae && sys.n>100
                    % eqn struct: system data
                    eqn=struct('A_',sys.A,'E_',sys.E,'B',sys.B,'C',sys.C,'type','N','haveE',sys.isDescriptor);

                    % opts struct: mess options
                    messOpts.adi=struct('shifts',struct('l0',20,'kp',50,'km',25,'b0',ones(sys.n,1),...
                        'info',0),'maxiter',300,'restol',0,'rctol',1e-12,...
                        'info',0,'norm','fro');

                    % user functions: default
                    if strcmp(Opts.lse,'gauss')
                        oper = operatormanager('default');
                    elseif strcmp(Opts.lse,'luChol')
                        if sys.isSym
                            oper = operatormanager('chol');
                        else
                            oper = operatormanager('lu');
                        end
                    end

                    if exist('q','var') %size of cholesky factor [sys.n x q] -> qmax=q
                        messOpts.adi.maxiter=q;
                        messOpts.adi.restol=0;
                        messOpts.adi.rctol=1e-30;
                    end

                    % get adi shifts
                    [messOpts.adi.shifts.p, eqn]=mess_para(eqn,messOpts,oper);

                    % low rank adi
                    [R,Rout,eqn]=mess_lradi(eqn,messOpts,oper);

                    if sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0
                        L=R;
                    else
                        eqn.type='T';
                        [L,Lout,eqn]=mess_lradi(eqn,messOpts,oper);
                    end

                    qmax=min([size(R,2),size(L,2)]); %make sure R and L have the same size
                    R=R(:,1:qmax);
                    L=L(:,1:qmax);

                    if exist('q','var') % warn user if rctol is satisfied before q_user
                        qminR=q;
                        qminL=q;
                        nStop=0;

                        % rctol is satisfied if rc<tol for 10 times consecutively
                        for i=1:length(Rout.rc)
                            if Rout.rc(i)<1e-9
                                nStop=nStop+1;
                            else
                                nStop=0;
                            end
                            if nStop==10
                                qminR=i;
                                break
                            end
                        end

                        if ~(sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0) && qminR<q
                            qminL=q;
                            nStop=0;
                            for i=1:length(Lout.rc)
                                if Lout.rc(i)<1e-9
                                    nStop=nStop+1;
                                else
                                    nStop=0;
                                end
                                if nStop==10
                                    qminL=i;
                                    break
                                end
                            end
                        elseif sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0
                            qminL=qminR;
                        end
                        q_min_in=max(qminR,qminL);

                        if q_min_in>0 && q_min_in < q && strcmp(Opts.warnOrError,'warn')
                            warning(['After q=', num2str(q_min_in,'%d'),...
                            ' the contribution of the ADI iterates was very small. Consider reducing the desired order accordingly.']);
                        end
                    end

                    % calculate balancing transformation and Hankel Singular Values
                    [K,S,M]=svd(full(L'*R));
                    expHsv = real(diag(S));
                    expTBalInv = R*M/diag(sqrt(expHsv));
                    expTBal = diag(sqrt(expHsv))\K'*L'/eqn.E_;

                    expV = expTBalInv(:,1:q);
                    expW = expTBal(1:q,:)';

                    expSysr=sss(expW'*sys.A*expV,expW'*sys.B,sys.C*expV,sys.D,expW'*sys.E*expV);

                    lyapOpts.type='adi';
                    lyapOpts.q=10;
                    [R,L]=lyapchol(sys,lyapOpts);


                    % calculate balancing transformation and Hankel Singular Values
                    [K,S,M]=svd(full(R*L'));
                    actHsv = real(diag(S));
                    actTBalInv = R'*K/diag(sqrt(actHsv));
                    actTBal = diag(sqrt(actHsv))\M'*L/sys.E;

                    actV = actTBalInv(:,1:q);
                    actW = actTBal(1:q,:)';

                    actSysr=sss(actW'*sys.A*actV,actW'*sys.B,sys.C*actV,sys.D, actW'*sys.E*actV);

                    expSolution={expHsv, abs(expTBal), abs(expTBalInv), cplxpair(eig(full(expSysr.A)))};
                    actSolution={actHsv, abs(actTBal), abs(actTBalInv), cplxpair(eig(full(actSysr.A)))};

                    verifyEqual(testCase, expSolution, actSolution,'RelTol',0.3,...
                         'Difference between actual and expected exceeds relative tolerance');
                end
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution, sysr)
       verifyEqual(testCase, actSolution, expSolution,'RelTol',0.3,...
            'Difference between actual and expected exceeds relative tolerance');
       verifyLessThanOrEqual(testCase, max(imag(sysr.A)), 0, ...
            'Ar is not purely real'); 
       verifyLessThanOrEqual(testCase, max(imag(sysr.E)), 0, ...
            'Er is not purely real'); 
       verifyEqual(testCase, rank(full(sysr.A)), length(sysr.B),...
            'Rank(Ar) is not full');
       verifyEqual(testCase, rank(full(sysr.E)), length(sysr.B),...
            'Rank(Er) is not full');
       verifyEqual(testCase, nnz(isinf(sysr.A)), 0, ...
            'Ar contains Inf');
       verifyEqual(testCase, nnz(isinf(sysr.E)), 0, ...
            'Er contains Inf');
       verifyEqual(testCase, nnz(isnan(sysr.A)), 0, ...
            'Ar contains Nan');
       verifyEqual(testCase, nnz(isnan(sysr.E)), 0, ...
            'Er contains Nan');
       verifyLessThan(testCase, max(real(eig(sysr))), 0, ...
            'Ar is not stable');
end