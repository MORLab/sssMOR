classdef testSSS < matlab.unittest.TestCase
%Nicht fertig/aktuell    
     
    properties
    end

 
    methods(Test)
         function testbode (testCase)
              %cell2mat + Stützpunkte + omega transponieren + size(1:1:...)
              %von mag und phase
           
              load('build.mat');
              
              [mag_sss, phase_sss, omega] = bode(sss(A,B,C),1:1000); %ss: [mag,phase,wout] = bode(sys)

              %Fehlt in Code:
              mag=zeros(1,1,1000);
              mag(1,:,:)=cell2mat(mag_sss);
              phase(1,:,:)=cell2mat(phase_sss);
       
              actSolution={mag, phase, omega'};
              
              [expmag, expphase, expomega]=bode(ss(full(A),full(B),full(C),0),1:1000);
              expSolution={expmag, expphase, expomega};
              
              verification (testCase, actSolution, expSolution);
              verifyInstanceOf(testCase, mag , 'double', 'Instances not matching');
              verifyInstanceOf(testCase, phase , 'double', 'Instances not matching');
              verifyInstanceOf(testCase, omega , 'double', 'Instances not matching');
              verifySize(testCase, full(mag), size(expmag), 'Size not matching');
              verifySize(testCase, full(phase), size(expphase), 'Size not matching');
              verifySize(testCase, full(omega'), size(expomega), 'Size not matching');
         end
         
%          function testdiag (testCase)         
%               %expsysd.B - calculation missing
%               %last entries of sysd.A differ from expsysd.A -
%               %disp(full(sysd.A)-expsysd.A);

%               load('build.mat');
%               
%               sysd = diag(sss(A, B, C));
%               [expsysd]=canon(ss(full(A), full(B), full(C),0));       
%              
%               actSolution={full(sysd.A)};
%               expSolution={expsysd.A};
%               
%               verification (testCase, actSolution, expSolution);
%               verifyInstanceOf(testCase, sysd.A , 'double', 'Instances not matching');
%               verifySize(testCase, full(sysd.A), size(A), 'Size not matching');
%          end

         function testfreqresp (testCase)         
              %cell2mat + dimension (1:1:x)
              
              load('build.mat');
              
              G_sss = freqresp(sss(A,B,C), [1+2i,3+7i, -4-8i]);
              expG = freqresp(ss(full(A), full(B), full(C),0),[1+2i,3+7i, -4-8i]);
              
              %Fehlt:
              G(1,:,:)=cell2mat(G_sss);
              
              actSolution={full(G)};
              expSolution={expG};
              
              verification (testCase, actSolution, expSolution);
              verifyInstanceOf(testCase, full(G) , 'double', 'Instances not matching');
              verifySize(testCase, full(G), size(expG), 'Size not matching');
         end
         
         function testimpulse (testCase)         
              %cell2mat + Ergebnis transponieren + Stützpunkte
              
              load('build.mat');

              [h, t] = impulse(sss(A,B,C),1:1000); %ss: [y,t] = impulse(sys,Tfinal) oder impulse(sys,t) 
              [exph,expt] = impulse(ss(full(A), full(B), full(C),0),1:1000);


              h=cell2mat(h);
              actSolution={h', t'};
              expSolution={exph, expt};
              
              verification (testCase, actSolution, expSolution);
              verifyInstanceOf(testCase, h' , 'double', 'Instances not matching');
              verifyInstanceOf(testCase, t' , 'double', 'Instances not matching');
              verifySize(testCase, h', size(exph), 'Size not matching'); 
              verifySize(testCase, t', size(expt), 'Size not matching');
         end
         
%          function testisstable (testCase)         
%               %Output von sss? -> Matlab-ss: B={1 wenn stabil, 0 sonst} +
%               %Header fehlt

%               load('build.mat');
% 
%               sys=is_stable(sss(A,B,C));
%               out = isstable(ss(full(A), full(B), full(C),0));
%               disp(sys.Stability);
% 
%               actSolution={(sys.Stability=='stable')};
%               expSolution={out};
%               
%               verification (testCase, actSolution, expSolution);
%          end

         function testnorm1 (testCase)         
              load('build.mat');

              nrm = norm(sss(A,B,C));
              expnrm = norm(ss(full(A), full(B), full(C),0),2);
             

              actSolution={nrm};
              expSolution={expnrm};
              
              verification (testCase, actSolution, expSolution);
              
              verifyInstanceOf(testCase, nrm , 'double', 'Instances not matching');
              verifySize(testCase, nrm, size(expnrm), 'Size not matching');
         end
         
         function testnorm2 (testCase)         
              load('build.mat');

              [nrm, fpeak] = norm(sss(A,B,C), Inf);
              [expnrm, expfpeak] = norm(ss(full(A), full(B), full(C),0),Inf);
             

              actSolution={nrm, fpeak};
              expSolution={expnrm, expfpeak};
              
              verification (testCase, actSolution, expSolution);
              
              verifyInstanceOf(testCase, nrm , 'double', 'Instances not matching');
              verifyInstanceOf(testCase, fpeak , 'double', 'Instances not matching');
              verifySize(testCase, nrm, size(expnrm), 'Size not matching');
              verifySize(testCase, fpeak, size(expfpeak), 'Size not matching');
         end
         
%          function testpzmap (testCase)  
%               % z = cell2mat(z)
%               % last entries of z and expz differ

%               load('build.mat');
% 
%               [p, z] = pzmap(sss(A,B,C));
%               [expp,expz] = pzmap(ss(full(A), full(B), full(C),0));
%               
%               z=cell2mat(z);
%               
% %               disp(abs(z-expz));
%               
%               actSolution={p,z};
%               expSolution={expp, expz};
%               
%               verification (testCase, actSolution, expSolution);
%               
%               verifyInstanceOf(testCase, p , 'double', 'Instances not matching');
%               verifyInstanceOf(testCase, z , 'double', 'Instances not matching');
%               verifySize(testCase, p, size(expp), 'Size not matching');
%               verifySize(testCase, z, size(expz), 'Size not matching');
%          end
         
         function testresidue (testCase)         
              %full(eig) -> evtl diag.m verwenden (sysd -> gleiches
              %Residuum, aber anderes A (Modalform?))
              
              load('build.mat');

              [r,p,d] = residue(sss(A,B,C));
              sysd = diag(sss(A,B,C));
              r=cell2mat(r);
              
              [T,lambda]=eig(full(A));   
              sysT = ss2ss(ss(full(A), full(B), full(C),0),inv(T));
              

              actSolution={r',p.',full(d), r'};
              expSolution={sysT.B.*sysT.C',diag(lambda),0, full(sysd.B.*sysd.C')};
              
              verification (testCase, actSolution, expSolution);
         end

         function testsigma (testCase)       
              %cell2mat + 2nd output (omega) transponieren + andere
              %frequenzvektoren zum plotten (wenn nicht vorgegeben)
              
              load('build.mat');

              [mag, omega] = sigma(sss(A,B,C), 1:100);
              [expmag,expomega] = sigma(ss(full(A), full(B), full(C),0),1:100);      
              mag=cell2mat(mag);

              actSolution={mag, omega'};
              expSolution={expmag, expomega};
              
              verification (testCase, actSolution, expSolution);
              verifyInstanceOf(testCase, mag , 'double', 'Instances not matching');
              verifyInstanceOf(testCase, omega , 'double', 'Instances not matching');
              verifySize(testCase, mag, size(expmag), 'Size not matching');
              verifySize(testCase, omega', size(expomega), 'Size not matching');
         end
         
         function teststep (testCase)       
              %cell2mat(h), h und t transponieren, Stützpunkte anders als
              %ss-step
              %Abs-Tol da sehr kleine Werte mit RelTol nicht gut
              %funktionieren (Fail obwohl Fehler im Bereich 10^-16)
              
              load('build.mat');

              [h, t] = step(sss(A,B,C),1:1000);
              [exph,expt] = step(ss(full(A), full(B), full(C),0),1:1000);    

              h=cell2mat(h);
              actSolution={h',t'};
              expSolution={exph, expt};
              
               verifyEqual(testCase, actSolution, expSolution, 'AbsTol', 0.000000001, ...
               'Difference between actual and expected exceeds relative tolerance');
              verifyInstanceOf(testCase, h , 'double', 'Instances not matching');
              verifyInstanceOf(testCase, t , 'double', 'Instances not matching');
              verifySize(testCase, h', size(exph), 'Size not matching');
              verifySize(testCase, t', size(expt), 'Size not matching');
         end
         
%          function testsim (testCase)       
%          % Header fehlt
% 
%          load('build.mat');
% %          TODO
%          end
    end
end

function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 0.3, ...
               'Difference between actual and expected exceeds relative tolerance');
end