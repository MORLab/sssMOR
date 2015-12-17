function [] = compareOrthogonalization()
% Compares the effectiveness of (re)orthogonalization with dkgs or gs by 
% evaluating the orthogonality of V and the required time (ONLY
% SISO-systems)
%
%   LOOP    -   END  (loop: after new direction is calculated, end: reorth)
% 1 DGKS    -    0
% 2 2xSGSN  -    0   (normalize after first GS-orthogonalization)
% 3 2xMGSN  -    0   (normalize after first GS-orthogonalization)
% 4 1xMGS   -  1xMGS
% 5 2xSGSO  -    0   (don't normalize after first GS-orthogonalization = DGKS with 2 iterations)
% 6 2xMGSO  -    0   (don't normalize after first GS-orthogonalization)

%choose type 'full' and a high number
loadBenchmarks; 

temp=load('benchmarksSysCell.mat');
sysCell=temp.benchmarksSysCell;
if isempty(sysCell)
    error('No benchmarks loaded.');
end

% add or remove s0 (require the same length)
s0=zeros(4,10);
s0(1,:)=[1,1,1,1,1,1,1,1,1,1]*10;
s0(2,:)=[1,1,1,1,1,1,1,1,1,1]*0;
s0(3,:)=[0, 0.1,  0.001, 0.0001, 0.000010, 0.000001, 0.0000001, 0.00000001, 0.000000001, 0.0000000001];
s0(4,:)=[1,2,3,4,5,6,7,8,9,10]/100000000;
s0(5,:)=[1,2,3,4,5,6,7,8,9,10]*100000;
s0(6,:)=[1,2,3,4,5,6,7,8,9,10];
s0(7,:)=[4+12.999i, 4-12.999i, 4+13i, 4-13i, 4+13.0001i, 4-13.0001i, 3.999+13i, 3.999-13i, 4.0001+13i, 4.0001-13i];

time1=zeros(length(sysCell),size(s0,1));
time2=zeros(length(sysCell),size(s0,1));
time3=zeros(length(sysCell),size(s0,1));
time4=zeros(length(sysCell),size(s0,1));
time5=zeros(length(sysCell),size(s0,1));
time6=zeros(length(sysCell),size(s0,1));
normV1=zeros(length(sysCell),size(s0,1));
normV2=zeros(length(sysCell),size(s0,1));
normV3=zeros(length(sysCell),size(s0,1));
normV4=zeros(length(sysCell),size(s0,1));
normV5=zeros(length(sysCell),size(s0,1));
normV6=zeros(length(sysCell),size(s0,1));

opts.krylov=0;
opts.dgksTol=5e-14;
for j=1:size(s0,1)
    disp([10,10, 's0(j):']);
    disp(s0(j,:));
    for i=1:length(sysCell)
        sys = sysCell{i};
        if(size(sys.B,2)==1)
            % 1 DGKS - 0
            opts.orth='dgks';
            opts.reorth=0;
            start1=tic();
            [V1] = arnoldi(sys.E,sys.A,sys.B,s0(j,:), @(x,y) (x'*y), opts);
            time1(i,j)=toc(start1);
            normV1(i,j)=norm(V1'*V1-speye(size(V1,2)),'fro');
            
            % 2 2xSGSN - 0
            opts.orth='2sgsn';
            opts.reorth=0;
            start=tic();
            [V2] = arnoldi(sys.E,sys.A,sys.B,s0(j,:), @(x,y) (x'*y), opts);
            time2(i,j)=toc(start);
            normV2(i,j)=norm(V2'*V2-speye(size(V2,2)),'fro');
            
            % 3 2xMGSN - 0
            opts.orth='2mgsn';
            opts.reorth=0;
            start=tic();
            [V3] = arnoldi(sys.E,sys.A,sys.B,s0(j,:), @(x,y) (x'*y), opts);
            time3(i,j)=toc(start);
            normV3(i,j)=norm(V3'*V3-speye(size(V3,2)),'fro');        
            
             % 4 1xMGS - 1xMGS
            opts.orth='mgs';
            opts.reorth='mgs';
            start=tic();
            [V4] = arnoldi(sys.E,sys.A,sys.B,s0(j,:), @(x,y) (x'*y), opts);
            time4(i,j)=toc(start);
            normV4(i,j)=norm(V4'*V4-speye(size(V4,2)),'fro');  
            
            % 5 2xMGSO - 0
            opts.orth='2sgso';
            opts.reorth=0;
            start=tic();
            [V5] = arnoldi(sys.E,sys.A,sys.B,s0(j,:), @(x,y) (x'*y), opts);
            time5(i,j)=toc(start);
            normV5(i,j)=norm(V5'*V5-speye(size(V5,2)),'fro');        
            
             % 4 2xSGSO - 0
            opts.orth='2mgso';
            opts.reorth=0;
            start=tic();
            [V6] = arnoldi(sys.E,sys.A,sys.B,s0(j,:), @(x,y) (x'*y), opts);
            time6(i,j)=toc(start);
            normV6(i,j)=norm(V6'*V6-speye(size(V6,2)),'fro');  
        end
    end
end

disp([10, 'Result:',10]);
disp('Best orthgonality: norm(V^T*V-eye,fro)');
disp([num2str(nnz(normV1<normV2 & normV1<normV3 & normV1<normV4 & normV1<normV5 & normV1<normV6 & normV1>zeros(length(sysCell),size(s0,1))),'%i') ' DGKS/0']);
disp([num2str(nnz(normV2<normV1 & normV2<normV3 & normV2<normV4 & normV2<normV5 & normV2<normV6 & normV2>zeros(length(sysCell),size(s0,1))),'%i') ' 2xSGSN/0']);
disp([num2str(nnz(normV3<normV1 & normV3<normV2 & normV3<normV4 & normV3<normV5 & normV3<normV6 & normV3>zeros(length(sysCell),size(s0,1))),'%i') ' 2xMGSN/0']);
disp([num2str(nnz(normV4<normV1 & normV4<normV2 & normV4<normV3 & normV4<normV5 & normV4<normV6 & normV4>zeros(length(sysCell),size(s0,1))),'%i') ' 1xMGS/1xMGS']);
disp([num2str(nnz(normV5<normV1 & normV5<normV2 & normV5<normV3 & normV5<normV4 & normV5<normV6 & normV5>zeros(length(sysCell),size(s0,1))),'%i') ' 2xSGSO/0']);
disp([num2str(nnz(normV6<normV1 & normV6<normV2 & normV6<normV3 & normV6<normV4 & normV6<normV5 & normV6>zeros(length(sysCell),size(s0,1))),'%i') ' 2xMGSO/0']);


disp([10,'Shortest time: tic-toc']);
disp([num2str(nnz(time1<time2 & time1<time3 & time1<time4 & time1<time5 & time1<time6 & time1>zeros(length(sysCell),size(s0,1))),'%i') ' DGKS/0']);
disp([num2str(nnz(time2<time1 & time2<time3 & time2<time4 & time2<time5 & time2<time6 & time2>zeros(length(sysCell),size(s0,1))),'%i') ' 2xSGSN/0']);
disp([num2str(nnz(time3<time1 & time3<time2 & time3<time4 & time3<time5 & time3<time6 & time3>zeros(length(sysCell),size(s0,1))),'%i') ' 2xMGSN/0']);
disp([num2str(nnz(time4<time1 & time4<time2 & time4<time3 & time4<time5 & time4<time6 & time4>zeros(length(sysCell),size(s0,1))),'%i') ' 1xMGS/1xMGS']);
disp([num2str(nnz(time5<time1 & time5<time2 & time5<time3 & time5<time4 & time5<time6 & time5>zeros(length(sysCell),size(s0,1))),'%i') ' 2xSGSO/0']);
disp([num2str(nnz(time6<time1 & time6<time2 & time6<time3 & time6<time4 & time6<time5 & time6>zeros(length(sysCell),size(s0,1))),'%i') ' 2xMGSO/0']);

for i=1:size(s0,1)
disp([10, 'Best orthgonality of shift ' num2str(i,'%i')]);
disp([num2str(nnz(normV1(:,i)<normV2(:,i) & normV1(:,i)<normV3(:,i) & normV1(:,i)<normV4(:,i) & normV1(:,i)<normV5(:,i) & normV1(:,i)<normV6(:,i) & normV1(:,i)>zeros(length(sysCell),1)),'%i') ' DGKS/0']);
disp([num2str(nnz(normV2(:,i)<normV1(:,i) & normV2(:,i)<normV3(:,i) & normV2(:,i)<normV4(:,i) & normV2(:,i)<normV5(:,i) & normV2(:,i)<normV6(:,i) & normV2(:,i)>zeros(length(sysCell),1)),'%i') ' 2xSGSN/0']);
disp([num2str(nnz(normV3(:,i)<normV1(:,i) & normV3(:,i)<normV2(:,i) & normV3(:,i)<normV4(:,i) & normV3(:,i)<normV5(:,i) & normV3(:,i)<normV6(:,i) & normV3(:,i)>zeros(length(sysCell),1)),'%i') ' 2xMGSN/0']);
disp([num2str(nnz(normV4(:,i)<normV1(:,i) & normV4(:,i)<normV2(:,i) & normV4(:,i)<normV3(:,i) & normV4(:,i)<normV5(:,i) & normV4(:,i)<normV6(:,i) & normV4(:,i)>zeros(length(sysCell),1)),'%i') ' 1xMGS/1xMGS']);
disp([num2str(nnz(normV5(:,i)<normV1(:,i) & normV5(:,i)<normV2(:,i) & normV5(:,i)<normV3(:,i) & normV5(:,i)<normV4(:,i) & normV5(:,i)<normV6(:,i) & normV5(:,i)>zeros(length(sysCell),1)),'%i') ' 2xSGSO/0']);
disp([num2str(nnz(normV6(:,i)<normV1(:,i) & normV6(:,i)<normV2(:,i) & normV6(:,i)<normV3(:,i) & normV6(:,i)<normV4(:,i) & normV6(:,i)<normV5(:,i) & normV6(:,i)>zeros(length(sysCell),1)),'%i') ' 2xMGSO/0']);
end

delete('benchmarksSysCell.mat');
end