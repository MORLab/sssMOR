function check_SD(sys)

% This function checks weather a given syste is in strictly dissipative
% form and if it is stable
%

%%  suppress warning backtrace
% s = warning('off','backtrace');

%%  Perform computations
fprintf('Check SD:\n'), 
% Er>0?
[~,p] = chol(sys.e); 
if ~(p>0)
   fprintf('- Er > 0\n'); 
else
   warning('Er ~> 0'); 
end
% mu2_r < 0?
% if eigs(sys.a+sys.a',1,'LA')< 0
[~,p] = chol(-sys.a-sys.a'); 
if ~(p>0)
    fprintf('- Ar + Ar'' < 0\n'); 
else
    warning('Ar is not SD'); 
end
%   Stable system?
lambda = eig(sys);
lambda = lambda(~isinf(lambda)); %get only finite eigenvalues
if max(real(lambda)) <= 0
    fprintf('- (Er,Ar) is a stable pair \n'); 
else
    warning('The reduced system is unstable!'); 
end


%%  Restore warning settings
% warning(s);