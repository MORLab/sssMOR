clear
clc


% Declare input 
sys = loadSss('iss');
s0 = [1     2-5i     3+4i     3-4i        2     2+5i];

Rt = [1     4     1+2i     1-2i        7-8i     7+8i;
      2     5     3+4i     3-4i        2-3i     2+3i;
      3     6     5+6i     5-6i        8-6i     8+6i];
  
Lt = [1     4     1+2i     1-2i        7-8i     7+8i;
      2     5     3+4i     3-4i        2-3i     2+3i;
      3     6     5+6i     5-6i        8-6i     8+6i];
  
sout = [1   2    7-4i     7+4i     2+5i        2-5i];
index = isreal(s0)
%[so_inp,tangential] = shiftmatrix(s0);

[sysr] = rk(sys, s0, Rt);




% function [snew] = adaptive_shift(Ar,s1,s2)
% 
% s_vect(1,1) = s1;       s_vect(1,2) = s2;
% 
% % compute eigenvalues of Ar
% eigAr = eig(Ar);
% 
% % build spectral set of subspace (spec_set is a vector) / eHpoints ist spec
% % set
% if isreal(s1) && isreal(s2)
%     spec_set = sort([s1 -s1 s2 -s2 (-eigAr)']);
% elseif isreal(s1) && ~isreal(s2)
%     spec_set = sort([s1 -s1 s2 conj(s2) (-eigAr)']);
% elseif ~isreal(s1) && isreal(s2)
%     spec_set = sort([s1 conj(s1) s2 -s2 (-eigAr)']);
% else
%     spec_set = sort([s1 conj(s1) s2 conj(s2) (-eigAr)']);
% end
% 
% % build convex hull or real poles
% if any(abs(imag(spec_set))) == 1
%     hull_conv = convhull(real(spec_set),imag(spec_set));
%     spec_set = hull_conv;
% end
% 
% % build residual; product das mit den z siehe formel 2.2
% for ii=1:1:length(spec_set)-1
%     v = linspace(spec_set(ii),spec_set(ii+1),500);
%     for jj=1:1:length(v)
%         res(jj) = prod((v(jj)-s1)./(v(jj)-eigAr));
%     end
%     [~,temp] = max(abs(res));
%     snew(ii) = v(temp);
% end
% for ii=1:1:length(snew)
%     res(ii) = prod((snew(ii)-s1)./(snew(ii)-eigAr));
% end
% [~,temp] = max(abs(res));
% snew = snew(temp);  
% s_vect(1,ii) = snew; 
% 
% 
% end










