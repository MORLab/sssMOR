function idx = ismemberf2(sNew,sOld,tol)
%
%   check which entries of sNew are in a circle of radius sOld*tol from any
%   entry of sOld

%   doing it with loops, there might be a faster version

idx = zeros(1,length(sNew));
if true
    %implementation 1: for loop
    for iS = 1:length(sNew)
        dist = abs(sOld-sNew(iS));
        idx(iS) = any(dist<abs(sOld)*tol);
    end
else
%     %implementation 2: vectorized
%     [Snew,Sold] = meshgrid(sNew,sOld);
%     dist = abs(Sold-Snew);
%     idx = any(dist<abs(sOld)*tol);
end