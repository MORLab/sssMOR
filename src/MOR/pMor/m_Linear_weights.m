function [L] = m_Linear_weights(param, pj)

if pj<min(param) || pj>max(param)
    error('parameter out of bounds')
elseif pj==min(param)
    L = eye(length(param),1);
    return
end

[p,x] = sort([pj,param]);
x = find(x==1,1,'first');

L = zeros(length(param), 1);

L(x-1) = p(x+1)-pj;
L(x) = pj-p(x-1);

L = L / (p(x+1)-p(x-1));
