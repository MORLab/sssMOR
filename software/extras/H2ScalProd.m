function ip = H2ScalProd(sys1, sys2)
% H2 inner product  - alpha-Version!!

if sys1.isDae || sys2.isDae
    error('this function does not work for DAEs') 
end

ip = trace(sys1.B'*lyap(sys1.A',sys2.A,sys1.C'*sys2.C)*sys2.B);

