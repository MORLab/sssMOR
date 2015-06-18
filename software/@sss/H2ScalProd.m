% H2 inner product  - alpha-Version!!

function ip = H2ScalProd(sys1, sys2)
if sys1.is_dae, error, end
if sys2.is_dae, error, end

ip = trace(sys1.B'*lyap(sys1.A',sys2.A,sys1.C'*sys2.C)*sys2.B);

