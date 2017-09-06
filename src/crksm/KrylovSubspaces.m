% testscript um die dimension von Krylov subspaces zu zeigen
%
%
% the aim of this symbolic calculation is to show, that the below rational Krylov
% subspaces build up the same spanning set:
%
%       K(s,r) = span{(A_s1)^-1*b, (A_s1)^-2*b, (A_s2)^-1*(A_s1)^-2*b, (A_s2)^-2*(A_s1)^-2*b}
%       and
%       K(s,r) = span{(A_s1)^-1*b, (A_s1)^-2*b, (A_s2)^-1*b, (A_s2)^-2*b}
%       with the shift vector s = [s1 s1 s2 s2] and the multiplicity r
%
% and more detailed, that span{(A_s2)^-1*(A_s1)^-2*b} = span{(A_s2)^-1*b} yields

clear
clc

% define symbolic variables
syms a11 a12 a13 a21 a22 a23 a31 a32 a33 v11 v21 v31 v12 v22 v32 t1 t2 real

%% SISO CASE
% symbolic matirx (A_s2)^-1
As2_inv = [a11 a12 a13; a21 a22 a23; a31 a32 a33];

% define the previous direction vector v as v = (A_s1)^-2*b
v_siso = [v11; v21; v31];

% calculate (A_s2)^-1*(A_s1)^-2*b = (A_s2)^-1*v for each column of As2_inv 
vnew1 = As2_inv(:,1).*v_siso(1,1);
vnew2 = As2_inv(:,2).*v_siso(2,1);
vnew3 = As2_inv(:,3).*v_siso(3,1);

% calculate the angle between the single vectors (scalar product)
alpha1 = dot(vnew1,vnew2)./(norm(vnew1).*norm(vnew2));
simplify(alpha1);

alpha2 = dot(vnew1,vnew3)./(norm(vnew1).*norm(vnew3));
simplify(alpha2);

alpha3 = dot(vnew2,vnew3)./(norm(vnew2).*norm(vnew3));
simplify(alpha3);

%% MIMO CASE
v_a1 = [a11; a21; a31];
v_a2 = [a12; a22; a32];
v_a3 = [a13; a23; a33];

% define the previous direction vector v as v = (A_s1)^-2*B
V_mimo = [v11 v12; v21 v22; v31 v32];

% calculate the matrix multiplication


%% MIMO TANGENTIAL

% tangential direction and tangential rhs
t = [t1; t2];       vt = V_mimo*t;

vt1 = v_a1*vt(1,:); 
vt2 = v_a2*vt(2,:); 
vt3 = v_a3*vt(3,:); 

% calculation scalar product
beta1 = dot(vt1,vt2)./(norm(vt1).*norm(vt2));
simplify(beta1)






