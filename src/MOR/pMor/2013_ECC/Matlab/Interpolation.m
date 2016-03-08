function [sys_int] = Interpolation(p_ref, p_int, sys_ref, Methode)

% Inhalt: Interpolation
% Autor: Matthias Geuﬂ

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = length(sys_ref{1}.A);
p = length(p_ref);

for i_Zeile = 1:q

     for i_Spalte = 1:q
         for i_ref = 1:p
             temp_E(i_ref) = sys_ref{i_ref}.E(i_Zeile,i_Spalte);                        
             temp_A(i_ref) = sys_ref{i_ref}.A(i_Zeile,i_Spalte);
         end
         E_int(i_Zeile, i_Spalte) = interp1(p_ref, temp_E, p_int, Methode);          
         A_int(i_Zeile, i_Spalte) = interp1(p_ref, temp_A, p_int, Methode);     
     end

     for i_ref = 1:p
         temp_B(i_ref) = sys_ref{i_ref}.B(i_Zeile);
         temp_C(i_ref) = sys_ref{i_ref}.C(i_Zeile);
     end
     B_int(i_Zeile, 1) = interp1(p_ref, temp_B, p_int, Methode);
     C_int(1, i_Zeile) = interp1(p_ref, temp_C, p_int, Methode);

end
      
sys_int = sss(A_int, B_int, C_int, 0, E_int);
