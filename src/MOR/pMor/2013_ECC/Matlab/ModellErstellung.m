
function [sys] = ModellErstellung(p)
% Inhalt: Erstellt Modelle Balken, Schwinger, ...
% Autor: Matthias Geuﬂ

% Erl‰uterungen:
% p: Wert des Parameters
% N: Anzahl der Elemente

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modellauswahl

Modell = 'Balken';
%Modell = 'Schwinger';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modellberechnung

if strcmp(Modell,'Balken')
    
    N = 200;                 % Anzahl finiter Elemente
    sys = fem_beam(p, N);
    

    %sys = sss(sys.E\sys.A, sys.E\sys.B, sys.C);


elseif strcmp(Modell,'Schwinger')
    
	                % Anzahl Feder-Massen-Systeme
    
end

