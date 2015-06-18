function sys=is_stable(sys)
% tests if the System is stable
% letzte Änderung:
% 11.01.11, Sylvia Cremer: handles werden nicht übergeben, zu der
% Entscheidung, ob das System stabil ist oder nicht, muss nur max(real(l))
% untersucht werden
try
    p=pzmap(sys); % versuchen pole zu berechnen
    clear('temp1')
    l=max(real(p));
catch exception
    if strcmp(exception.identifier,'MATLAB:nomem') || strcmp(exception.identifier,'MATLAB:eig:matrixWithNaNInf') %system war zu groß
        try
            l=eigs(sys.A,sys.E,1,'lr'); %versuchen eigenwert mit größtem realteil zu berechnen, konvergiert oft nicht
        catch exception2
            if strcmp(exception2.identifier,'MATLAB:eigs:ARPACKroutineErrorMinus14')
                %eigs ist nicht konvergiert
                opts.tol=1e-4;
                try
                    l=eigs(sys.A,sys.E,1,'lr',opts); % mit geringerer Toleranz eigs probieren
                catch exception3
                    if strcmp(exception3.identifier,'MATLAB:eigs:ARPACKroutineErrorMinus14')
                        l=[];
                    else
                        assignin('base','LastError',exception3)
                        throw exception3
                    end
                end
            else % anderer Fehler
                assignin('base','LastError',exception2)
                throw exception2
            end
        end
    else
        assignin('base','LastError',exception)
        throw(exception)
    end
end
if max(real(l))<0 %Sytem ist stabil
    sys.Stability='stable';
elseif max(real(l))>=0
    sys.Stability='unstable';
else
    sys.Stability='indecisive';
end