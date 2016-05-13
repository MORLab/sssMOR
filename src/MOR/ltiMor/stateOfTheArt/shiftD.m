function sysD = shiftD(sys,D,L,R)

switch class(sys)
    case 'sss'
        sysD = sss(sys.A-L.'*D*R, sys.B-L.'*D, sys.C-D*R, D, sys.E);
    case 'ss'
        sysD = dss(sys.A-L.'*D*R, sys.B-L.'*D, sys.C-D*R, D, sys.E);
    otherwise
        error('sssMOR:shiftD:notSupportedClass','The sys object provided is not supported')
end
