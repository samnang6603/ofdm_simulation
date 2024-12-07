function FEC = fec_init(SIM)
if ~SIM.FEC.Toggle
    FEC.DemapOutputType = 'bit';
    return
end
switch lower(SIM.FEC.Type)
    case 'ldpc'
        FEC.Type = 'LDPC';
        FEC.P = [
    16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
    25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
    25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
     9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
    24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
     2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0
    ];
        FEC.Z = 27;
        FEC.PCMatrix = ldpcQuasiCyclicMatrix(FEC.Z,FEC.P);
        FEC.Ecfg = ldpcEncoderConfig(FEC.PCMatrix);
        FEC.Dcfg = ldpcDecoderConfig(FEC.PCMatrix,'layered-bp');
        FEC.DemapOutputType = 'approxllr';
        FEC.Maxnumiter = 10;
    case {'conv','convolutional'}
        FEC.Type = 'Convolutional';
        FEC.Constraint = 7;
        FEC.Trellis = poly2trellis(FEC.Constraint,[171 133]);
        FEC.VitDec.TraceBackDepth = 5*(FEC.Constraint-1); % approx for constraint
        FEC.VitDec.OpMode = 'trunc';
        FEC.VitDec.DecType = 'hard';
        FEC.DemapOutputType = 'bit';
    otherwise
        error('Unsupported FEC')
end
end