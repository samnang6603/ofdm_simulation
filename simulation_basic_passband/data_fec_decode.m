function x = data_fec_decode(cw,FEC)
switch FEC.Type
    case 'Convolutional'
        x = vitdec(cw,FEC.Trellis,FEC.VitDec.TraceBackDepth,...
                          FEC.VitDec.OpMode,FEC.VitDec.DecType);
    case 'LDPC'
        % Decoder Reverse the data segmentation
        llr = cw; % variable naming to denote input to LDPC is LLR
        nseg = FEC.NumSegment;
        ninfo = FEC.Ecfg.NumInformationBits;
        cwinfo = FEC.Ecfg.BlockLength;
        x = zeros(FEC.SegInputLength,1,'int8');
        for k = 1:nseg
            % tmp index array
            tmp_idx_x = ninfo*(k-1)+1:ninfo*(k-1) + ninfo;
            tmp_idx_cw = cwinfo*(k-1)+1:cwinfo*(k-1) + cwinfo;
            % decode
            x(tmp_idx_x) = ldpcDecode(llr(tmp_idx_cw),FEC.Dcfg,FEC.Maxnumiter);
        end
        % Remove paddded zeros
        x(end-FEC.NumZeroPadded+1:end) = [];
end
end