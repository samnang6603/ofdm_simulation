function OFDM = data_fec_decode(OFDM,FEC)
OFDM.RxDemod = vitdec(OFDM.RxDemod,FEC.Trellis,FEC.VitDec.TraceBackDepth,...
    FEC.VitDec.OpMode,FEC.VitDec.DecType);
% Re-update OFDM param after decoding
OFDM.NumBitsPerFrame = length(OFDM.data_frame);
end