function OFDM = data_fec_encode(OFDM,FEC)
OFDM.data_frame_codeword = convenc(OFDM.data_frame,FEC.Trellis);
% Update OFDM parameter based on FEC parameter
OFDM.NumBitsPerFrame = length(OFDM.data_frame_codeword);
OFDM.NumCarriersPerFrame = OFDM.NumBitsPerFrame/OFDM.BitsPerSymbol;
OFDM.NumPilotPerFrame = OFDM.NumCarriersPerFrame/OFDM.NumPilotSpacing;
OFDM.NumCarrierPilotPerFrame = OFDM.NumCarriersPerFrame + OFDM.NumPilotPerFrame;
OFDM.NumCyclicSymsPerFrame = floor(OFDM.NumCarriersPerFrame*0.25);
OFDM.NumCyclicPilotSymsPerFrame = floor(OFDM.NumCarrierPilotPerFrame*0.25);
end