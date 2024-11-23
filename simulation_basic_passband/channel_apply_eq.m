function [OFDM,CHANNEL] = channel_apply_eq(CHANNEL,OFDM,SIM)
H = diag(CHANNEL.EstFreqResponse);
CHANNEL.EqFreqResponse = (H'*H + (1/SIM.SNR*eye(OFDM.IdMatrixDiagLength)))\H';
OFDM.RxFFT = CHANNEL.EqFreqResponse*OFDM.RxFFT;
if SIM.ChannelEstimation % remove pilot signal if channel est being used
    OFDM.RxFFT(OFDM.PilotSignalLocation) = [];
end
end