function [OFDM,DATA] = ofdm_receive(DATA,OFDM,RF,CHANNEL,SIM,FEC,frame)

if RF.PassBandProcessingToggle
    % downconvert to baseband
    [OFDM,RF] = rf_downconvert(OFDM,RF);
end

cext_rem = OFDM.TxAirChannel;

if SIM.ChannelEstimation
    cext_rem(1:OFDM.NumCyclicPilotSymsPerFrame) = [];
    OFDM.RxFFT = fft(cext_rem);
    [OFDM,CHANNEL] = channel_estimate(CHANNEL,OFDM,SIM);
else
    % no estimation algorithm, no pilot signal
    % ideal, take the channel itself, no pilot signal
    cext_rem(1:OFDM.NumCyclicSymsPerFrame) = [];
    OFDM.RxFFT = fft(cext_rem);
    CHANNEL.EstFreqResponse = fft(CHANNEL.ImpulseResponse,OFDM.NumCarriersPerFrame);
    OFDM.IdMatrixDiagLength = OFDM.NumCarriersPerFrame;
end

if CHANNEL.Equalization
    [OFDM,CHANNEL] = channel_apply_eq(CHANNEL,OFDM,SIM);
end

OFDM.RxDemod = qamdemod(OFDM.RxFFT,OFDM.M,OutputType='bit',UnitAveragePower=1);
if SIM.FECToggle
    OFDM = data_fec_decode(OFDM,FEC);
end

if SIM.Interleave
    OFDM.RxDemod = data_deinterleave(OFDM.RxDemod,FEC);
end

DATA.BitDataHat((frame-1)*OFDM.NumBitsPerFrame+1:...
    (frame-1)*OFDM.NumBitsPerFrame+OFDM.NumBitsPerFrame) = OFDM.RxDemod;
end