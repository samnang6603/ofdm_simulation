function [OFDM,DATA] = ofdm_receive(DATA,OFDM,RF,CHANNEL,SIM,FEC,frame)

if RF.PassBandProcessingToggle
    % downconvert to baseband
    [OFDM,RF] = rf_downconvert(OFDM,RF,SIM);
end

cext_rem = OFDM.TxAirChannel;

% IQ imbalance correction and mitigation
if RF.IMPAIRMENT.IQImbalance.Mitigation.Toggle
    iqimbal_input = cext_rem;
    cext_rem = rf_downconvert_impairment_IQimb_mitigate(iqimbal_input,RF);
end

% If pilot signal is active
if SIM.PilotSignalToggle
    % remove the cyclic prefix length affected by addition of pilot signal
    cext_rem(1:OFDM.NumCyclicPilotSymsPerFrame) = [];
else 
    % remove unaffected length of cyclic prefix
    cext_rem(1:OFDM.NumCyclicSymsPerFrame) = [];
end

OFDM.RxFFT = fft(cext_rem);

% If channel estimation is active
if SIM.ChannelEstimation
     [OFDM,CHANNEL] = channel_estimate(CHANNEL,OFDM,SIM);
else
    % no estimation algorithm, no pilot signal
    % ideal, take the channel itself, no pilot signal
    CHANNEL.EstFreqResponse = fft(CHANNEL.ImpulseResponse,OFDM.NumCarriersPerFrame);
    OFDM.IdMatrixDiagLength = OFDM.NumCarriersPerFrame;
end

if CHANNEL.Equalization
    [OFDM,CHANNEL] = channel_apply_eq(CHANNEL,OFDM,SIM);
end

% remove pilot signal
if SIM.PilotSignalToggle % remove pilot signal if active
    OFDM.RxFFT(OFDM.PilotSignalLocation) = [];
end

% symbols de-mapping
OFDM.RxDemod = qamdemod(OFDM.RxFFT,OFDM.M,OutputType='bit',UnitAveragePower=1);

% FEC decode
if SIM.FECToggle
    OFDM = data_fec_decode(OFDM,FEC);
end

% Deinterleave
if SIM.Interleave
    OFDM.RxDemod = data_deinterleave(OFDM.RxDemod,FEC);
end

% bit accumulation
DATA.BitDataHat((frame-1)*OFDM.NumBitsPerFrame+1:...
    (frame-1)*OFDM.NumBitsPerFrame+OFDM.NumBitsPerFrame) = OFDM.RxDemod;
end