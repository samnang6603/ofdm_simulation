function [OFDM,RF] = rf_downconvert(OFDM,RF,SIM)
% downconvert to baseband

% if CFO is present and PLL toggle is ON
if RF.ANTENNA.RX.PLL.Toggle
    RF = rf_downconvert_carrier_recover(OFDM,RF);
else
    RF.PassBandSignalDownconvertGeneration = RF.PassBandSignalUpconvertGeneration;
end


% downconvert passband to baseband by multiplying by carrier again
baseband_recovery = RF.ANTENNA.RX.Gain*conj(RF.PassBandSignalDownconvertGeneration)...
                    .*OFDM.RFTxAirChannel;

%baseband_recovery = rf_antenna_agc(baseband_recovery,RF.ANTENNA.RX.AGC.TargetPower,...
 %                   RF.ANTENNA.RX.Gain,0.1,RF.ANTENNA.RX.AGC.MinMaxGain);

% create low pass filter
% fpass = 100e3; % Passband at 100 kHz to allow baseband recovery
% fstop = 200e3; % Stopband at 200 kHz, ensuring higher-frequency components are attenuated well
% h = rf_downconvert_lpf(RF.SamplingFrequency, fpass, fstop); % Update filter design
load FIR_blackman_harris_60_taps.mat h

% apply low pass filter (doubles as both baseband lpf and anti-aliasing
% filter) only if fading and AWGN are present
if SIM.Fading
    baseband_recovery = 4*filter(h,1,baseband_recovery);
end

% downsample to original OFDM symbol length
OFDM.TxAirChannel = baseband_recovery(1:RF.CarrierSamplesPerSyms:end);
