function [OFDM,RF] = rf_downconvert(OFDM,RF)

if RF.IMPAIRMENT.DOPPLER.Toggle
    RF = rf_downconvert_carrier_recover(OFDM,RF);
end


% multiply by carrier again
baseband_recovery = RF.ANTENNA.RX.Gain*conj(RF.PassBandSignalDownconvertGeneration)...
                    .*OFDM.RFTxAirChannel;

% create low pass filter
% fpass = 100e3; % Passband at 100 kHz to allow baseband recovery
% fstop = 200e3; % Stopband at 200 kHz, ensuring higher-frequency components are attenuated well
% h = rf_downconvert_lpf(RF.SamplingFrequency, fpass, fstop); % Update filter design
load FIR_blackman_harris_60_taps.mat h

% apply low pass filter (doubles as both baseband lpf and anti-aliasing
% filter)
baseband_recovery_lpf = filter(h,1,baseband_recovery);

% downsample to original OFDM symbol length
OFDM.TxAirChannel = baseband_recovery_lpf(1:RF.CarrierSamplesPerSyms:end);