function [OFDM,RF] = rf_upconvert(OFDM,RF)
OFDM.RFTransmitSymsLength = length(OFDM.TxAir); % Length of OFDM symbols
OFDM.SymsDurationPerFrame = RF.SamplingPeriod*OFDM.RFTransmitSymsLength; % OFDM symbols period or duration
RF.CarrierSamplesPerSyms = RF.SamplingFrequency/RF.CarrierFrequency;
RF.PassBandSamples = OFDM.RFTransmitSymsLength*RF.CarrierSamplesPerSyms;
RF.TimeGeneration = (0:RF.PassBandSamples-1)*RF.SamplingPeriod; % get the time vector to generate passband carrier
RF.TimeGeneration = RF.TimeGeneration(:);  

% do ZOH to spread the OFDM samples per frame to match the passband samples
% perframe with passband sampling frequency
ofdm_syms_zoh = repmat(OFDM.TxAir,1,RF.CarrierSamplesPerSyms).'; 
ofdm_syms_zoh = ofdm_syms_zoh(:);

% now generate passband carrier signal
RF.PassBandSignalUpconvertGeneration = exp(2*1i*pi*RF.CarrierFrequency*RF.TimeGeneration);

% mismatch length checking
if length(ofdm_syms_zoh) ~= length(RF.PassBandSignalUpconvertGeneration)
    error('Mismatch between the length of OFDM symbols and passband signal.');
end

% upconvert the baseband to passband by multiplying by carrier
OFDM.RFTxAir = RF.ANTENNA.TX.Gain*ofdm_syms_zoh.*RF.PassBandSignalUpconvertGeneration;

%OFDM.RFTxAir = rf_antenna_agc(OFDM.RFTxAir,RF.ANTENNA.TX.AGC.TargetPower,...
                        %RF.ANTENNA.TX.Gain,0.1,RF.ANTENNA.TX.AGC.MinMaxGain);

end