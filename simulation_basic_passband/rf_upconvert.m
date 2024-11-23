function [OFDM,RF] = rf_upconvert(OFDM,RF)
OFDM.RFTransmitSymsLength = length(OFDM.TxAir); % Length of OFDM symbols
OFDM.SymsDurationPerFrame = RF.SamplingPeriod*OFDM.RFTransmitSymsLength; % OFDM symbols period or duration
RF.CarrierSamplesPerSyms = RF.SamplingFrequency/RF.CarrierFrequency;
RF.PassBandSamples = OFDM.RFTransmitSymsLength*RF.CarrierSamplesPerSyms;
