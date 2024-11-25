function s = ofdm_transmit_subcarrier_map(OFDM,SIM)
% Subcarrier mapping and allocation
% OFDM - OFDM structure
% SIM - Simulation control parameter
x = OFDM.TxSymbols;

if SIM.ChannelEstimation % case for pilot signal insertion


else
    
end


