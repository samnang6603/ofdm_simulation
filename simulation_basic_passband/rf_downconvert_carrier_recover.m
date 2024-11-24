function RF = rf_downconvert_carrier_recover(OFDM,RF)
% implements Phase Locked Loop to recover frequency of incoming signal
pll_signal_in = OFDM.RFTxAirChannel;
f_hat = rf_downconvert_dpll(pll_signal_in,RF);
end
