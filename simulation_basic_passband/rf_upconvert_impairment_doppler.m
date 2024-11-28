function y = rf_upconvert_impairment_doppler(x,RF)
% Doppler effect causes Carrier Frequency Offset (CFO)
dopplershift = RF.IMPAIRMENT.DOPPLER.FrequencyShift; % Get the Doppler frequency shift
y = x.*exp(2j*pi*dopplershift*RF.TimeGeneration);
end