function y = rf_downconvert_impairment_IQimb_mitigate(x,tx_pilot,pilot_index)
% IQ Imbalance mitigation by pilot signal estimation
% x: baseband signal after downconversion
% tx_pilot: ideal pilot signal
% pilot_index: index of pilot symbols of OFDM signal

% Extract pilot symbols
x_length = length(x);
x_index_array = 1:x_length;
rx_pilot = x(pilot_index);

% Adjustment correction fraction
% phase_correction_fraction = 0.1;  % 50% correction
% amplitude_correction_fraction = 0.1;  % 50% correction

% Estimate imbalance for Amplitude
amplitude_imbalance_I = abs(real(rx_pilot)) ./ abs(real(tx_pilot));  % I channel imbalance
amplitude_imbalance_Q = abs(imag(rx_pilot)) ./ abs(imag(tx_pilot));  % Q channel imbalance

% Imbalance correction
amplitude_correction_I = 1./amplitude_imbalance_I;
amplitude_correction_Q = 1./amplitude_imbalance_Q;

% Interpolate amplitude correction for the entire signal (I and Q separately)
amplitude_correction_I = interp1(pilot_index, amplitude_correction_I, x_index_array, 'spline', 'extrap').';
amplitude_correction_Q = interp1(pilot_index, amplitude_correction_Q, x_index_array, 'spline', 'extrap').';

% % Apply separate amplitude correction 
% amplitude_imbalance = abs(rx_pilot)./abs(tx_pilot);
% amplitude_correction = 1./amplitude_imbalance;
% amplitude_correction_full = interp1(pilot_index,amplitude_correction,x_index_array,'spline','extrap').';

% Estimate imbalance for Phase
phase_imbalance = angle(rx_pilot) - angle(tx_pilot);
phase_correction = phase_imbalance;
phase_correction_full = interp1(pilot_index,phase_correction,x_index_array,'spline','extrap').';

% Apply compensation to the received signal
% Compensate for amplitude and phase imbalance
% y = x.*exp(-1j*phase_correction_full*phase_correction_fraction)...
%     .*amplitude_correction_full*amplitude_correction_fraction;

y = real(x).*amplitude_correction_I + ...
    1j*imag(x).*amplitude_correction_Q;  % Apply amplitude correction
y = y.*exp(-1j * phase_correction_full);
end
