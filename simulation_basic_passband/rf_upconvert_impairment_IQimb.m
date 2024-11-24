function y = rf_upconvert_impairment_IQimb(x,RF)
% Apply gain imbalance: Modify I and Q components by gain factor
I_imb = real(x) * (1 + RF.IMPAIRMENT.IGainImbalance);
Q_imb = imag(x) * (1 - RF.IMPAIRMENT.QGainImbalance);

% Apply phase imbalance: Rotate Q by phase imbalance
phase_imb_rad = deg2rad(RF.IMPAIRMENT.PhaseImbalance);  % Convert phase imbalance to radians
Q_imb = Q_imb * cos(phase_imb_rad) - I_imb * sin(phase_imb_rad); % Rotate the Q component

% Reconstruct the baseband signal with imbalanced I/Q components
y = I_imb + 1i*Q_imb;
end