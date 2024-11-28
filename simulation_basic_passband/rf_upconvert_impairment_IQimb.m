function y = rf_upconvert_impairment_IQimb(x,RF)
% Apply gain imbalance: Modify I and Q components by gain factor and phase
% shift factor
I_imb = real(x)*(1 + RF.IMPAIRMENT.IQImbalance.IGain);
Q_imb = imag(x)*(1 + RF.IMPAIRMENT.IQImbalance.QGain);

% Apply phase imbalance: Rotate Q by phase imbalance
phase_imb_rad = deg2rad(RF.IMPAIRMENT.IQImbalance.Phase);  % Convert phase imbalance to radians
%Q_imb = Q_imb*cos(phase_imb_rad) - I_imb*sin(phase_imb_rad); % Rotate the Q component
Q_imb = Q_imb.*exp(1i*phase_imb_rad);

% Reconstruct the baseband signal with imbalanced I/Q components
y = I_imb + 1i*Q_imb;
end