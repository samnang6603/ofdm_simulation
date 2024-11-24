function f_hat = rf_downconvert_dpll(x,RF)
% apply DPLL to incoming signal to figure out frequency of incoming signal

% Sampling Parameters
nIterations = RF.ANTENNA.RX.PLL.Iterations;
fs          = RF.SamplingFrequency;   % Nominal sampling frequency
y_ppm       = 50;    % Fractional frequency offset in ppm
f0          = RF.ANTENNA.RX.PLL.NominalClockFrquency;   % Nominal clock frequency
pn_var      = RF.ANTENNA.RX.PLL.PhaseNoiseVariance;  % Phase noise variance
Kp          = RF.ANTENNA.RX.PLL.Kp;  % Proportional Constant
Ki          = RF.ANTENNA.RX.PLL.Ki;  % Integral Constant
pd_choice   = 0;     % Phase Detector Choice (0 -> arctan{.}; 1 -> Im{.})

y = y_ppm*1e-6;
f_offset = y*f0;
F_offset = f_offset/fs;




