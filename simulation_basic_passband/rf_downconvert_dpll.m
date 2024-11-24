function f_hat = rf_downconvert_dpll(x,RF)
% apply DPLL to incoming signal to figure out frequency of incoming signal
% x is input signal
% f_hat is the approximated frequency of the carrier signal
%
% Selected Bibliography:
% [1] Rice, Michael. Digital Communications: A Discrete-Time Approach.
% Appendix C.

% Sampling Parameters
nIterations = RF.ANTENNA.RX.PLL.Iterations;
fs          = RF.SamplingFrequency;   % Nominal sampling frequency
y_ppm       = 50;    % Fractional frequency offset in ppm
f0          = RF.ANTENNA.RX.PLL.NominalClockFrquency;   % Nominal clock frequency
Kp          = RF.ANTENNA.RX.PLL.Kp;  % Proportional Constant
Ki          = RF.ANTENNA.RX.PLL.Ki;  % Integral Constant
pd_choice   = 0;     % Phase Detector Choice (0 -> arctan{.}; 1 -> Im{.})

y = y_ppm*1e-6; % Fractional frequency offset in ppm
f_offset = y*f0; % Absolute freq error (in Hz)
F_offset = f_offset/fs; % Normalized absolute frequency offset

% Nominal Phase Increment (used in the loop phase accumulator)
delta_phi_c = 2*pi*f0/fs;

% PLL pull range
% Normalize frequency offset must be smaller than Kp*pi, namely:
%
%   F_offset < Kp*pi
if (F_offset > Kp*pi)
    warning('Frequency offset is not within the pull range\n');
end

% Preallocate
s_loop             = zeros(nIterations, 1); % DDS Complex Value
phi_loop           = zeros(nIterations, 1); % DDS Phase Accumulator
dds_mult           = zeros(nIterations, 1); % Conjugate Product
phi_error          = zeros(nIterations, 1); % Phase Error
phi_error_filtered = zeros(nIterations, 1); % Filtered Phase Error

% Initialize the loop DDS to a random phase
phi_loop(1)  = sqrt(pi)*randn;

% Initialize integral filter output
integral_out = 0;

for k = 1:nIterations
    % Loop DDS:
    s_loop(k) = exp(1j*phi_loop(k));

    % Phase Detector
    % Multiply the input complex exponential by the conjugate of the loop DDS:
    dds_mult(k) = x(k) * conj(s_loop(k));

    % Phase error:
    if (pd_choice)
        % Phase detector choice: Im{.}
        phi_error(k) = imag(dds_mult(k));
    else
        % Phase detector choice: arctan{.}
        phi_error(k) = angle(dds_mult(k));
    end

% Loop Filter
% The loop filter consists of a Proportional+Integral (PI) controller

% Proportional term
proportional_out = phi_error(k)*Kp;

% Integral term
integral_out     = phi_error(k)*Ki + integral_out;

% PI Filter output:
phi_error_filtered(k) = proportional_out + integral_out;

% Update the phase accumulator of the loop DDS
% The angle of the DDS complex exponential in the next clock cycle results
% from the sum between the nominal phase increment and the filtered phase
% error term:
phi_loop(k+1) = phi_loop(k) + delta_phi_c + phi_error_filtered(k);

end

f_hat = phi_loop(end)*fs/(2*pi)/RF.ANTENNA.RX.PLL.Iterations

f_impaired = RF.CarrierFrequency + RF.IMPAIRMENT.DOPPLER.FrequencyShift % carrier frequency


% Expected filter steady-state value
% The values to which the phase error and the filtered phase error converge
% depend on the loop order. For a second-order loop (with a PI controller),
% assuming an input with fixed frequency offset, the phase error converges
% to 0. In contrast, for a first-order loop (with a Proportional controller
% only), the phase error converges to 2*pi*(f_offset/fs)/Kp. In both cases,
% the filter output converges to 2*pi*(f_offset/fs).
if (Ki == 0)
    phi_error_ss_expected = 2*pi * F_offset / Kp;
else
    phi_error_ss_expected = 0;
end

phi_error_filtered_ss_expected = 2*pi * F_offset;
