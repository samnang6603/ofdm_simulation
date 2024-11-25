function [y,agcgain] = rf_antenna_agc(x,target_power,current_gain,alpha,gain_limits)
% Calculate RMS power of the complex signal
power = sqrt(mean(abs(x).^2)); % RMS power

% Compute the required gain
desired_gain = sqrt(target_power/power);

% Smoothly update the gain
new_gain = current_gain + alpha*(desired_gain - current_gain);

% Apply gain limits
new_gain = max(min(new_gain, gain_limits(2)), gain_limits(1));

% Adjust the complex signal
y = new_gain*x;

% Output the updated gain
agcgain = new_gain;
end
