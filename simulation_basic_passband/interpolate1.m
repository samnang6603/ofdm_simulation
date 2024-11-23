function H_interpolated = interpolate1(H_est,pilot_loc,Nfft)
% First, linearly interpolate pilot signal and H_est by 1 sample to avoid
% NaN output by interp1
if pilot_loc(end) < Nfft
    slope = (H_est(end)-H_est(end-1))/(pilot_loc(end)-pilot_loc(end-1));
    H_est = [H_est; H_est(end)+slope*(Nfft-pilot_loc(end))];
    pilot_loc = [pilot_loc, Nfft];
end
H_interpolated = interp1(pilot_loc,H_est,(1:Nfft)','spline');
end