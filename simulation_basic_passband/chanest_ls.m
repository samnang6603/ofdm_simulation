function H_LS = chanest_ls(Y,Xp,pilot_loc,Nfft)
% reference: Cho MIMO-OFDM Wireless Communications with MATLAB
Np=length(pilot_loc); %Nfft/Nps; 
k=1:Np;
LS_est(k) = Y(pilot_loc(k))./Xp(k); % LS channel estimation
LS_est = LS_est(:);
% Linear/Spline interpolation
H_LS = interpolate1(LS_est,pilot_loc,Nfft);
end