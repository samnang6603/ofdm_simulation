%% Simulation of OFDM Fading Channel with Channel Estimation and EQ
%{ 
    Comparison of OFDM with Fading Channel and Channel estimation 
    performance using MMSE or LS with EQ at the Rx.

    Author: Sam An
    Date: 08/19/2024
%}


close all; clear; clc

seed = randi([1 1e7],1,1);

% Modulation 16-QAM
NumQam=16;
K=log2(NumQam);
EbN0 = 25;
SNR = EbN0 + 10*log10(sqrt(10));    % Convert EbN0 to SNR
chanlen = 6;

%% Pilot data
NumPilotSpacing = 16;

%% Data
x = imread('cameraman2.tif');
Lx = numel(x);
xhat = zeros(Lx,1,'uint8');
xhat_eq = xhat;
xhat_LSest_eq = xhat;
xhat_MMSEest_eq = xhat;
xbin = num2cell(dec2bin(x(:)));
total_bits = 8*numel(x);
data = zeros(total_bits,1,'uint8');
for m = 0:Lx-1
    str = sprintf(' %s',xbin{m+1,:});
    bit_num = sscanf(str,'%d');
    data(8*m+1:8*m+8) = bit_num;
end
clear xbin
data_hat = zeros(size(data),'uint8');
data_hat_eq = data_hat;
data_hat_LS_eq = data_hat;
data_hat_MMSE_eq = data_hat;

%% Generating and coding data
NumBits = total_bits; % # bits. Can modify this, but has to be power of two
NumBitsPerFrame=256*4;
NumFrames = NumBits/NumBitsPerFrame;
NumCarriersPerFrame=NumBitsPerFrame/K;
NumPilotPerFrame = NumCarriersPerFrame/NumPilotSpacing;
NumCarrierPilotPerFrame = NumCarriersPerFrame + NumPilotPerFrame;
NumCyclicSymsPerFrame = floor(NumCarriersPerFrame*0.25);
NumCyclicPilotSymsPerFrame = floor(NumCarrierPilotPerFrame*0.25);

%% 
for frame = 1:NumFrames
    data_frame = data((frame-1)*NumBitsPerFrame+1:(frame-1)*NumBitsPerFrame+NumBitsPerFrame);
    y = qammod(data_frame,NumQam,InputType='bit',UnitAveragePower=1);

    % Add pilot
    yp = zeros(NumCarrierPilotPerFrame,1);
    Xp = 2*randi([0 1],NumPilotPerFrame,1)-1;
    spc = NumPilotSpacing;
    pilot_loc = 1:spc+1:NumCarrierPilotPerFrame;
    yp(pilot_loc) = Xp;
    for k = 0:NumPilotPerFrame-1
        yp(2+k*(spc+1):1+k*(spc+1)+spc) = y(k*spc+1:k*spc+spc);
    end

    % transmitter no pilot
    ifft_sig = ifft(y,NumCarriersPerFrame);
    cyclic_idx = NumCarriersPerFrame-NumCyclicSymsPerFrame+1:NumCarriersPerFrame;
    cext_data = [ifft_sig(cyclic_idx); ifft_sig];

    % transmitter with pilot
    ifft_sig_pilot = ifft(yp,NumCarrierPilotPerFrame);
    cyclic_idx = NumCarrierPilotPerFrame-NumCyclicPilotSymsPerFrame+1:NumCarrierPilotPerFrame;
    cext_data_pilot = [ifft_sig_pilot(cyclic_idx); ifft_sig_pilot];
    
    % channel without pilot
    h = jake(120,chanlen);
    %h = 1/sqrt(2)*(randn(6,1) + 1j*randn(6,1)); 

    Hvec = fft(h,NumCarriersPerFrame);
    H = diag(Hvec);
    ofdm_sig = filter(h,1,cext_data);
    %rng(seed)  % fix random seed forced seed for comparison
    ofdm_sig = awgn(ofdm_sig,SNR,'measured');

    % channel with pilot
    Hvec_pilot = fft(h,NumCarrierPilotPerFrame);
    H_pilot = diag(Hvec_pilot);
    ofdm_sig_pilot = filter(h,1,cext_data_pilot);
    %rng(seed)
    ofdm_sig_pilot = awgn(ofdm_sig_pilot,SNR,'measured');

    % receiver
    %rng('shuffle')
    cext_rem = ofdm_sig;
    cext_rem(1:NumCyclicSymsPerFrame) = [];
    fft_sig = fft(cext_rem);
    Heq = (H'*H+(1/SNR*eye(NumCarriersPerFrame)))\H'; % Linear MMSE EQ using direct H info
    fft_sig_eq = Heq*fft_sig;

    %receiver with channel estimation
    cext_pilot_rem = ofdm_sig_pilot;
    cext_pilot_rem(1:NumCyclicPilotSymsPerFrame) = [];
    fft_sig_pilot = fft(cext_pilot_rem);

    % channel estimation LS
    H_est_LS = LS_CE1(fft_sig_pilot,Xp,pilot_loc,NumCarrierPilotPerFrame,'spline');
    h_est_LS = ifft(H_est_LS);
    h_est_LS_DFT_chanlen = h_est_LS(1:chanlen);
    H_est_LS_DFT_chanlen = fft(h_est_LS_DFT_chanlen,NumCarrierPilotPerFrame);
    H_est_LS = diag(H_est_LS_DFT_chanlen);
    Heq_LS = (H_est_LS'*H_est_LS+(1/SNR*eye(NumCarrierPilotPerFrame)))\H_est_LS'; % Linear MMSE EQ
    fft_sig_LS_eq = Heq_LS*fft_sig_pilot;  % Apply EQ

    H_est_MMSE = MMSE_CE1(fft_sig_pilot,Xp,pilot_loc,NumCarrierPilotPerFrame,NumPilotSpacing,h.',SNR);
    h_est_MMSE = ifft(H_est_MMSE);
    h_est_MMSE_DFT_chanlen = h_est_MMSE(1:chanlen);
    H_est_MMSE_DFT_chanlen = fft(h_est_MMSE_DFT_chanlen,NumCarrierPilotPerFrame);
    H_est_MMSE = diag(H_est_MMSE_DFT_chanlen);
    Heq_MMSE = (H_est_MMSE'*H_est_MMSE+(1/SNR*eye(NumCarrierPilotPerFrame)))\H_est_MMSE';
    fft_sig_MMSE_eq = Heq_MMSE*fft_sig_pilot;

    % Remove symbol at pilot signal location
    fft_sig_LS_eq(pilot_loc) = [];
    fft_sig_MMSE_eq(pilot_loc) = [];

    data_hat_frame = qamdemod(fft_sig,NumQam,OutputType='bit',UnitAveragePower=1);
    data_hat_eq_frame = qamdemod(fft_sig_eq,NumQam,OutputType='bit',UnitAveragePower=1);
    data_hat_LS_eq_frame = qamdemod(fft_sig_LS_eq,NumQam,OutputType='bit',UnitAveragePower=1);
    data_hat_MMSE_eq_frame = qamdemod(fft_sig_MMSE_eq,NumQam,OutputType='bit',UnitAveragePower=1);

    data_hat((frame-1)*NumBitsPerFrame+1:(frame-1)*NumBitsPerFrame+NumBitsPerFrame) = data_hat_frame;
    data_hat_eq((frame-1)*NumBitsPerFrame+1:(frame-1)*NumBitsPerFrame+NumBitsPerFrame) = data_hat_eq_frame;
    data_hat_LS_eq((frame-1)*NumBitsPerFrame+1:(frame-1)*NumBitsPerFrame+NumBitsPerFrame) = data_hat_LS_eq_frame;
    data_hat_MMSE_eq((frame-1)*NumBitsPerFrame+1:(frame-1)*NumBitsPerFrame+NumBitsPerFrame) = data_hat_MMSE_eq_frame;
end
BER = sum(data~=data_hat)/total_bits;
BER_eq = sum(data~=data_hat_eq)/total_bits;
BER_LSest_eq = sum(data~=data_hat_LS_eq)/total_bits;
BER_MMSEest_eq = sum(data~=data_hat_MMSE_eq)/total_bits;

%% Decode data back into image
data_hat_reshaped = reshape(data_hat,8,NumBits/8).';
data_hat_decode = mat2str(data_hat_reshaped);
data_hat_decode(1) = [];
data_hat_decode(end) = [];

data_hat_eq_reshaped = reshape(data_hat_eq,8,NumBits/8).';
data_hat_eq_decode = mat2str(data_hat_eq_reshaped);
data_hat_eq_decode(1) = [];
data_hat_eq_decode(end) = [];

data_hat_LSest_eq_reshaped = reshape(data_hat_LS_eq,8,NumBits/8).';
data_hat_LSest_eq_decode = mat2str(data_hat_LSest_eq_reshaped);
data_hat_LSest_eq_decode(1) = [];
data_hat_LSest_eq_decode(end) = [];

data_hat_MMSEest_eq_reshaped = reshape(data_hat_MMSE_eq,8,NumBits/8).';
data_hat_MMSEest_eq_decode = mat2str(data_hat_MMSEest_eq_reshaped);
data_hat_MMSEest_eq_decode(1) = [];
data_hat_MMSEest_eq_decode(end) = [];

for m = 0:Lx-1
    start_idx = NumQam*m+1;
    end_idx = NumQam*m+NumQam-1;
    xhat(m+1) = uint8(bin2dec(data_hat_decode(start_idx:end_idx)));
    xhat_eq(m+1) = uint8(bin2dec(data_hat_eq_decode(start_idx:end_idx)));
    xhat_LSest_eq(m+1) = uint8(bin2dec(data_hat_LSest_eq_decode(start_idx:end_idx)));
    xhat_MMSEest_eq(m+1) = uint8(bin2dec(data_hat_MMSEest_eq_decode(start_idx:end_idx)));
end
xhat = reshape(xhat,size(x));
xhat_eq = reshape(xhat_eq,size(x));
xhat_LSest_eq = reshape(xhat_LSest_eq,size(x));
xhat_MMSEest_eq = reshape(xhat_MMSEest_eq,size(x));
figure
subplot(231), imshow(x), title('Original Image')
set(gca,FontSize=14)
subplot(232), imshow(xhat), title('Received Image Rayleigh Fading No EQ')
set(gca,FontSize=14)
subplot(233), imshow(xhat_eq), title('Received Image Rayleigh Fading Perfect LMMSE EQ')
set(gca,FontSize=14)
subplot(234), imshow(xhat_LSest_eq), title('Received Image Rayleigh Fading LS-ChanEst LMMSE EQ')
set(gca,FontSize=14)
subplot(235), imshow(xhat_MMSEest_eq), title('Received Image Rayleigh Fading MMSE-ChanEst LMMSE EQ')
set(gca,FontSize=14)

figure
stem(abs(Hvec),'filled'), hold on,
title('Channel Estimation Comparison')
stem(abs(H_est_LS_DFT_chanlen),'filled')
stem(abs(H_est_MMSE_DFT_chanlen),'filled')
legend('True Channel','LS-DFT','MMSE-DFT')

function H_interpolated = interpolate1(H_est,pilot_loc,Nfft,method)
% First, linearly interpolate pilot signal and H_est by 1 sample to avoid
% NaN output by interp1
if pilot_loc(end) < Nfft
    slope = (H_est(end)-H_est(end-1))/(pilot_loc(end)-pilot_loc(end-1));
    H_est = [H_est; H_est(end)+slope*(Nfft-pilot_loc(end))];
    pilot_loc = [pilot_loc, Nfft];
end
H_interpolated = interp1(pilot_loc,H_est,(1:Nfft)',method);
end

function H_LS = LS_CE1(Y,Xp,pilot_loc,Nfft,int_opt)
Np=length(pilot_loc); %Nfft/Nps; 
k=1:Np;
LS_est(k) = Y(pilot_loc(k))./Xp(k); % LS channel estimation
LS_est = LS_est(:);

if lower(int_opt) == 'l'
    method = 'linear'; 
else 
    method = 'spline'; 
end

% Linear/Spline interpolation
H_LS = interpolate1(LS_est,pilot_loc,Nfft,method);
end

function H_MMSE = MMSE_CE1(Y,Xp,pilot_loc,Nfft,Nps,h,SNR)
snr = 10^(SNR*0.1); 
Np=length(pilot_loc); %floor(Nfft/Nps); 
k=1:Np;
H_tilde = Y(pilot_loc(k)).'./Xp(k).'; % LS estimate Eq. (6.12) or (6.8)
k = 0:length(h)-1; % k_ts = k*ts
hh = h*h';
tmp = h.*conj(h).*k; % tmp = h.*conj(h).*k_ts
r = sum(tmp)/hh; 
r2 = tmp*k.'/hh; % r2 = tmp*k_ts.' /hh;
tau_rms = sqrt(r2-r^2); % rms delay
df = 1/Nfft; % 1/(ts*Nfft);
j2pi_tau_df = 1i*2*pi*tau_rms*df;
K1 = repmat((0:Nfft-1).', 1, Np); 
K2 = repmat((0:Np-1),Nfft,1);
rf = 1./(1+j2pi_tau_df*(K1-K2*Nps)); % Eq. (6.17a)
K3 = repmat((0:Np-1).',1,Np); 
K4 = repmat((0:Np-1),Np,1);
rf2 = 1./(1+j2pi_tau_df*Nps*(K3-K4)); % Eq. (6.17a)
Rhp = rf;
Rpp = rf2 + eye(length(H_tilde),length(H_tilde))/snr; % Eq. (6.14)
H_MMSE = transpose(Rhp*inv(Rpp)*H_tilde.'); % MMSE estimate Eq. (6.15)
end