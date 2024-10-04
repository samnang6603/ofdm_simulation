%% Simulation of OFDM Fading Channel with Channel Estimation and EQ
%{ 
    Simulation of OFDM with Fading Channel and Channel estimation using
    MMSE or LS with EQ at the Rx.

    No Global Variables
    Some Parameters are adjustable

    Author: Sam An
    Date: 09/01/2024
%}

%% 
clear, clc, close all

%% Simulation parameter
rng('shuffle');
SIM.EbN0 = 15;
SIM.SNR = SIM.EbN0 + 10*log10(sqrt(10));
SIM.Fading = true;
SIM.AWGN = true;
SIM.ChannelEstimation = true;
if ~SIM.Fading % if no fading, don't estimate
    SIM.ChannelEstimation = false;
end
SIM.FECToggle = 1;

%% Data parameter
DATA.Data = imread('cameraman2.tif');
DATA = image_bitencode(DATA);

%% OFDM parameter
OFDM.M = 16;
OFDM.BitsPerSymbol = log2(OFDM.M);
OFDM.NumPilotSpacing = 8;
OFDM.NumBits = DATA.TotalBits;
OFDM.NumBitsPerFrame = 256;
OFDM.NumFrames = OFDM.NumBits/OFDM.NumBitsPerFrame;
OFDM.NumCarriersPerFrame = OFDM.NumBitsPerFrame/OFDM.BitsPerSymbol;
OFDM.NumPilotPerFrame = OFDM.NumCarriersPerFrame/OFDM.NumPilotSpacing;
OFDM.NumCarrierPilotPerFrame = OFDM.NumCarriersPerFrame + OFDM.NumPilotPerFrame;
OFDM.NumCyclicSymsPerFrame = floor(OFDM.NumCarriersPerFrame*0.25);
OFDM.NumCyclicPilotSymsPerFrame = floor(OFDM.NumCarrierPilotPerFrame*0.25);

%% Channel parameter
CHANNEL.TapLength = 6;
CHANNEL.Delay = 0;
CHANNEL.DelayObj = dsp.Delay(CHANNEL.Delay);
CHANNEL.Equalization = true;
CHANNEL.EstimationType = 'ls';
if ~SIM.Fading % if no fading, don't do eq
    CHANNEL.Equalization = false;
end

%% FEC parameter
FEC.Constraint = 7;
FEC.Trellis = poly2trellis(FEC.Constraint,[171 133]);
FEC.VitDec.TraceBackDepth = 5*(FEC.Constraint-1); % approx for constraint
FEC.VitDec.OpMode = 'trunc';
FEC.VitDec.DecType = 'hard';

%% Start Sim
for frame = 1:OFDM.NumFrames
    OFDM = ofdm_transmit(DATA,OFDM,SIM,FEC,frame);
    [OFDM,CHANNEL] = channel_apply(CHANNEL,OFDM,SIM);
    [OFDM,DATA] = ofdm_receive(DATA,OFDM,CHANNEL,SIM,FEC,frame);
end
DATA = image_bitdecode(DATA,OFDM);
figure,
subplot(121)
imshow(DATA.Data), title('Original Image')
subplot(122)
imshow(DATA.DataHat), title('Received Image')

%% Local Functions
function DATA = image_bitencode(DATA)
% encode grayscale image to binary representation
x = DATA.Data;
DATA.Size = size(x);
DATA.NumElement = numel(x);
cbin = num2cell(dec2bin(x(:)));
DATA.TotalBits = 8*DATA.NumElement;
DATA.BitData = zeros(DATA.TotalBits,1,'uint8');
for m = 0:DATA.NumElement-1
    str = sprintf(' %s',cbin{m+1,:});
    bit_num = sscanf(str,'%d');
    DATA.BitData(8*m+1:8*m+8) =  bit_num;
end
DATA.DataHatVec = zeros(DATA.TotalBits,1,'uint8'); % placeholder for decoded data
DATA.DataHat = zeros(DATA.NumElement,1,'uint8');
end

function DATA = image_bitdecode(DATA,OFDM)
% decode grayscale
data_hat_reshape = reshape(DATA.DataHatVec,8,DATA.TotalBits/8).';
data_hat_decode = mat2str(data_hat_reshape);
data_hat_decode([1 end]) = [];
for m = 0:DATA.NumElement-1
    si = OFDM.M*m + 1;
    ei = OFDM.M*m + OFDM.M-1;
    DATA.DataHat(m+1) = uint8(bin2dec(data_hat_decode(si:ei)));
end
DATA.DataHat = reshape(DATA.DataHat,DATA.Size);
end

function OFDM = ofdm_transmit(DATA,OFDM,SIM,FEC,frame)
OFDM.data_frame = DATA.BitData((frame-1)*OFDM.NumBitsPerFrame+1:...
                  (frame-1)*OFDM.NumBitsPerFrame+OFDM.NumBitsPerFrame);
OFDM.data_frame_codeword = OFDM.data_frame;
if SIM.FECToggle
    OFDM.data_frame_codeword = convenc(OFDM.data_frame,FEC.Trellis);
    % Update OFDM parameter based on FEC parameter
    OFDM.NumBitsPerFrame = length(OFDM.data_frame_codeword); 
    OFDM.NumCarriersPerFrame = OFDM.NumBitsPerFrame/OFDM.BitsPerSymbol;
    OFDM.NumPilotPerFrame = OFDM.NumCarriersPerFrame/OFDM.NumPilotSpacing;
    OFDM.NumCarrierPilotPerFrame = OFDM.NumCarriersPerFrame + OFDM.NumPilotPerFrame;
    OFDM.NumCyclicSymsPerFrame = floor(OFDM.NumCarriersPerFrame*0.25);
    OFDM.NumCyclicPilotSymsPerFrame = floor(OFDM.NumCarrierPilotPerFrame*0.25);
end

OFDM.TxSymbols = qammod(OFDM.data_frame_codeword,OFDM.M,InputType='bit',UnitAveragePower=1);

if SIM.ChannelEstimation
    % Add pilot
    yp = zeros(OFDM.NumCarrierPilotPerFrame,1);
    Xp = 2*randi([0 1],OFDM.NumPilotPerFrame,1)-1;
    spc = OFDM.NumPilotSpacing;
    pilot_loc = 1:spc+1:OFDM.NumCarrierPilotPerFrame;
    yp(pilot_loc) = Xp;
    for k = 0:OFDM.NumPilotPerFrame-1
        yp(2+k*(spc+1):1+k*(spc+1)+spc) = OFDM.TxSymbols(k*spc+1:k*spc+spc);
    end
    % transmitter with pilot
    ifft_sig_pilot = ifft(yp,OFDM.NumCarrierPilotPerFrame);
    cyclic_idx = OFDM.NumCarrierPilotPerFrame-...
        OFDM.NumCyclicPilotSymsPerFrame+1:OFDM.NumCarrierPilotPerFrame;
    cext_data = [ifft_sig_pilot(cyclic_idx); ifft_sig_pilot];

    OFDM.PilotSignal = Xp;
    OFDM.PilotSignalLocation = pilot_loc;
else
    % transmitter no pilot
    ifft_sig = ifft(OFDM.TxSymbols,OFDM.NumCarriersPerFrame);
    cyclic_idx = OFDM.NumCarriersPerFrame-...
        OFDM.NumCyclicSymsPerFrame+1:OFDM.NumCarriersPerFrame;
    cext_data = [ifft_sig(cyclic_idx); ifft_sig];
end
OFDM.TxAir = cext_data;
end

function [OFDM,CHANNEL] = channel_apply(CHANNEL,OFDM,SIM)
if SIM.Fading
    CHANNEL.ImpulseResponse = jake(120,CHANNEL.TapLength);
    CHANNEL.ImpulseResponse = 1/sqrt(2)*(randn(6,1) + 1j*randn(6,1)); 
else
    CHANNEL.ImpulseResponse = 1;
end
OFDM.TxAirChannel = filter(CHANNEL.ImpulseResponse,1,OFDM.TxAir);
OFDM.TxAirChannel = CHANNEL.DelayObj(OFDM.TxAirChannel);
if SIM.AWGN
    OFDM.TxAirChannel = awgn(OFDM.TxAirChannel,SIM.SNR,'measured');
end
end

function [OFDM,CHANNEL] = channel_estimate(CHANNEL,OFDM,SIM)
switch lower(CHANNEL.EstimationType)
    case 'mmse'
        % Minimum-Mean-Squared-Error method
        CHANNEL.EstFreqResponse = chanest_mmse(OFDM.RxFFT,OFDM.PilotSignal,...
                                               OFDM.PilotSignalLocation,...
                                               OFDM.NumCarrierPilotPerFrame,...
                                               OFDM.NumPilotSpacing,...
                                               CHANNEL.ImpulseResponse.',SIM.SNR);
        CHANNEL = chanest_dft_enhance(CHANNEL,OFDM.NumCarrierPilotPerFrame);
    otherwise
        % Least-Square method default
        CHANNEL.EstFreqResponse = chanest_ls(OFDM.RxFFT,OFDM.PilotSignal,...
                                             OFDM.PilotSignalLocation,...
                                             OFDM.NumCarrierPilotPerFrame);
        CHANNEL = chanest_dft_enhance(CHANNEL,OFDM.NumCarrierPilotPerFrame);
end
OFDM.IdMatrixDiagLength = OFDM.NumCarrierPilotPerFrame;
end

function [OFDM,DATA] = ofdm_receive(DATA,OFDM,CHANNEL,SIM,FEC,frame)
cext_rem = OFDM.TxAirChannel;
if SIM.ChannelEstimation
    cext_rem(1:OFDM.NumCyclicPilotSymsPerFrame) = [];
    OFDM.RxFFT = fft(cext_rem);
    [OFDM,CHANNEL] = channel_estimate(CHANNEL,OFDM,SIM);
else
    % no estimation algorithm, no pilot signal
    % ideal, take the channel itself, no pilot signal
    cext_rem(1:OFDM.NumCyclicSymsPerFrame) = [];
    OFDM.RxFFT = fft(cext_rem);
    CHANNEL.EstFreqResponse = fft(CHANNEL.ImpulseResponse,OFDM.NumCarriersPerFrame);
    OFDM.IdMatrixDiagLength = OFDM.NumCarriersPerFrame;
end

if CHANNEL.Equalization
    [OFDM,CHANNEL] = channel_apply_eq(CHANNEL,OFDM,SIM);
end

OFDM.RxDemod = qamdemod(OFDM.RxFFT,OFDM.M,OutputType='bit',UnitAveragePower=1);

if SIM.FECToggle
    OFDM.RxDemod = vitdec(OFDM.RxDemod,FEC.Trellis,FEC.VitDec.TraceBackDepth,...
                          FEC.VitDec.OpMode,FEC.VitDec.DecType);
    % Re-update OFDM param after decoding
    OFDM.NumBitsPerFrame = length(OFDM.data_frame);
end

DATA.DataHatVec((frame-1)*OFDM.NumBitsPerFrame+1:...
    (frame-1)*OFDM.NumBitsPerFrame+OFDM.NumBitsPerFrame) = OFDM.RxDemod;
end

function [OFDM,CHANNEL] = channel_apply_eq(CHANNEL,OFDM,SIM)
H = diag(CHANNEL.EstFreqResponse);
CHANNEL.EqFreqResponse = (H'*H + (1/SIM.SNR*eye(OFDM.IdMatrixDiagLength)))\H';
OFDM.RxFFT = CHANNEL.EqFreqResponse*OFDM.RxFFT;
if SIM.ChannelEstimation % remove pilot signal if channel est being used
    OFDM.RxFFT(OFDM.PilotSignalLocation) = [];
end
end

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

function H_LS = chanest_ls(Y,Xp,pilot_loc,Nfft)
% reference: Cho MIMO-OFDM Wireless Communications with MATLAB
Np=length(pilot_loc); %Nfft/Nps; 
k=1:Np;
LS_est(k) = Y(pilot_loc(k))./Xp(k); % LS channel estimation
LS_est = LS_est(:);
% Linear/Spline interpolation
H_LS = interpolate1(LS_est,pilot_loc,Nfft);
end

function H_MMSE = chanest_mmse(Y,Xp,pilot_loc,Nfft,Nps,h,SNR)
% reference: Cho MIMO-OFDM Wireless Communications with MATLAB
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

function CHANNEL = chanest_dft_enhance(CHANNEL,nfft)
h_est = ifft(CHANNEL.EstFreqResponse);
h_est = h_est(1:CHANNEL.TapLength);
CHANNEL.EstFreqResponse = fft(h_est,nfft);
end


