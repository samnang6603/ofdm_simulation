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
clear

%% Simulation parameter
rng('shuffle');
SIM.EbN0 = 15;
SIM.SNR = SIM.EbN0 + 10*log10(sqrt(10));
SIM.Fading = true;
SIM.AWGN = true;
SIM.ChannelEstimation = false;
if ~SIM.Fading % if no fading, don't estimate
    SIM.ChannelEstimation = false;
end
SIM.FECToggle = 1;
SIM.Interleave = 1;

%% Data parameter
main_path = fileparts(cd);
path_parts = strsplit(main_path,filesep);
main_folder_ind = find(strcmp(path_parts,'ofdm_simulation'));
if length(path_parts) > main_folder_ind
    path_parts(end:main_folder_ind+1) = [];
end
data_file_path = fullfile(strjoin(path_parts,filesep),'data','cameraman2.tif');
DATA.Data = imread(data_file_path);
DATA = image_bitencode(DATA);

%% OFDM parameter
OFDM.M = 16;
OFDM.BitsPerSymbol = log2(OFDM.M);
OFDM.NumPilotSpacing = 8;
OFDM.NumBits = DATA.TotalBits;
OFDM.NumBitsPerFrame = 256*6;
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
FEC = fec_init('ldpc',SIM);


%% Start Sim
for frame = 1:OFDM.NumFrames
    [OFDM,FEC] = ofdm_transmit(DATA,OFDM,SIM,FEC,frame);
    [OFDM,CHANNEL] = channel_apply(CHANNEL,OFDM,SIM);
    [OFDM,DATA] = ofdm_receive(DATA,OFDM,CHANNEL,SIM,FEC,frame);
end
DATA = image_bitdecode(DATA);
error_percentage = 100*sum(DATA.BitData == DATA.BitDataHat)/length(DATA.BitData);
ber = sum(DATA.BitData ~= DATA.BitDataHat)/length(DATA.BitData);
fprintf('Correct bit recovered = %.6f%% \n',error_percentage);
fprintf('BER = %.6f \n',ber);
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
DATA.BitData = zeros(DATA.TotalBits,1,'int8');
for m = 0:DATA.NumElement-1
    str = sprintf(' %s',cbin{m+1,:});
    bit_num = sscanf(str,'%d');
    DATA.BitData(8*m+1:8*m+8) =  bit_num;
end
DATA.BitDataHat = zeros(DATA.TotalBits,1,'int8'); % placeholder for decoded data
DATA.DataHat = zeros(DATA.NumElement,1,'uint8');
end

function DATA = image_bitdecode(DATA)
% decode grayscale
data_hat_reshape = reshape(DATA.BitDataHat,8,DATA.TotalBits/8).';
data_hat_decode = mat2str(data_hat_reshape);
data_hat_decode([1 end]) = []; % remove brackets caused by mat2str()
str_len = 16; % num of char including bit and delimiter
for m = 0:DATA.NumElement-1
    si = str_len*m + 1;
    ei = str_len*m + str_len-1;
    DATA.DataHat(m+1) = uint8(bin2dec(data_hat_decode(si:ei)));
end
DATA.DataHat = reshape(DATA.DataHat,DATA.Size);
end

function [OFDM,FEC] = ofdm_transmit(DATA,OFDM,SIM,FEC,frame)
OFDM.data_frame = DATA.BitData((frame-1)*OFDM.NumBitsPerFrame+1:...
                  (frame-1)*OFDM.NumBitsPerFrame+OFDM.NumBitsPerFrame);
if SIM.Interleave
    FEC.IntrlvElem = randperm(length(OFDM.data_frame));
    OFDM.data_frame = intrlv(OFDM.data_frame,FEC.IntrlvElem);
end
OFDM.data_frame_codeword = OFDM.data_frame;
if SIM.FECToggle
    [OFDM.data_frame_codeword,FEC] = fec_encode(OFDM.data_frame,FEC);
    % Update OFDM parameter based on FEC parameter
    OFDM.NumBitsPerFrame = length(OFDM.data_frame_codeword); 
    OFDM.NumCarriersPerFrame = OFDM.NumBitsPerFrame/OFDM.BitsPerSymbol;
    OFDM.NumPilotPerFrame = OFDM.NumCarriersPerFrame/OFDM.NumPilotSpacing;
    OFDM.NumCarrierPilotPerFrame = OFDM.NumCarriersPerFrame + OFDM.NumPilotPerFrame;
    OFDM.NumCyclicSymsPerFrame = floor(OFDM.NumCarriersPerFrame*0.25);
    OFDM.NumCyclicPilotSymsPerFrame = floor(OFDM.NumCarrierPilotPerFrame*0.25);
end

[OFDM.data_frame_codeword,OFDM] = ofdm_transmit_validate_mapper_input(...
                        OFDM.data_frame_codeword,OFDM);

OFDM.TxSymbols = qammod(OFDM.data_frame_codeword,OFDM.M,InputType='bit');

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
    %CHANNEL.ImpulseResponse = jake(120,CHANNEL.TapLength);
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

OFDM.RxDemod = qamdemod(OFDM.RxFFT,OFDM.M,OutputType=FEC.DemapOutputType);

% Remove any padded zeros by Mapper at Tx
OFDM.RxDemod(end-OFDM.MapperZeroPadded+1:end) = [];

if SIM.FECToggle
    % OFDM.RxDemod = vitdec(OFDM.RxDemod,FEC.Trellis,FEC.VitDec.TraceBackDepth,...
    %                       FEC.VitDec.OpMode,FEC.VitDec.DecType);
    OFDM.RxDemod = fec_decode(OFDM.RxDemod,FEC);
    % Re-update OFDM param after decoding
    OFDM.NumBitsPerFrame = length(OFDM.data_frame);
end
if SIM.Interleave
    OFDM.RxDemod = deintrlv(OFDM.RxDemod,FEC.IntrlvElem);
end

DATA.BitDataHat((frame-1)*OFDM.NumBitsPerFrame+1:...
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

function FEC = fec_init(fec_type,SIM)
if ~SIM.FECToggle
    FEC.DemapOutputType = 'bit';
    return
end
switch lower(fec_type)
    case 'ldpc'
        FEC.Type = 'LDPC';
        FEC.P = [
    16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
    25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
    25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
     9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
    24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
     2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0
    ];
        FEC.Z = 27;
        FEC.PCMatrix = ldpcQuasiCyclicMatrix(FEC.Z,FEC.P);
        FEC.Ecfg = ldpcEncoderConfig(FEC.PCMatrix);
        FEC.Dcfg = ldpcDecoderConfig(FEC.PCMatrix,'layered-bp');
        FEC.DemapOutputType = 'approxllr';
        FEC.Maxnumiter = 10;
    case {'conv','convolutional'}
        FEC.Type = 'Convolutional';
        FEC.Constraint = 7;
        FEC.Trellis = poly2trellis(FEC.Constraint,[171 133]);
        FEC.VitDec.TraceBackDepth = 5*(FEC.Constraint-1); % approx for constraint
        FEC.VitDec.OpMode = 'trunc';
        FEC.VitDec.DecType = 'hard';
        FEC.DemapOutputType = 'bit';
    otherwise
        error('Unsupported FEC')
end
end

function [cw,FEC] = fec_encode(x,FEC)
switch FEC.Type
    case 'Convolutional'
        cw = convenc(x,FEC.Trellis);
    case 'LDPC'
        % Data segmentation to make LDPC input valid
        [x,FEC] = fec_encode_ldpc_segment(x,FEC); 
        nseg = FEC.NumSegment; % LDPC number of segment
        ninfo = FEC.Ecfg.NumInformationBits; % LDPC input bit lenght
        cwinfo = FEC.Ecfg.BlockLength; % LDPC output codeword length
        cw = zeros(FEC.SegCodeWordLength,1,'like',x); % allocate codeword
        for k = 1:nseg
            tmp_idx_x = ninfo*(k-1)+1:ninfo*(k-1) + ninfo;
            tmp_idx_cw = cwinfo*(k-1)+1:cwinfo*(k-1) + cwinfo;
            cw(tmp_idx_cw) = ldpcEncode(x(tmp_idx_x),FEC.Ecfg);
        end
end
end

function x = fec_decode(cw,FEC)
switch FEC.Type
    case 'Convolutional'
        x = vitdec(cw,FEC.Trellis,FEC.VitDec.TraceBackDepth,...
                          FEC.VitDec.OpMode,FEC.VitDec.DecType);
    case 'LDPC'
        % Decoder Reverse the data segmentation
        llr = cw; % variable naming to denote input to LDPC is LLR
        nseg = FEC.NumSegment;
        ninfo = FEC.Ecfg.NumInformationBits;
        cwinfo = FEC.Ecfg.BlockLength;
        x = zeros(FEC.SegInputLength,1,'int8');
        for k = 1:nseg
            tmp_idx_x = ninfo*(k-1)+1:ninfo*(k-1) + ninfo;
            tmp_idx_cw = cwinfo*(k-1)+1:cwinfo*(k-1) + cwinfo;
            x(tmp_idx_x) = ldpcDecode(llr(tmp_idx_cw),FEC.Dcfg,FEC.Maxnumiter);
        end
        % Remove paddded zeros
        x(end-FEC.NumZeroPadded+1:end) = [];
end
end

function [xseg,FEC] = fec_encode_ldpc_segment(x,FEC)
data_len = length(x);
ldpc_input_len = FEC.Ecfg.NumInformationBits;
if data_len == ldpc_input_len
    xseg = x;
    FEC.NumSegment = 1;
    return
end
% get remainder length
remainder_len = rem(data_len, ldpc_input_len); 

% find the length that is an integer multiple of fit length
fit_len = data_len - remainder_len; 

% Find number of full segments (segment that don't need padding)
num_full_segments = fit_len/ldpc_input_len;

% Total number of segment is number of full segments + 1 (segment that
% needs padding)
num_segments = num_full_segments + 1;

% Find number of zero pad for the segment that needs padding
num_zero_pad = ldpc_input_len - remainder_len;

% concatenate the zeros
xseg = [x
        zeros(num_zero_pad,1,'like',x)];

FEC.NumSegment = num_segments;

% Final length of LDPC input including all segment
FEC.SegInputLength = num_segments*ldpc_input_len; 

% Final length of codeword output after padding
FEC.SegCodeWordLength = num_segments*FEC.Ecfg.BlockLength;

% Save the number of zeros being padded so receiver can remove
FEC.NumZeroPadded = num_zero_pad;

end

function [y,OFDM] = ofdm_transmit_validate_mapper_input(x,OFDM)
% validate and adjust the length of output y to see if valid for mapper
% input.
% Then Find next value of the length of x that is divisible by n then pad zero
% by num_zero_pad to that value to output y
%
l = length(x);
n = OFDM.BitsPerSymbol;
mapper_len_cond = rem(length(x),n) ~= 0;
if mapper_len_cond
    nextval = l + (n - rem(l,n));
    num_zero_pad = nextval-l;
    y = [x(:);zeros(num_zero_pad,1)];
    OFDM.MapperZeroPadded = num_zero_pad;
else
    y = x;
    OFDM.MapperZeroPadded = 0;
end
end