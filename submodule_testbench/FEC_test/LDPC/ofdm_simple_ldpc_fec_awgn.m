%% LDPC/Convolutional code testbench

close all; clear; clc

rng shuffle  % fix random seed

% Modulation 16-QAM
NumQam=16;
K=log2(NumQam);
EbN0 = 5;
SNR = EbN0 + 10*log10(sqrt(10));    % Convert EbN0 to SNR

%% LDPC parameter
if ~isdeployed
    addpath('./vodafone-chair-5g-nr-ldpc-master/codes');
end
blksize = 256;
coderate = '1/2';
LDPC = ldpcGet(blksize,coderate);
LDPCToggle = 0;
if LDPCToggle
    qamdemodout = 'approxllr';
else
    qamdemodout = 'bit';
end
Hsp = sparse(logical(LDPC.H));
ecfg = ldpcEncoderConfig(Hsp);
dcfg = ldpcDecoderConfig(Hsp);

%% Data
total_bits = 10000*LDPC.numInfBits;
data = randi([0 1],total_bits,1);
data_hat = zeros(total_bits,1);

%% OFDM Param
NumBits = total_bits; 
NumBitsPerFrame=LDPC.numInfBits;
NumFrames = NumBits/NumBitsPerFrame;
NumCarriersPerFrame=floor(NumBitsPerFrame/K);
NumCyclicSymsPerFrame=floor(NumCarriersPerFrame*0.25);
NumSyms = NumBits/K;
NumSymsCyclic = NumCarriersPerFrame + NumCyclicSymsPerFrame;


%% Trellis Definition and Convolutional Encode Parameter
constraintL = 7;
trellis = poly2trellis(7,[171 133]); % Define trellis
tbdepth = 5*(constraintL-1); % coderate = 1/2;
opmode = 'trunc';
dectype = 'hard';


%% 
tic
for frame = 1:NumFrames
    data_frame = data((frame-1)*NumBitsPerFrame+1:(frame-1)*NumBitsPerFrame+NumBitsPerFrame);
    
    %data_frame_codeword = convenc(data_frame,trellis);
    if LDPCToggle
        %data_frame_codeword = ldpcEncode1(data_frame',LDPC);
        data_frame_codeword = ldpcEncode(data_frame,ecfg);
    else
        data_frame_codeword = data_frame;
    end

    mapper_len_cond = rem(length(data_frame_codeword),K) ~= 0;

    if mapper_len_cond % ensure proper length of qam
        [data_frame_codeword,num_zero_pad] = ...
            nextdivpadzero(data_frame_codeword,K);
    end

    % modulation
    y = qammod(data_frame_codeword(:),NumQam,InputType='bit',UnitAveragePower=1);

    % IFFT and add cyclic prefix
    nfft = length(y); %NumCarriersPerFrame + m*NumBitsPerFrame/k;
    NumCyclicSymsPerFrameEncode = floor(nfft*0.25); 
    ifft_sig = ifft(y);
    cyclic_idx = nfft-NumCyclicSymsPerFrameEncode+1:nfft;
    cext_data = [ifft_sig(cyclic_idx); ifft_sig];

    % AWGN at Receiver
    ofdm_sig = awgn(cext_data,SNR,'measured');
    cext_rem = ofdm_sig;
    cext_rem(1:NumCyclicSymsPerFrameEncode) = [];

    % FFT
    fft_sig = fft(cext_rem);

    data_hat_codeword = qamdemod(fft_sig(:),...
        NumQam,OutputType=qamdemodout,UnitAveragePower=1);

    % if additional zero pad to make length proper for constellation mapper
    % then remove the extra padded zeros
    if mapper_len_cond 
        data_hat_codeword(end-num_zero_pad+1:end) = [];
    end

    % demodulation
    if LDPCToggle
        %data_hat_frame = ldpcDecode1(data_hat_codeword',LDPC);
        data_hat_frame = ldpcDecode(data_hat_codeword,dcfg,10);
    else
        data_hat_frame = data_hat_codeword;
    end

    % reshape for decode
    data_hat((frame-1)*NumBitsPerFrame+1:(frame-1)*NumBitsPerFrame+NumBitsPerFrame) = data_hat_frame;
end
BER = sum(data~=data_hat)/total_bits
toc

function [y,num_zero_pad] = nextdivpadzero(x,n)
% find next value of the length of x that is divisible by n then pad zero
% by num_zero_pad to that value to output y
% 
l = length(x);
nextval = l + (n - rem(l,n));
num_zero_pad = nextval-l;
y = [x(:);zeros(num_zero_pad,1)];
end
