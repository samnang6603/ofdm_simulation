%% FEC Convolutional code testbench

close all; clear; clc

rng shuffle  % fix random seed

% Modulation 16-QAM
NumQam=16;
K=log2(NumQam);
EbN0 = 10;
SNR = EbN0 + 10*log10(sqrt(10));    % Convert EbN0 to SNR

%% Data
total_bits = 256*40;
data = randi([0 1],total_bits,1);
data_hat = zeros(total_bits,1);

%% Generating and data
NumBits = total_bits; % # bits. Can modify this, but has to be power of two
NumBitsPerFrame=256*4;
NumFrames = NumBits/NumBitsPerFrame;
NumCarriersPerFrame=NumBitsPerFrame/K;
NumCyclicSymsPerFrame=NumCarriersPerFrame*0.25;
NumSyms = NumBits/K;
NumSymsCyclic = NumCarriersPerFrame + NumCyclicSymsPerFrame;

%% Trellis Definition and Convolutional Encode Parameter
constraintL = 7;
trellis = poly2trellis(7,[171 133]); % Define trellis
tbdepth = 5*(constraintL-1); % coderate = 1/2;
opmode = 'trunc';
dectype = 'hard';


%% 
for frame = 1:NumFrames
    data_frame = data((frame-1)*NumBitsPerFrame+1:(frame-1)*NumBitsPerFrame+NumBitsPerFrame);
    
    data_frame_codeword = convenc(data_frame,trellis);

    % modulation
    y = qammod(data_frame_codeword,NumQam,InputType='bit',UnitAveragePower=1);

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

    % demodulation
    data_hat_codeword = qamdemod(fft_sig,NumQam,OutputType='bit',UnitAveragePower=1);

    data_hat_frame = vitdec(data_hat_codeword,trellis,tbdepth,opmode,dectype);

    % reshape for decode
    data_hat((frame-1)*NumBitsPerFrame+1:(frame-1)*NumBitsPerFrame+NumBitsPerFrame) = data_hat_frame;
end
BER = sum(data~=data_hat)/total_bits;