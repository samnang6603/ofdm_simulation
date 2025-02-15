%% CFO correction testbench

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
NumBitsPerFrame=256*1;
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

%% Preamble
NumPreambleSpacing = 4;
NumPreamblePerFrame = NumCarriersPerFrame/NumPreambleSpacing;
NumCarriersPreamblePerFrame = NumPreamblePerFrame + NumCarriersPerFrame;
NumCyclicSymsPreamblePerFrame= ceil(NumCarriersPreamblePerFrame*0.25);

%% CFO
cfo = 0.15;

%% 
for frame = 1:NumFrames
    data_frame = data((frame-1)*NumBitsPerFrame+1:(frame-1)*NumBitsPerFrame+NumBitsPerFrame);
    
    %data_frame_codeword = convenc(data_frame,trellis);
    data_frame_codeword = data_frame;

    % modulation
    y = qammod(data_frame_codeword,NumQam,InputType='bit');

    % add preamble for CFO correction
    cpf = NumCarriersPreamblePerFrame;
    spc = NumPreambleSpacing;
    np = NumPreamblePerFrame;
    yp = zeros(cpf,1);
    % Generate preamble using quadratic phase sequence
    Xp = exp(1i*pi*(0:NumPreamblePerFrame-1).^2/np); 
    pream_loc = 1:spc+1:cpf;
    yp(pream_loc) = Xp;
    for b = 0:NumPreamblePerFrame-1
        yp(2+b*(spc+1):1+b*(spc+1)+spc) = y(b*spc+1:b*spc+spc);
    end
    
    % IFFT and add cyclic prefix
    nfft = length(y); %NumCarriersPerFrame + m*NumBitsPerFrame/k;
    NumCyclicSymsPerFrameEncode = floor(nfft*0.25); 
    ifft_sig = ifft(yp);
    cyclic_idx = nfft-NumCyclicSymsPerFrameEncode+1:nfft;
    cext_data = [ifft_sig(cyclic_idx); ifft_sig];

    % AWGN at Receiver
    ofdm_sig = awgn(cext_data,SNR,'measured');

    % Add CFO and correct
    ofdm_sig_cfo = ofdm_sig.*exp(1i*2*pi*cfo*(0:length(ofdm_sig)-1)'/nfft);

    % CFO correction
    cfo_est = cfo_cp(ofdm_sig_cfo,nfft,NumCyclicSymsPerFrameEncode);
    ofdm_sig_cfo_corr = ofdm_sig_cfo.*exp(-1i*2*pi*cfo_est*(0:length(ofdm_sig)-1)'/nfft);

    % cyclic prefix remove
    cext_rem = ofdm_sig;
    cext_rem_cfo = ofdm_sig_cfo;
    cext_rem_cfo_corr = ofdm_sig_cfo_corr;
    cext_rem(1:NumCyclicSymsPerFrameEncode) = [];
    cext_rem_cfo(1:NumCyclicSymsPerFrameEncode) = [];
    cext_rem_cfo_corr(1:NumCyclicSymsPerFrameEncode) = [];

    % FFT
    fft_sig = fft(cext_rem);
    fft_sig_cfo = fft(cext_rem_cfo);
    fft_sig_cfo_corr = fft(cext_rem_cfo_corr);

    fft_sig(pream_loc) = [];
    fft_sig_cfo(pream_loc) = [];
    fft_sig_cfo_corr(pream_loc) = [];

    % demodulation
    data_hat_codeword = qamdemod(fft_sig_cfo_corr,NumQam,OutputType='bit');

    data_hat_frame = data_hat_codeword;

    %data_hat_frame = vitdec(data_hat_codeword,trellis,tbdepth,opmode,dectype);

    % reshape for decode
    data_hat((frame-1)*NumBitsPerFrame+1:(frame-1)*NumBitsPerFrame+NumBitsPerFrame) = data_hat_frame;
end
BER = sum(data~=data_hat)/total_bits;

function cfo_est = cfo_cp(x,nfft,ng)
x = x.';
nn = 1:ng;
cfo_est = angle(x(nn+nfft)*x(nn)')/(2*pi);
end