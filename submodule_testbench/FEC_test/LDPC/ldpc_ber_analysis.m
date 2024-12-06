%% LDPC vs No LDPC BER
close all; clear; clc

rng shuffle  % fix random seed

% Modulation 16-QAM
NumQam=16;
K=log2(NumQam);
EbN0 = 0:3:21;
SNR = EbN0 + 10*log10(sqrt(10));    % Convert EbN0 to SNR

%% LDPC Parameters
if ~isdeployed
    addpath('./vodafone-chair-5g-nr-ldpc-master/codes');
end
blksize = 256;
coderate = '1/2';
LDPC = ldpcGet(blksize,coderate);
Hsp = sparse(logical(LDPC.H));
ecfg = ldpcEncoderConfig(Hsp);
dcfg = ldpcDecoderConfig(Hsp);
numiter = 15;

%% Data Parameters
total_bits = 1000*LDPC.numInfBits;
data = randi([0,1],total_bits,1);
data_hat = zeros(total_bits,1);
data_hat_ldpc = zeros(total_bits,1);

%% OFDM Parameters
NumBits = total_bits;
NumBitsPerFrame=LDPC.numInfBits;
NumFrames = NumBits/NumBitsPerFrame;
NumCarriersPerFrame=floor(NumBitsPerFrame/K);
NumCyclicSymsPerFrame=floor(NumCarriersPerFrame*0.25);
NumSyms = NumBits/K;
NumSymsCyclic = NumCarriersPerFrame + NumCyclicSymsPerFrame;

%% Monte Carlo Parameters
monte_carlo_iter = 5;
BER_mc = zeros(monte_carlo_iter,1);
BER_mc_ldpc = BER_mc;

%% BER for each SNR allocation
BER = zeros(length(SNR),1);
BER_ldpc = BER;

%% Channel
Chann_tap = 2;

%% Loop
for snr_idx = 1:length(SNR)
    for mc_idx = 1:monte_carlo_iter
        for frame = 1:NumFrames
            data_frame = data((frame-1)*NumBitsPerFrame+1:(frame-1)*NumBitsPerFrame+NumBitsPerFrame);

            data_frame_ldpc = ldpcEncode(data_frame,ecfg);

            mapper_len_cond = rem(length(data_frame_ldpc),K) ~= 0;

            if mapper_len_cond % ensure proper length of qam
                [data_frame_ldpc,num_zero_pad] = ...
                    nextdivpadzero(data_frame_ldpc,K);
            end

            % modulation
            y = qammod(data_frame,NumQam,InputType='bit',UnitAveragePower=1);
            y_ldpc = qammod(data_frame_ldpc,NumQam,InputType='bit',UnitAveragePower=1);

            % IFFT and add cyclic prefix for no LDPC
            nfft = length(y); %NumCarriersPerFrame + m*NumBitsPerFrame/k;
            NumCyclicSymsPerFrame = floor(nfft*0.25);
            ifft_sig = ifft(y);
            cyclic_idx = nfft-NumCyclicSymsPerFrame+1:nfft;
            cext_data = [ifft_sig(cyclic_idx); ifft_sig];

            % IFFT and add cyclic prefix for LDPC
            nfft = length(y_ldpc); %NumCarriersPerFrame + m*NumBitsPerFrame/k;
            NumCyclicSymsPerFrame_ldpc = floor(nfft*0.25);
            ifft_sig = ifft(y_ldpc);
            cyclic_idx = nfft-NumCyclicSymsPerFrame_ldpc+1:nfft;
            cext_data_ldpc = [ifft_sig(cyclic_idx); ifft_sig];

            % Fading Rayleigh
            ray_fading = (randn(Chann_tap,1) + 1i*randn(Chann_tap,1))/sqrt(2);
            cext_data = filter(ray_fading,1,cext_data);
            cext_data_ldpc = filter(ray_fading,1,cext_data_ldpc);

            % AWGN at Receiver for no LDPC
            ofdm_sig = awgn(cext_data,SNR(snr_idx),'measured');
            cext_rem = ofdm_sig;
            cext_rem(1:NumCyclicSymsPerFrame) = [];

            % AWGN at Receiveer for LDPC
            ofdm_sig = awgn(cext_data_ldpc,SNR(snr_idx),'measured');
            cext_rem_ldpc = ofdm_sig;
            cext_rem_ldpc(1:NumCyclicSymsPerFrame_ldpc) = [];

            % FFT
            fft_sig = fft(cext_rem);
            fft_sig_ldpc = fft(cext_rem_ldpc);

            % demodulation no LDPC
            y_hat = qamdemod(fft_sig(:),...
                NumQam,OutputType='bit',UnitAveragePower=1);
            % demodulation LDPC
            y_hat_ldpc = qamdemod(fft_sig_ldpc(:),...
                NumQam,OutputType='approxllr',UnitAveragePower=1);

            % if additional zero padded to make length proper for 
            % constellation mapper then remove the extra padded zeros
            if mapper_len_cond
                y_hat_ldpc(end-num_zero_pad+1:end) = [];
            end
            
            data_hat_frame = y_hat;
            data_hat_frame_ldpc = ldpcDecode(y_hat_ldpc,dcfg,numiter);

            % reshape for decode
            idx = (frame-1)*NumBitsPerFrame+1:(frame-1)*NumBitsPerFrame+NumBitsPerFrame;
            data_hat(idx) = data_hat_frame;
            data_hat_ldpc(idx) = data_hat_frame_ldpc;
        end
        BER_mc(mc_idx) = sum(data~=data_hat)/total_bits;
        BER_mc_ldpc(mc_idx) = sum(data~=data_hat_ldpc)/total_bits;
    end
    BER(snr_idx) = mean(BER_mc);
    BER_ldpc(snr_idx) = mean(BER_mc_ldpc);
end
figure
semilogy(SNR,BER,'r-*'), hold on, grid on
semilogy(SNR,BER_ldpc,'b-^')
legend('No LDPC','LDPC')


%% Local fcn
function [y,num_zero_pad] = nextdivpadzero(x,n)
% find next value of the length of x that is divisible by n then pad zero
% by num_zero_pad to that value to output y
%
l = length(x);
nextval = l + (n - rem(l,n));
num_zero_pad = nextval-l;
y = [x(:);zeros(num_zero_pad,1)];
end
