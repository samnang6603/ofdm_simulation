%% LDPC vs No LDPC BER
close all; clear; clc

rng shuffle  % fix random seed

% Modulation 16-QAM
NumQam=16;
K=log2(NumQam);
EbN0 = 0:3:21;
SNR = EbN0 + 10*log10(sqrt(10));    % Convert EbN0 to SNR
SNR = 0:3:30;

%% LDPC Parameters
% Using code from MATLAB example
P = [
    16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
    25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
    25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
     9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
    24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
     2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0
    ];
blockSize = 27;
pcmatrix = ldpcQuasiCyclicMatrix(blockSize,P);
ecfg = ldpcEncoderConfig(pcmatrix);
dcfg_bp = ldpcDecoderConfig(pcmatrix);
dcfg_lbp = ldpcDecoderConfig(pcmatrix,'layered-bp');
dcfg_nms = ldpcDecoderConfig(pcmatrix,'norm-min-sum');
dcfg_oms = ldpcDecoderConfig(pcmatrix,'offset-min-sum');
maxnumiter = 10;

%% Data Parameters
total_bits = 100*ecfg.NumInformationBits;
data = randi([0,1],total_bits,1);
data_hat = zeros(total_bits,1);
data_hat_ldpc_bp = zeros(total_bits,1);
data_hat_ldpc_lbp = zeros(total_bits,1);
data_hat_ldpc_nms = zeros(total_bits,1);
data_hat_ldpc_oms = zeros(total_bits,1);

%% OFDM Parameters
NumBits = total_bits;
NumBitsPerFrame=ecfg.NumInformationBits;
NumFrames = NumBits/NumBitsPerFrame;
NumCarriersPerFrame=floor(NumBitsPerFrame/K);
NumCyclicSymsPerFrame=floor(NumCarriersPerFrame*0.25);
NumSyms = NumBits/K;
NumSymsCyclic = NumCarriersPerFrame + NumCyclicSymsPerFrame;

%% Monte Carlo Parameters
monte_carlo_iter = 5;
BER_mc = zeros(monte_carlo_iter,1);
BER_mc_ldpc_bp = BER_mc;
BER_mc_ldpc_lbp = BER_mc;
BER_mc_ldpc_nms = BER_mc;
BER_mc_ldpc_oms = BER_mc;

%% BER for each SNR allocation
BER = zeros(length(SNR),1);
BER_ldpc_bp = BER;
BER_ldpc_lbp = BER;
BER_ldpc_nms = BER;
BER_ldpc_oms = BER;

%% Channel
Chann_tap = 0;

%% Loop
for snr_idx = 1:length(SNR)
    for mc_idx = 1:monte_carlo_iter
        for frame = 1:NumFrames
            data_frame = data((frame-1)*NumBitsPerFrame+1:(frame-1)*NumBitsPerFrame+NumBitsPerFrame);

            data_frame_ldpc = ldpcEncode(data_frame,ecfg);

            mapper_len_cond = rem(length(data_frame),K) ~= 0;
            mapper_len_cond_ldpc = rem(length(data_frame_ldpc),K) ~= 0;

            if mapper_len_cond % ensure proper length of qam for no LDPC case
                [data_frame,num_zero_pad] = ...
                    nextdivpadzero(data_frame,K);
            end
            if mapper_len_cond_ldpc % ensure proper length of qam for LDPC case
                [data_frame_ldpc,num_zero_pad_ldpc] = ...
                    nextdivpadzero(data_frame_ldpc,K);
            end

            % modulation
            y = qammod(data_frame,NumQam,InputType='bit',UnitAveragePower=0);
            y_ldpc = qammod(data_frame_ldpc,NumQam,InputType='bit',UnitAveragePower=0);

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
            % ray_fading = (randn(Chann_tap,1) + 1i*randn(Chann_tap,1))/sqrt(2);
            % cext_data = filter(ray_fading,1,cext_data);
            % cext_data_ldpc = filter(ray_fading,1,cext_data_ldpc);

            % AWGN at Receiver for no LDPC
            [ofdm_sig,noisevar] = awgn(cext_data,SNR(snr_idx),"measured");
            cext_rem = ofdm_sig;
            cext_rem(1:NumCyclicSymsPerFrame) = [];

            % AWGN at Receiveer for LDPC
            [ofdm_sig,noisevar_ldpc] = awgn(cext_data_ldpc,SNR(snr_idx),"measured");
            cext_rem_ldpc = ofdm_sig;
            cext_rem_ldpc(1:NumCyclicSymsPerFrame_ldpc) = [];

            % FFT
            fft_sig = fft(cext_rem);
            fft_sig_ldpc = fft(cext_rem_ldpc);

            % demodulation no LDPC
            y_hat = qamdemod(fft_sig(:),...
                NumQam,OutputType='bit',UnitAveragePower=0);
            % demodulation LDPC
            y_hat_ldpc = qamdemod(fft_sig_ldpc(:),...
                NumQam,OutputType='approxllr',UnitAveragePower=0);

            % y_hat = qamdemod(awgn(y,SNR(snr_idx),'measured'),NumQam,...
            %     OutputType="bit",UnitAveragePower=0);
            % y_hat_ldpc = qamdemod(awgn(y_ldpc,SNR(snr_idx),'measured'),NumQam,...
            %     OutputType="approxllr",UnitAveragePower=0);

            % if additional zero padded to make length proper for
            % constellation mapper then remove the extra padded zeros
            if mapper_len_cond
                y_hat(end-num_zero_pad+1:end) = [];
            end
            if mapper_len_cond_ldpc
                y_hat_ldpc(end-num_zero_pad_ldpc+1:end) = [];
            end

            data_hat_frame = y_hat;

            data_hat_frame_ldpc_bp = ldpcDecode(y_hat_ldpc,dcfg_bp,maxnumiter);
            data_hat_frame_ldpc_lbp = ldpcDecode(y_hat_ldpc,dcfg_lbp,maxnumiter);
            data_hat_frame_ldpc_nms = ldpcDecode(y_hat_ldpc,dcfg_nms,maxnumiter);
            data_hat_frame_ldpc_oms = ldpcDecode(y_hat_ldpc,dcfg_oms,maxnumiter);


            % reshape for decode
            idx = (frame-1)*NumBitsPerFrame+1:(frame-1)*NumBitsPerFrame+NumBitsPerFrame;
            data_hat(idx) = data_hat_frame;
            data_hat_ldpc_bp(idx) = data_hat_frame_ldpc_bp;
            data_hat_ldpc_lbp(idx) = data_hat_frame_ldpc_lbp;
            data_hat_ldpc_nms(idx) = data_hat_frame_ldpc_nms;
            data_hat_ldpc_oms(idx) = data_hat_frame_ldpc_oms;

        end
        BER_mc(mc_idx) = sum(data~=data_hat)/total_bits;
        BER_mc_ldpc_bp(mc_idx) = sum(data~=data_hat_ldpc_bp)/total_bits;
        BER_mc_ldpc_lbp(mc_idx) = sum(data~=data_hat_ldpc_lbp)/total_bits;
        BER_mc_ldpc_nms(mc_idx) = sum(data~=data_hat_ldpc_nms)/total_bits;
        BER_mc_ldpc_oms(mc_idx) = sum(data~=data_hat_ldpc_oms)/total_bits;
    end
    BER(snr_idx) = mean(BER_mc);
    BER_ldpc_bp(snr_idx) = mean(BER_mc_ldpc_bp);
    BER_ldpc_lbp(snr_idx) = mean(BER_mc_ldpc_lbp);
    BER_ldpc_nms(snr_idx) = mean(BER_mc_ldpc_nms);
    BER_ldpc_oms(snr_idx) = mean(BER_mc_ldpc_oms);
end
figure
semilogy(SNR,BER,'r-*'), hold on, grid on
semilogy(SNR,BER_ldpc_bp,'b-^')
semilogy(SNR,BER_ldpc_lbp,'g-*')
semilogy(SNR,BER_ldpc_nms,'k-o')
semilogy(SNR,BER_ldpc_oms,'m-v')
legend('No LDPC','LDPC-BP','LDPC-LBP','LDPC-NMS','LDPC-OMS')
SNR = SNR(:);
compare = [SNR,BER,BER_ldpc_bp,BER_ldpc_lbp,BER_ldpc_nms,BER_ldpc_oms];
T = table(SNR,BER,BER_ldpc_bp,BER_ldpc_lbp,BER_ldpc_nms,BER_ldpc_oms);
disp(T)


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
