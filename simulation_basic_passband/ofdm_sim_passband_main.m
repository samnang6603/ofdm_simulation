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
clear, clc

%% Simulation parameters
rng('shuffle');
SIM.EbN0 = 10;
SIM.SNR = SIM.EbN0 + 10*log10(sqrt(10));
SIM.Fading = true;
SIM.AWGN = true;
SIM.ChannelEstimation = true;
if ~SIM.Fading % if no fading, don't estimate
    SIM.ChannelEstimation = false;
end
SIM.FECToggle = 1;
SIM.Interleave = 1;

%% Data parameters
main_path = fileparts(cd);
path_parts = strsplit(main_path,filesep);
main_folder_ind = find(strcmp(path_parts,'ofdm_simulation'));
if length(path_parts) > main_folder_ind
    path_parts(end:main_folder_ind+1) = [];
end
data_file_path = fullfile(strjoin(path_parts,filesep),'data','cameraman2.tif');
DATA.Data = imread(data_file_path);
DATA = image_bitencode(DATA);

%% OFDM parameters
OFDM.M = 16;
OFDM.BitsPerSymbol = log2(OFDM.M);
OFDM.NumPilotSpacing = 8;
OFDM.NumBits = DATA.TotalBits;
OFDM.NumBitsPerFrame = 256*4;
OFDM.NumFrames = OFDM.NumBits/OFDM.NumBitsPerFrame;
OFDM.NumCarriersPerFrame = OFDM.NumBitsPerFrame/OFDM.BitsPerSymbol;
OFDM.NumPilotPerFrame = OFDM.NumCarriersPerFrame/OFDM.NumPilotSpacing;
OFDM.NumCarrierPilotPerFrame = OFDM.NumCarriersPerFrame + OFDM.NumPilotPerFrame;
OFDM.NumCyclicSymsPerFrame = floor(OFDM.NumCarriersPerFrame*0.25);
OFDM.NumCyclicPilotSymsPerFrame = floor(OFDM.NumCarrierPilotPerFrame*0.25);

%% Passband Upconverter parameters
RF.PassBandProcessingToggle = 1;
RF.CarrierFrequency = 100e3; % carrier frequency
RF.SamplingFrequency = 20*RF.CarrierFrequency; % sampling frequency
RF.SamplingPeriod = 1/RF.SamplingFrequency;

%% Channel parameters
CHANNEL.TapLength = 6;
CHANNEL.Delay = 0;
CHANNEL.DelayObj = dsp.Delay(CHANNEL.Delay);
CHANNEL.Equalization = true;
CHANNEL.EstimationType = 'ls';
if ~SIM.Fading % if no fading, don't do eq
    CHANNEL.Equalization = false;
end

%% FEC parameters
FEC.Type = 'Convolutional';
FEC.Constraint = 7;
FEC.Trellis = poly2trellis(FEC.Constraint,[171 133]);
FEC.VitDec.TraceBackDepth = 5*(FEC.Constraint-1); % approx for constraint
FEC.VitDec.OpMode = 'trunc';
FEC.VitDec.DecType = 'hard';

%% Start Sim
for frame = 1:OFDM.NumFrames
    [OFDM,RF,FEC] = ofdm_transmit(DATA,OFDM,RF,SIM,FEC,frame);
    [OFDM,CHANNEL] = channel_apply(CHANNEL,OFDM,RF,SIM);
    [OFDM,DATA] = ofdm_receive(DATA,OFDM,RF,CHANNEL,SIM,FEC,frame);
end
DATA = image_bitdecode(DATA,OFDM);
error_percentage = 100*sum(DATA.BitData == DATA.BitDataHat)/length(DATA.BitData);
ber = sum(DATA.BitData ~= DATA.BitDataHat)/length(DATA.BitData);
fprintf('Correct bit recovered = %.6f%% \n',error_percentage);
fprintf('BER = %.6f \n',ber);
figure,
subplot(121)
imshow(DATA.Data), title('Original Image')
subplot(122)
imshow(DATA.DataHat), title('Received Image')