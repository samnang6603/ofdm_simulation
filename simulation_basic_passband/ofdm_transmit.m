function [OFDM,RF,FEC] = ofdm_transmit(DATA,OFDM,RF,SIM,FEC,frame)
% Data per frame
OFDM.data_frame = DATA.BitData((frame-1)*OFDM.NumBitsPerFrame+1:...
                  (frame-1)*OFDM.NumBitsPerFrame+OFDM.NumBitsPerFrame);

if SIM.FEC.InterleaveToggle
    [OFDM.data_frame,FEC] = data_interleave(OFDM.data_frame,FEC);
end
OFDM.data_frame_codeword = OFDM.data_frame;

if SIM.FEC.Toggle
    OFDM.data_frame_codeword = data_fec_encode(OFDM.data_frame_codeword,FEC);
end

[OFDM.data_frame_codeword,OFDM] = ofdm_transmit_adjust_mapper_input(...
                        OFDM.data_frame_codeword,OFDM);

OFDM.TxSymbols = qammod(OFDM.data_frame_codeword,OFDM.M,InputType='bit');

if SIM.CHANNEL.PilotSignalToggle
    % Create Pilot signal
    yp = zeros(OFDM.NumCarrierPilotPerFrame,1);
    xp = 2*randi([0 1],OFDM.NumPilotPerFrame,1)-1;
    spc = OFDM.NumPilotSpacing;
    pilot_loc = 1:spc+1:OFDM.NumCarrierPilotPerFrame;
    yp(pilot_loc) = xp;

    % Insert Pilot signal
    for k = 0:OFDM.NumPilotPerFrame-1
        yp(2+k*(spc+1):1+k*(spc+1)+spc) = OFDM.TxSymbols(k*spc+1:k*spc+spc);
    end

    % transmitter with pilot
    ifft_sig_pilot = ifft(yp,OFDM.NumCarrierPilotPerFrame);

    % cylic prefix
    cyclic_idx = OFDM.NumCarrierPilotPerFrame-...
        OFDM.NumCyclicPilotSymsPerFrame+1:OFDM.NumCarrierPilotPerFrame;
    cext_data = [ifft_sig_pilot(cyclic_idx); ifft_sig_pilot];

    OFDM.PilotSignal = xp;
    OFDM.PilotSignalIFFT = ifft_sig_pilot(pilot_loc);
    OFDM.PilotSignalLocation = pilot_loc;
else
    % transmitter no pilot
    ifft_sig = ifft(OFDM.TxSymbols,OFDM.NumCarriersPerFrame);

    % cyclic prefix
    cyclic_idx = OFDM.NumCarriersPerFrame-...
        OFDM.NumCyclicSymsPerFrame+1:OFDM.NumCarriersPerFrame;
    cext_data = [ifft_sig(cyclic_idx); ifft_sig];
end
OFDM.TxAir = cext_data;

% If IQimbalance impairment is included
if SIM.RF.IMPAIRMENT.IQImbalanceToggle
    OFDM.TxAir = rf_upconvert_impairment_IQimb(OFDM.TxAir,RF);
end

if SIM.RF.PassBandProcessingToggle
    % upconvert to passband carrrier
    [OFDM,RF] = rf_upconvert(OFDM,RF);
end

% simulate doppler shift if toggle
if SIM.RF.IMPAIRMENT.DopplerToggle
    OFDM.RFTxAir = rf_upconvert_impairment_doppler(OFDM.RFTxAir, RF);
end

end

function [y,OFDM] = ofdm_transmit_adjust_mapper_input(x,OFDM)
% validate and adjust the length of output y to see if valid for mapper
% input.
% Then Find next value of the length of x that is divisible by n then pad zero
% by num_zero_pad to that value to output y
%
l = length(x);
n = OFDM.BitsPerSymbol;
mapper_len_cond = rem(l,n) ~= 0;
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