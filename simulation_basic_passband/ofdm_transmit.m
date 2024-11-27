function [OFDM,RF,FEC] = ofdm_transmit(DATA,OFDM,RF,SIM,FEC,frame)
% Data per frame
OFDM.data_frame = DATA.BitData((frame-1)*OFDM.NumBitsPerFrame+1:...
                  (frame-1)*OFDM.NumBitsPerFrame+OFDM.NumBitsPerFrame);

if SIM.Interleave
    [OFDM.data_frame,FEC] = data_interleave(OFDM.data_frame,FEC);
end
OFDM.data_frame_codeword = OFDM.data_frame;

if SIM.FECToggle
    OFDM = data_fec_encode(OFDM,FEC);
end

OFDM.TxSymbols = qammod(OFDM.data_frame_codeword,OFDM.M,InputType='bit',UnitAveragePower=1);

if SIM.PilotSignalToggle
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
if RF.IMPAIRMENT.IQImbalance.Toggle
    OFDM.TxAir = rf_upconvert_impairment_IQimb(OFDM.TxAir,RF);
end

% if SIM.ChannelEstimation
%     % Add pilot
%     yp = zeros(OFDM.NumCarrierPilotPerFrame,1);
%     Xp = 2*randi([0 1],OFDM.NumPilotPerFrame,1)-1;
%     spc = OFDM.NumPilotSpacing;
%     pilot_loc = 1:spc+1:OFDM.NumCarrierPilotPerFrame;
%     yp(pilot_loc) = Xp;
%     for k = 0:OFDM.NumPilotPerFrame-1
%         yp(2+k*(spc+1):1+k*(spc+1)+spc) = OFDM.TxSymbols(k*spc+1:k*spc+spc);
%     end
% 
%     % transmitter with pilot
%     ifft_sig_pilot = ifft(yp,OFDM.NumCarrierPilotPerFrame);
% 
%     % cylic prefix
%     cyclic_idx = OFDM.NumCarrierPilotPerFrame-...
%         OFDM.NumCyclicPilotSymsPerFrame+1:OFDM.NumCarrierPilotPerFrame;
%     cext_data = [ifft_sig_pilot(cyclic_idx); ifft_sig_pilot];
% 
%     OFDM.PilotSignal = Xp;
%     OFDM.PilotSignalLocation = pilot_loc;
% else
%     % transmitter no pilot
%     ifft_sig = ifft(OFDM.TxSymbols,OFDM.NumCarriersPerFrame);
% 
%     % cyclic prefix
%     cyclic_idx = OFDM.NumCarriersPerFrame-...
%         OFDM.NumCyclicSymsPerFrame+1:OFDM.NumCarriersPerFrame;
%     cext_data = [ifft_sig(cyclic_idx); ifft_sig];
% end
% OFDM.TxAir = cext_data;

if RF.PassBandProcessingToggle
    % upconvert to passband carrrier
    [OFDM,RF] = rf_upconvert(OFDM,RF);
end

end