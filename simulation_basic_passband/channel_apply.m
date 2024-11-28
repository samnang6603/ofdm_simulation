function [OFDM,CHANNEL] = channel_apply(CHANNEL,OFDM,RF,SIM)
% Apply channel based on fading and AWGN conditions

% Handle fading (Jake model or custom fading model)
if SIM.CHANNEL.Fading
    CHANNEL.ImpulseResponse = jake(120,CHANNEL.TapLength);
    % CHANNEL.ImpulseResponse = 1/sqrt(2)*(randn(CHANNEL.TapLength,1) +...
    %                           1j*randn(CHANNEL.TapLength,1)); 
else
    % No fading, simple channel
    CHANNEL.ImpulseResponse = 1;
end

if SIM.RF.PassBandProcessingToggle
    signal_in = OFDM.RFTxAir;
else
    signal_in = OFDM.TxAir;
end

% Apply the channel impulse response to the passband signal (upconverted)
signal_out = filter(CHANNEL.ImpulseResponse,1,signal_in);

% Apply delay if applicable
signal_out = CHANNEL.DelayObj(signal_out);

if SIM.CHANNEL.AWGN
    if SIM.RF.PassBandProcessingToggle
        % Adding AWGN to the passband signal at the defined SNR
        OFDM.RFTxAirChannel = awgn(signal_out,SIM.SNR,'measured');
    else
        OFDM.TxAirChannel = awgn(signal_out,SIM.SNR,'measured');
    end
else
    if SIM.RF.PassBandProcessingToggle
        % Adding AWGN to the passband signal at the defined SNR
        OFDM.RFTxAirChannel = signal_out;
    else
        OFDM.TxAirChannel = signal_out;
    end
end
end

