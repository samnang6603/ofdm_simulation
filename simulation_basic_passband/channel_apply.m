function [OFDM,CHANNEL] = channel_apply(CHANNEL,OFDM,SIM)
if SIM.Fading
    CHANNEL.ImpulseResponse = jake(120,CHANNEL.TapLength);
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