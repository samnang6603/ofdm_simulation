function CHANNEL = chanest_dft_enhance(CHANNEL,nfft)
h_est = ifft(CHANNEL.EstFreqResponse);
h_est = h_est(1:CHANNEL.TapLength);
CHANNEL.EstFreqResponse = fft(h_est,nfft);
end