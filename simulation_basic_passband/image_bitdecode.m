function DATA = image_bitdecode(DATA,OFDM)
% decode grayscale
data_hat_reshape = reshape(DATA.BitDataHat,8,DATA.TotalBits/8).';
data_hat_decode = mat2str(data_hat_reshape);
data_hat_decode([1 end]) = [];
for m = 0:DATA.NumElement-1
    si = OFDM.M*m + 1;
    ei = OFDM.M*m + OFDM.M-1;
    DATA.DataHat(m+1) = uint8(bin2dec(data_hat_decode(si:ei)));
end
DATA.DataHat = reshape(DATA.DataHat,DATA.Size);
end