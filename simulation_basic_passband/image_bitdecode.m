function DATA = image_bitdecode(DATA)
% remove padded 0 from frame handling non-integer number of frames from
% DATA.BitData and BitDataHat
DATA.BitData(end-DATA.FrameNumZeroPadded+1:end) = [];
DATA.BitDataHat(end-DATA.FrameNumZeroPadded+1:end) = [];
% decode grayscale
data_hat_reshape = reshape(DATA.BitDataHat,8,DATA.TotalBits/8).';
data_hat_decode = mat2str(data_hat_reshape);
data_hat_decode([1 end]) = []; % remove brackets caused by mat2str()
str_len = 16; % num of char including bit and delimiter
for m = 0:DATA.NumElement-1
    si = str_len*m + 1;
    ei = str_len*m + str_len-1;
    DATA.DataHat(m+1) = uint8(bin2dec(data_hat_decode(si:ei)));
end
DATA.DataHat = reshape(DATA.DataHat,DATA.Size);
end