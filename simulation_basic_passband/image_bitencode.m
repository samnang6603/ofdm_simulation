function DATA = image_bitencode(DATA)
% encode grayscale image to binary representation
x = DATA.Data;
DATA.Size = size(x);
DATA.NumElement = numel(x);
cbin = num2cell(dec2bin(x(:)));
DATA.TotalBits = 8*DATA.NumElement;
DATA.BitData = zeros(DATA.TotalBits,1,'uint8');
for m = 0:DATA.NumElement-1
    str = sprintf(' %s',cbin{m+1,:});
    bit_num = sscanf(str,'%d');
    DATA.BitData(8*m+1:8*m+8) =  bit_num;
end
DATA.BitDataHat = zeros(DATA.TotalBits,1,'uint8'); % placeholder for decoded data
DATA.DataHat = zeros(DATA.NumElement,1,'uint8');
end