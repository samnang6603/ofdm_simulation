function [cw,FEC] = data_fec_encode(x,FEC)
switch FEC.Type
    case 'Convolutional'
        cw = convenc(x,FEC.Trellis);
    case 'LDPC'
        % Data segmentation to make LDPC input valid
        [x,FEC] = fec_encode_ldpc_segment(x,FEC); 
        nseg = FEC.NumSegment; % LDPC number of segment
        ninfo = FEC.Ecfg.NumInformationBits; % LDPC input bit lenght
        cwinfo = FEC.Ecfg.BlockLength; % LDPC output codeword length
        cw = zeros(FEC.SegCodeWordLength,1,'like',x); % allocate codeword
        for k = 1:nseg
            tmp_idx_x = ninfo*(k-1)+1:ninfo*(k-1) + ninfo;
            tmp_idx_cw = cwinfo*(k-1)+1:cwinfo*(k-1) + cwinfo;
            cw(tmp_idx_cw) = ldpcEncode(x(tmp_idx_x),FEC.Ecfg);
        end
end
end

function [xseg,FEC] = fec_encode_ldpc_segment(x,FEC)
data_len = length(x);
ldpc_input_len = FEC.Ecfg.NumInformationBits;
if data_len == ldpc_input_len
    xseg = x;
    FEC.NumSegment = 1;
    return
end
% get remainder length
remainder_len = rem(data_len, ldpc_input_len); 

% find the length that is an integer multiple of fit length
fit_len = data_len - remainder_len; 

% Find number of full segments (segment that don't need padding)
num_full_segments = fit_len/ldpc_input_len;

% Total number of segment is number of full segments + 1 (segment that
% needs padding)
num_segments = num_full_segments + 1;

% Find number of zero pad for the segment that needs padding
num_zero_pad = ldpc_input_len - remainder_len;

% concatenate the zeros
xseg = [x
        zeros(num_zero_pad,1,'like',x)];

FEC.NumSegment = num_segments;

% Final length of LDPC input including all segment
FEC.SegInputLength = num_segments*ldpc_input_len; 

% Final length of codeword output after padding
FEC.SegCodeWordLength = num_segments*FEC.Ecfg.BlockLength;

% Save the number of zeros being padded so receiver can remove
FEC.NumZeroPadded = num_zero_pad;

end