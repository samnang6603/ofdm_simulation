function [DATA,OFDM] = data_validate_frame_fit(DATA,OFDM)
% Handle non-integer frame size. Pad data with zero for the fractional i-th
% frame to make data_frame input valid length
if rem(OFDM.NumFrames,1) == 0 % if frame is integer do nothing
    return
end

% Find remainder of data
remainder_data = DATA.TotalBits - floor(OFDM.NumFrames)*OFDM.NumBitsPerFrame;

% To make one full frame, we pad zeros to DATA.BitData make remainder data 
% to one full frame = OFDM.NumBitsPerFrame
DATA.FrameNumZeroPadded = OFDM.NumBitsPerFrame-remainder_data;
DATA.BitData = [DATA.BitData
                zeros(DATA.FrameNumZeroPadded,1,'like',...
                      DATA.BitData)];

% Also update BitDataHat length
DATA.BitDataHat = zeros(length(DATA.BitData),1,'like',DATA.BitData);

% Ceil up the number of frame to next largest integer
OFDM.NumFrames = ceil(OFDM.NumFrames);

end