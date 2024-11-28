function [y,FEC] = data_deinterleave(x,FEC)
% deinterleave from random permutation elements
y = deintrlv(x,FEC.IntrlvElem);
end