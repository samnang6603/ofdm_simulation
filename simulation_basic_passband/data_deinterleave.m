function [y,FEC] = data_deinterleave(x,FEC)
y = deintrlv(x,FEC.IntrlvElem);
end