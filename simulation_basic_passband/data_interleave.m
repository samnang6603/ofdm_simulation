function [y,FEC] = data_interleave(x,FEC)
FEC.IntrlvElem = randperm(length(x));
y = intrlv(x,FEC.IntrlvElem);
end