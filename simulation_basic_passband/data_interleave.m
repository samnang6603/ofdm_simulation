function [y,FEC] = data_interleave(x,FEC)
% interleave from random permutation elements
FEC.IntrlvElem = randperm(length(x));
y = intrlv(x,FEC.IntrlvElem);
end