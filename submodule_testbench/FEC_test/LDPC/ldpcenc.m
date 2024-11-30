function y = ldpcenc(x,LDPC)
%LDPCENCODE Encode the input bit vector x using the LDPC code specified in
% struct 'LDPC'.
% Based largely on Vodafone-Chair-5g-nr-ldpc at github.com

% Check if input has correct size
assert(size(x,2)==LDPC.numInfBits, 'Error: The given input vector has length %d, correct input length is %d!', length(x), LDPC.numInfBits)

% Encode 
y = mod(x * LDPC.G, 2); % Note: This is very inefficient.

% tic 
% Gsp = sparse(logical(LDPC.G));
% xsp = sparse(logical(x));
% yy = mod(xsp*Gsp,2);
% toc
% ----> Worse performance

% Puncturing of 2*Z first systematic bits
y = y(:,2*LDPC.Z+1:end);