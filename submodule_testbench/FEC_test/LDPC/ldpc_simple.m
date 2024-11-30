% Steps to LDPC process
%{
1. Create parity check matrix H (sparse matrix)
2. Make sure the last n-k columns of H is invertible in GF(2)
3. Get generator matrix G
4. Do a modulo-2 matrix multiplication input*G: codeword c = mod2(input*G)
5. Ensure that the contraint H*c_T = 0 is met
6. Pass codeword c to the rest of the Transmit processing chain
7. At receiver, after some receiver processing chain, produces the received codeword c_hat
8. Decode using Bit-Flipping or Belief Propagation (Sum-Product) algorithm
9. At same decoder, check if c_hat*H_T = 0 (all zero vector) is met
10. At same decoder, If condition of number 9 is met, no bits are in error, output the decoded message bit
11. At same decoder, if condition of number 9 is not met, check the index that are not 0 and do correction to the received codeword to meet condition 9.
    After iteratively correct to produce condition 9, output the final
    decoded message bit.
%}

% Define the parity-check matrix (H) in systematic form
H = [1 0 1 1;  % Parity-check matrix in systematic form
     0 1 1 0];

% Extract I and P from H (H = [I | P])
I = H(:, 1:2); % Identity matrix (first 2 columns)
P = H(:, 3:4); % Parity submatrix (last 2 columns)

% Define the input bits (infoBits)
infoBits = [1; 0]; % Column vector

% Compute the parity bits
parityBits = mod(P * infoBits, 2); % Perform binary matrix multiplication mod 2

% Form the codeword (concatenate infoBits and parityBits)
codeword = [infoBits; parityBits];

% Display the results
disp('Parity-check matrix (H):');
disp(H);

disp('Input bits (infoBits):');
disp(infoBits');

disp('Parity bits (computed):');
disp(parityBits');

disp('Encoded codeword:');
disp(codeword');


% MATLAB function ldpcEncode()
cfg = ldpcEncoderConfig();