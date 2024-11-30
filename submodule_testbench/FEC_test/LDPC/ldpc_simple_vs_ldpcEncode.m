% Define the parity-check matrix H for a (7, 4) LDPC code
H = [1 0 0 1 1 0 0
0 1 0 1 0 1 0
0 0 1 0 1 1 1
];

% Extract P from H (H is of the form [P | I])
P = H(:, 1:4);  % The parity part of the matrix (first 3 columns)
I = eye(4);     % Identity matrix (4x4)

% Transpose P to get P^T 
P_T = P';

% Create the generator matrix G (systematic form)
G = [I P_T];  % Concatenate I and P^T to get a 4x7 generator matrix

% Define the input message (4 bits)
inputMessage = [1 0 1 0]; % Example message bits

% Step 1: Manual Encoding using Matrix Multiplication (mod 2)
% Encode the message by multiplying with the generator matrix (mod 2)
encodedMessageManual = mod(inputMessage * G, 2);

% Step 2: MATLAB's built-in LDPC encoding
% Define the LDPC code (this is a simple example code with rate 1/2)

cfg = ldpcEncoderConfig(sparse(logical(H)));

% Use MATLAB's ldpcEncode function to encode the input message
encodedMessageMATLAB = ldpcEncode(inputMessage', cfg);

% Display the results
disp('Manual Encoding Output:');
disp(encodedMessageManual); % Output from manual encoding

disp('MATLAB ldpcEncode() Output:');
disp(encodedMessageMATLAB); % Output from MATLAB's ldpcEncode

% Compare the results
if isequal(encodedMessageManual, encodedMessageMATLAB)
    disp('The manual encoding and MATLAB ldpcEncode() output match!');
else
    disp('The manual encoding and MATLAB ldpcEncode() output do not match.');
end