function y = ldpcdec_bp(llr, LDPC)
% Decode LDPC using BP(Belief-Propagation)-&SP(Sum-Product) Algorithm
% llr : Log-Likelihood Ratio (LLR)
% alg_type : BA-SPA, or BF 
% Author: Sam A. 11/30/2024
% Based largely on Vodafone-Chair-5g-nr-ldpc at github.com
%{ 
References
[1] E. Sharon, S. Litsyn and J. Goldberger, "An efficient message-passing 
schedule for LDPC decoding," 2004 23rd IEEE Convention of Electrical and 
Electronics Engineers in Israel, Tel-Aviv, Israel, 2004, pp. 223-226.  

[2] J. G. Proakis, M. Saleh, "Fundamentals of Wireless Communication
Systems," 2nd Edition, Pearson Education, Prentice Hall, 2014, pp. 746-747
%}

%{
BP-SPA Algorithm Description:
Notation:
- v and c as variable and check node index, respectively
- datahat is the detected codeword
- Pv as the LLR of the received codeword after constellation demapping:
    Pv = ln(pv(0)/pv(1)) where:
    pv(0) is the probability that the true bit is 0
    pv(1) is the probability that the true bit is 1
- Qvc as message passing from variable node v to check node c
- Rcv as message passing from check node c to variable node v
- Qv as the total computed bits at variable nodes v
- SUM(x) as the capital Sigma notation for Sum of x elements
- PROD(x) as the capital Pi notation for Product x elements
- N(v)\c as the set of variable nodes v connected to check nodes c
- N(c)\v as the set of check nodes c connected to variable nodes v
- phi(x) as the log domain of tanh(x) defined as:
    phi(x) = (sign(x), -ln(tanh(|x|/2)) 
    and its inverse:
    iphi(sign(x),x) = (-1)^sign * -ln(tanh(x/2)), 
    where sign(x) is the indicator function 1_{x<0} [1]

The algorithm relies on the messge-passing via the Tanner Graph between the
variable nodes (number of column of parity check matrix) and the check
nodes (number of row of parity check matrix). The sparse matrix's ones are
the non-zeros that are the edges, that is, the nodes with connections.

1. The BP-SPA algorithm begins by passing the message llr from variable
nodes v to neighboring check nodes c:
    
    Qvc = Pv +  SUM(Rcv) 

    Since, initially, that is, at time t = 0, Rcv = 0, this reduces to:

    Qvc = Pv
             /pv(0)\
        = ln|-------|
             \pv(1)/

2. The check nodes received the message and use the below relation to send
message back to its neighboring variable nodes v:

    Rcv = 2*atanh(PROD(tanh(Qvc/2))) such that v is the subset of N(c)\v

3. Using the log domain of tanh(x), where the PROD becomes SUM as in the
 log arguments, the response message from c to v can be expressed as:

    Rcv = iphi(SUM(phi(Qvc))) such that v is the subset of N(c)\v

4. The variable nodes received the response message and send back another
message to the neighboring check nodes to start another iteration but this
time the Qvc is updated with:

    Qvc = Pv + SUM(Rcv) such that c is the subset of N(v)\c

5. The process repeats and the bits detection is based on computing the
total likelihoods at the variable nodes:

    Qv_total = Qv + SUM(Rcv) such athat c is the subser of N(v)\c

Then a detected bit datahat can be computed with hard-decision threshold

    datahat = 1 if Qv_total < 0; 0 if else

6. The process repeats until the constraint condition or convergence metric
 is met:
    H*datahat_T = 0 where
    - H is the parity check matrix
    - datahat is the decoded bits
    - _T is the matrix transpose operation
 or until it meets the iteration limit preset.

7. Additional convergence metric can be measured by monitoring changes in
Qvc. If minimal suggest algorithm has stabilized or is stabilizing to a
convergence point.

8. As additional consideration for numerical stability, the ln() will be 
implemented using tol = 2.2251e-308, which is the spacing between
64-bits floating point 1.0 to the next larger 64-bits floating point, as
the smallest possible argument.

This is the basic idea of BP-SPA algorithm implementation. This algorithm 
is the original vanilla algorithm that is based on Gallager's work called
the flooding-scheduling. Flooding-scheduling is when in each iteration,
messages from all variable nodes are sent to check nodes and vice versa.
This is by no means optimal in terms of efficiency in computation and
resources.

However this code will implement the more efficient variant like in [1]. It
is called dual-scheduling. Instead of passing all messages from variable
nodes to check nodes and vice versa, the check nodes are gone over in some
order and for each node, messages are sent into and out from:

    1. Compute Qvc for each v subset of N(c): send all Qvc messages into
    check node c.
    2. Compute Rcv for each v subset of N(c): send all Rcv messages out
    from check node c.

%}

% Initialization
% Pad the punctured llr with zeros
llr = [zeros(1,2*LDPC.Z), llr];

% Find number of non-zeros entries in H (the edges)
numEntries = nnz(LDPC.H);

% Possible smallest value
tol = realmin;

% The received LLR (using Qv notation as in [1])
Qv = llr;

% Allocate sparse message passing matrix
Rcv = spalloc(LDPC.numParBits, LDPC.numTotBits + 2*LDPC.Z, numEntries);

% Decode
for iter = 1:LDPC.iterations % loop over max iteration
    for checkIdx = 1:LDPC.numParBits % loop over all check nodes

        % Find all neighboring variable nodes of current check node
        nbVarNodes = find(LDPC.H(checkIdx,:)==1);
        
        % Tmp update llr for the current check node 
        tmpllr = Qv(nbVarNodes) - full(Rcv(checkIdx,nbVarNodes));

        % Compute SUM(phi(Qvc)) with 
        % phi(x) = (sign(x), -ln(tanh(|x|/2)) 
        % Let Ssign = sign(x), Smag = SUM(-ln(tanh(|x|/2)), x = tmpllr
        Smag = sum(-log(tol+tanh(abs(tmpllr)/2)));

        % Sign depends on the number of negative elements
        if mod(sum(tmpllr<0),2) == 0
            Ssign = +1;
        else
            Ssign = -1;
        end

        % Loop all neighboring variable nodes
        for varIter = 1:length(nbVarNodes)

            varIdx = nbVarNodes(varIter);
            Qtmp = Qv(varIdx) - Rcv(checkIdx, varIdx);
            QtmpMag = -log(tol+tanh(abs(Qtmp)/2));
            % Note: +minVal in order to deal with llr=0;
            % implementation can be improved
            QtmpSign = sign(Qtmp+tol);

            % Update message passing matrix
            % From reference: Rcv = phi^-1(S-phi(Qtmp))
            %Rcv(checkIdx, varIdx) = Ssign*QtmpSign * (-log(minVal+tanh(abs(Smag-QtmpMag)/2)));
            % To avoid sparse matrix indexing
            tmpRcv = Ssign*QtmpSign * (-log(tol+tanh(abs(Smag-QtmpMag)/2)));

            % Update Qv. From reference: Qv = Qtmp + Rcv
            %Qv(varIdx)  = Qtmp + Rcv(checkIdx, varIdx);
            Qv(varIdx)  = Qtmp + tmpRcv;
        end

    end

end
% Convert Qv to decoded bits
y = zeros(1,LDPC.numInfBits);
y(Qv(1:LDPC.numInfBits)<0) = 1;

    