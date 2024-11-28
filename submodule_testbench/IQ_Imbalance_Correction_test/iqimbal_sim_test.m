N = 4096;
sps = 16;
x = randi([0 sps-1],N,1);
txSig = qammod(x,sps,UnitAveragePower=true);
%txSig = pskmod(x,sps,pi/4);
scatterplot(txSig,sps), title('TX Signal Constellation')

ampImb_dB = 5;
phImb_deg = 15;

ampImb = db2mag(ampImb_dB);
phImb = deg2rad(phImb_deg);

Igain = 1;%10.^(0.5*ampImb_dB/20).*exp(-0.5i*phImb);
Qgain = 10.^(-0.5*ampImb_dB/20).*exp(0.5i*phImb);

rxSig = real(txSig).*Igain + 1i*imag(txSig).*Qgain;

scatterplot(rxSig,sps), title('RX Signal Constellation')

compSig = complex(zeros(length(rxSig),1),zeros(length(rxSig),1));
mu = 1e-5;

w = complex(0,0);
for n = 1:1000
    for m = 1:length(rxSig)
        compSig(m) = rxSig(m) + w*conj(rxSig(m));
        w(:) = w - mu*compSig(m)*compSig(m);
    end
end

% 4-tap w and 4-points average
% w = complex(zeros(4,1),zeros(4,1));
% for n = 1:10
%     for m = 1:4:length(rxSig)
%         compSig(m:m+3) = rxSig(m:m+3) + filter(w,1,conj(rxSig(m:m+3)));
%         w(:) = w - mu*compSig(m:m+3).*compSig(m:m+3);
%     end
% end

%scatterplot(compSig,sps), title('RX Comp Signal Constellation')

scatterplot(compSig(N - 1000:end),sps), title('Last 1000 samples')