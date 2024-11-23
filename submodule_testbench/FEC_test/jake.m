function Impulse_response = jake(v,N)
% Jakes Rayleigh fading
% Reference: Khaled Ramadan, UDEMY, The Complete MATLAB Course for Wireless
% Communication Systems

c=3*10^8; % light speed
f=2*10^9;
Tc=24.4*10^-6;
v=v/3.6; % velocity
wM=2*pi*f*v/c;
theta=2*pi*rand(1,200);
t=0:0.6*Tc:Tc;
No=16;

A=hadamard(No);

for j=1:N
    T=zeros(1,length(t));
    for n=1:No
    
    An=2*pi*(n-0.5)/(4*No);
    wn=wM*cos(An);
    Bn=pi*n/No;
    T=A(j,n)*sqrt(2/No)*(cos(Bn)+1i*sin(Bn))*cos(wn*t+theta(n))+T;
    end
    T1(j,:)=T;
end
    
for z1=1:N
T1(z1,:)=(T1(z1,:)-mean(T1(z1,:)))/std(T1(z1,:));
end
PDP_dB=[0;-1;-9;-10;-15;-20];
PDP=10.^(PDP_dB/10);
tp=sqrt(PDP)/norm(sqrt(PDP));
Impulse_response=tp.*T1(:,1);






    

