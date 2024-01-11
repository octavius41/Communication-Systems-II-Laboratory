close all; clear all; clc;

fs = 3000; d = 0.4; N = 4; Tb = 0.1; dt = 1/fs;
t = dt:dt:0.4; tp = dt:dt:0.1; tp0 = tp(1:75); tp1 = tp(1:100);

s0t = [zeros(1,length(tp0)) ones(1,length(tp0)) ones(1,length(tp0)) zeros(1,length(tp0))];
s1t = [ones(1,length(tp1))*-1 zeros(1,length(tp1)) ones(1,length(tp1))];

figure(1)
subplot(211)
plot(tp,s0t);xlabel("time(s)");ylabel("amplitude");title("0 bit corresponding signal");legend("s0(t)");
subplot(212)
plot(tp,s1t);xlabel("time(s)");ylabel("amplitude");title("1 bit corresponding signal");legend("s1(t)");

b = [1 0 1 0];

st = [s1t s0t s1t s0t];

figure(2)
plot(t,st);xlabel("time(s)");ylabel("amplitude");title("transmitted signal s(t) for bit sequence b = 1 0 1 0");legend("s(t)");

sp = sum(abs(st).^2)/length(st);
snr1 = 15;
snr2 = 0;

snrlin1 = 10^(0.1*snr1);
snrlin2 = 10^(0.1*snr2);

varn1 = sp / snrlin1;
varn2 = sp / snrlin2;

n1t = sqrt(varn1).*randn(1,length(st));
n2t = sqrt(varn2).*randn(1,length(st));

r1t = st + n1t;
r2t = st + n2t;

wb = Tb/dt;
k = 1:1:N;

figure(3)

plot(t,st);
hold on;
plot(t,r1t);xlabel("time(s)");ylabel("amplitude");title("transmitted and received signals for dB = 15");legend("s1(t)","r1(t)");
hold off;

figure(4)
plot(t,st);
hold on;
plot(t,r2t);xlabel("time(s)");ylabel("amplitude");title("transmitted and received signals for dB = 0");legend("s2(t)","r2(t)");
hold off;

r0k_1 = zeros(1,N);
r1k_1 = zeros(1,N);

for k = 1:1:N
   n = (k-1)*wb+1:k*wb;
   r0k_1(k) = sum(r1t(n).*s0t((n-(k-1)*wb)));
   r1k_1(k) = sum(r1t(n).*s1t((n-(k-1)*wb)));
end

r0k_2 = zeros(1,N);
r1k_2 = zeros(1,N);

for k = 1:1:N
   n = (k-1)*wb+1:k*wb;
   r0k_2(k) = sum(r2t(n).*s0t((n-(k-1)*wb)));
   r1k_2(k) = sum(r2t(n).*s1t((n-(k-1)*wb)));
end

figure(5)
scatter([1:N],r0k_1);
hold on;
scatter([1:N],r1k_1);xlabel("number of symbols vector N=4");ylabel("0 - 1 bit chance");title("corellator output for dB = 15");legend("r0_1(k)","r1_1(k)");
hold off;

figure(6)
scatter([1:N],r0k_2);
hold on;
scatter([1:N],r1k_2);xlabel("number of symbols vector N=4");ylabel("0 - 1 bit chance");title("corellator output for dB = 0");legend("r0_2(k)","r1_2(k)");
hold off;
