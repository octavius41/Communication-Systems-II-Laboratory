close all;

fs = 10^4; dt = 1/fs; fc = 2500; bitd = 0.02;bit_s = bitd/dt; A = 5;
b = [0 1 0 1 1 0 1 0 1 0 1 1 1 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 1 1 1 0 1 0 1 0 0 0 0 0 0 0 1 1 1 1 0 0 0 1];

b_oi = [];
for i=1:bit_s
    b_oi = [b_oi;b];
end
b_oi = b_oi(:)';

d = length(b)*bitd;
t = 0:dt:(d-dt);
tp = 0:dt:(bitd-dt);

bask_md = [];

for i=1:length(b)
    if b(i) == 0
        bask_md = [bask_md 0*tp];
    else
        bask_md = [bask_md A*cos(2*pi*fc*(((i-1)*bitd)+tp))];
    end
end

bfsk_md = [];
f0 = 500;
f1 = 750;
for i=1:length(b)
    if b(i) == 0
        bfsk_md = [bfsk_md A*cos(2*pi*f0*(((i-1)*bitd)+tp))];
    else
        bfsk_md = [bfsk_md A*cos(2*pi*f1*(((i-1)*bitd)+tp))];
    end
end

na = awgn(bask_md,10);
nf = awgn(bfsk_md,10);

randoms = randn(1,length(bask_md))/log(10);
spa = sum(abs(bask_md).^2)/length(bask_md);
spf = sum(abs(bfsk_md).^2)/length(bfsk_md);

snr = -10;

snrlin = 10^(0.1*snr);

varna = spa / snrlin;
varnf = spf / snrlin;

nta = sqrt(varna).*randn(1,length(bask_md));
ntf = sqrt(varnf).*randn(1,length(bfsk_md));

bask_n = bask_md + nta;
bfsk_n = bfsk_md + ntf;

s_ask_0 = 0*tp;
s_ask_1 = A*cos(2*pi*fc*(tp));

s_fsk_0 = A*cos(2*pi*f0*(tp));
s_fsk_1 = A*cos(2*pi*f1*(tp));

bask_dmd = [];
bfsk_dmd = [];

for i = 1:length(b)
    L1 = xcorr(bask_n((1:bit_s)+bit_s*(i-1)),A*cos(2*pi*fc*((bitd*(i-1)+tp))),0);
    L0 = xcorr(bask_n((1:bit_s)+bit_s*(i-1)),0*tp,0);
    th = L1-L0
    if th>20
    bask_dmd(i) = 1; 
    else
    bask_dmd(i) = 0;
    end
end

for i = 1:length(b)
    L1 = xcorr(bfsk_n((1:bit_s)+bit_s*(i-1)),s_fsk_1,0);
    L0 = xcorr(bfsk_n((1:bit_s)+bit_s*(i-1)),s_fsk_0,0);
    th = L1-L0;
    if th>0
    bfsk_dmd(i) = 1; 
    else
    bfsk_dmd(i) = 0;
    end
end

bask_dmd_p = [];
for i=1:bit_s
    bask_dmd_p = [bask_dmd_p;bask_dmd];
end
bask_dmd_p = bask_dmd_p(:)';

bfsk_dmd_p = [];
for i=1:bit_s
    bfsk_dmd_p = [bfsk_dmd_p;bfsk_dmd];
end
bfsk_dmd_p = bfsk_dmd_p(:)';

figure(1)
subplot(411)
plot(t,b_oi);title("message signal");xlabel("time(s)");ylabel("amplitude");legend("m(t)");

subplot(412)
plot(t,bask_md);title("bask modulated signal");xlabel("time(s)");ylabel("amplitude");legend("bask(t)");

subplot(413)
plot(t,bask_n);title("noisy bask signal");xlabel("time(s)");ylabel("amplitude");legend("bask_in_noise(t)");ylim([-5 5]);

subplot(414)
plot(t,bask_dmd_p);title("demodulated bask signal");xlabel("time(s)");ylabel("amplitude");legend("demod(t)");

figure(2)
subplot(411)
plot(t,b_oi);title("message signal");xlabel("time(s)");ylabel("amplitude");legend("m(t)");

subplot(412)
plot(t,bfsk_md);title("bfsk modulated signal");xlabel("time(s)");ylabel("amplitude");legend("bfsk(t)");

subplot(413)
plot(t,bfsk_n);title("noisy bfsk signal");xlabel("time(s)");ylabel("amplitude");legend("bfsk_in_noise(t)");ylim([-5 5]);

subplot(414)
plot(t,bfsk_dmd_p);title("demodulated bfsk signal");xlabel("time(s)");ylabel("amplitude");legend("demod(t)");





