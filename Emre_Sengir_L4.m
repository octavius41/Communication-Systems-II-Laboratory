clear all; close all; clc;

fs = 27*10^3; Rb = 18*10^3; N = 1000; dum = 1:1:1000; dt = 1/fs; k = 2;

s = randi(4,1,N);

for i = 1:1:1000
    if s(i) == 4
        s(i) = 3;
    elseif s(i) == 3
        s(i) = 1;
    elseif s(i) == 2
        s(i) = -1;
    else
        s(i) = -3;
    end
end

Rs = Rb/k;  %% symbol rate 
symps = fs/Rs;  %% number of samples per symbol
roloff = [0 0.5 1];
trunc = 10;

rcos1 = rcosdesign(roloff(1),trunc,symps,"sqrt");
rcos2 = rcosdesign(roloff(2),trunc,symps,"sqrt");
rcos3 = rcosdesign(roloff(3),trunc,symps,"sqrt");

Rcosf1 = fftshift(fft(rcos1,N));
Rcosf2 = fftshift(fft(rcos2,N));
Rcosf3 = fftshift(fft(rcos3,N));

w = linspace(-fs/2,fs/2,N);

figure(1)
subplot(231)
plot(rcos1);title("Raised cosine with B = 0 in time domain");xlabel("time(s)");ylabel("amplitude");legend("rcos1(t)");

subplot(232)
plot(rcos2);title("Raised cosine with B = 0.5 in time domain");xlabel("time(s)");ylabel("amplitude");legend("rcos2(t)");


subplot(233)
plot(rcos3);title("Raised cosine with B = 1 in time domain");xlabel("time(s)");ylabel("amplitude");legend("rcos3(t)");


subplot(234)
plot(w,abs(Rcosf1)/N);title("Raised cosine with B = 0 in frequency domain");xlabel("frequency(f)");ylabel("amplitude");legend("Rcos1(f)");


subplot(235)
plot(w,abs(Rcosf2)/N);title("Raised cosine with B = 0.5 in frequency domain");xlabel("frequency(f)");ylabel("amplitude");legend("Rcos2(f)");

subplot(236)
plot(w,abs(Rcosf3)/N);title("Raised cosine with B = 1 in frequency domain");xlabel("frequency(f)");ylabel("amplitude");legend("Rcos3(f)");

snr = 20;
%%
tran1 = upfirdn(s,rcos1,symps);

r1n = awgn(tran1,snr);

y1_out = upfirdn(r1n,rcos1,1,symps);
y1_out = y1_out(trunc:(end-(trunc+1)));

eyediagram(y1_out,symps,dt,0);  title("eyediagram of output with B = 0, trunc = 10 samples on each plot"); xlabel("time(s)"); ylabel("amplitude");
%%

tran2 = upfirdn(s,rcos2,symps);

r2n = awgn(tran2,snr);

y2_out = upfirdn(r2n,rcos2,1,symps);
y2_out = y2_out(trunc:(end-(trunc+1)));

eyediagram(y2_out,symps,dt,0); title("eyediagram of output with B = 0.5, trunc = 10 samples on each plot"); xlabel("time(s)"); ylabel("amplitude");
%%

tran3 = upfirdn(s,rcos3,symps);

r3n = awgn(tran3,snr);

y3_out = upfirdn(r3n,rcos3,1,symps);
y3_out = y3_out(trunc:(end-(trunc+1)));

eyediagram(y3_out,symps,dt,0);  title("eyediagram of output with B = 1, trunc = 10 samples on each plot"); xlabel("time(s)"); ylabel("amplitude");
%%

rolexc = [0.25 0.33 0.50 0.67 0.75 1];

bw = 6000;

symrate = zeros(1,6);

for i = 1:1:6
    symrate(i) = bw/(1+rolexc(i));
end

figure(5)
plot(rolexc,symrate); title("excess bandwidth with rolloff factors and symbol rate relation"); xlabel("rolloff factor (excess bandwidth percentage)"); ylabel("symbol rate"); legend("ratio");


q = UnitQuaternion(rotx(0.2));




