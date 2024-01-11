clear all ;close all; clc;

%% SAMPLING

fs = 80*10^3;
f1 = 2.5*10^3; f2 = 6*10^3;
dt = 1/fs;
d = 0.008;
n = 0:dt:d;

x1 = cos(2*pi*f1*n);
x2 = sin(2*pi*f2*n);
x = x1 + x2;

xd4 = [];
xd12 = [];
nd4 = [];
nd12 = [];

i = 1;
j = 1;

for i=1:4:length(n);
    xd4(j) = x(i);
    nd4(j) = n(i);
    j = j+1;
end

i=1;
j=1;

for i=1:12:length(n);
    xd12(j) = x(i);
    nd12(j) = n(i);
    j = j+1;
end

figure(1)
subplot(311);
stem(n,x);title("original signal");xlabel("discrete time points n (s)");ylabel("amplitude");legend("x[n]");

subplot(312);
stem(nd4,xd4);title("factor 4 downsampled signal");xlabel("discrete time points n (s)");ylabel("amplitude");legend("xd4[n]");

subplot(313);
stem(nd12,xd12);title("factor 12 downsampled signal");xlabel("discrete time points n (s)");ylabel("amplitude");legend("xd12[n]");

N = length(n);
Nd4 = length(nd4);
Nd12 = length(nd12);

w = linspace(-fs/2,fs/2,N);
wd4 = linspace((-fs/2)/4,(fs/2)/4,N);
wd12 = linspace((-fs/2)/4,(fs/2)/4,N);

XD = fftshift(fft(x,N));
XD4 = fftshift(fft(xd4,N));
XD12 = fftshift(fft(xd12,N));

figure(2)
subplot(311)
plot(w,abs(XD)/N);title("mag. freq. response of original signal");xlabel("frequency (f)");ylabel("amplitude");legend("XD[f]");

subplot(312)
plot(wd4,abs(XD4)/N);title("mag. freq. response of downsampled signal (4)");xlabel("frequency (f)");ylabel("amplitude");legend("XD4[f]");

subplot(313)
plot(wd12,abs(XD12)/N);title("mag. freq. response of downsampled signal (12)");xlabel("frequency (f)");ylabel("amplitude");legend("XD12[f]");

xd4lin = interp1(nd4,xd4,n,'linear');
xd4cub = interp1(nd4,xd4,n,'PCHIP');

xd12lin = interp1(nd12,xd12,n,'linear');
xd12cub = interp1(nd12,xd12,n,'PCHIP');

figure(3)
subplot(211)
plot(n,x,"green");
hold on;
plot(n,xd4lin,"red");
hold on;
plot(n,xd4cub,"blue");title("original, linear and cubic interpolated signal (4)");xlabel("time (s)");ylabel("amplitude");legend("x(t)","x4lin_ip(t)","x4cub_ip(t)");
hold off;

subplot(212)
plot(n,x,"green");
hold on;
plot(n,xd12lin,"red");
hold on;
plot(n,xd12cub,"blue");title("original, linear and cubic interpolated signal (12)");xlabel("time (s)");ylabel("amplitude");legend("x(t)","x12lin_ip(t)","x12cub_ip(t)");
hold off;


%% QUANTIZATION

fsp = 80*10^3; f1p = 2000; f2p = 400;
t = 0:dt:d;

xt = 3*cos(2*pi*f1p*t)+sin(2*pi*f2p*t);
nq1 = 4; nq2 = 6; b = max(xt); a = min(xt);


xq4 = floor(((xt-a)/(b-a))*((2^nq1) - 1)) * (((b-a)/((2^nq1) - 1))) + a;

xq12 = floor(((xt-a)/(b-a))*((2^nq2) - 1)) * (((b-a)/((2^nq2) - 1))) + a;


%xq4 = (xt-a)/(xt-b)*((2^nq1) - 1) * (((b-a)/((2^nq1) - 1)) + a);
%xq12 = (xt-a)/(xt-b)*((2^nq2) - 1) * (((b-a)/((2^nq2) - 1)) + a);

figure(4)
subplot(211)
plot(t,xt);
hold on;
plot(t,xq4);title("original signal and 4 bit-quantized signal");xlabel("time (s)");ylabel("amplitude");legend("x(t)","x4q(t)");
hold off;

subplot(212)
plot(t,xt);
hold on;
plot(t,xq12);title("original signal and 6 bit-quantized signal");xlabel("time (s)");ylabel("amplitude");legend("x(t)","x6q(t)");

% I had a amplitude miscalculation on quantization and ı couldnt find the
% paramaeter ı mistook so ı continued in calculatin sqnr.

e4 = xt - xq4;
e12 = xt - xq12;

xvar = var(xt);
xq4var = var(xq4);
xq12var = var(xq12);
e4var = var(e4);
e12var = var(e12);

sqnr4 = 10*log10((xvar^2) / (e4var^2))

sqnr6 = 10*log10((xvar^2) / (e12var^2))



