clear all; close all; clc;
%% DM_MOD
d = 0.05; fs =2*10^3; f1 = 60; f2 = 90; A1 = 0.8; A2 = 0.7; dt = 1/fs; t = 0:dt:d;
mt = A1*cos(2*pi*f1*t) + A2*sin(2*pi*f2*t);
dmin = (A1*2*pi*f1/fs) + (A2*2*pi*f2/fs);

size_mt = length(mt);
mq = zeros(0,size_mt);
e = zeros(0,size_mt);
eq = zeros(0,size_mt);
encoded = zeros(0,size_mt);
mq(1) = 0;
e(1) = mt(1);

encoded(1) = 0;
for i=1:size_mt-1;
    e(i) = mt(i+1) - mq(i);
    eq(i) = dmin*sign(e(i));
    mq(i+1) = mq(i) + eq(i);

    if eq(i) > 0;
        encoded(i) = 1;
    else 
        encoded(i) = 0;
    end
       
end

size_en = length(encoded);
mqp = zeros(1,size_en);
mqp(1) = 0;

for i=1:size_en;
    if encoded(i) == 0;
        mqp(i+1) = mqp(i) - dmin;
    else
        mqp(i+1) = mqp(i) + dmin;
    end
end

cf = 0.1;
n = 100;
b = fir1(n,cf);
out = conv2(mqp,b,'same');

%{
figure(2)
plot(t,mt);
hold on;
plot(t,out);
hold off;
%}

figure(1)
subplot(211)
plot(t,mt);
hold on;
plot(t,mqp);xlabel("time(s)");ylabel("Amplitude");title("Message and decoder output signal w.r.t. dmin = 0.3487");legend("m(t)","mqp(t)");
hold off;

subplot(212)
plot(t,mt);
hold on;
plot(t,out);xlabel("time(s)");ylabel("Amplitude");title("Message and LPF output signal w.r.t. dmin = 0.3487");legend("m(t)","out(t)");
hold off;



mq(1) = 0;
e(1) = mt(1);
dmin2 = 0.1;
for i=1:size_mt-1;
    e(i) = mt(i+1) - mq(i);
    eq(i) = dmin2*sign(e(i));
    mq(i+1) = mq(i) + eq(i);

    if eq(i)>0
        encoded(i) = 1;
    else 
        encoded(i) = 0;
    end
end

size_en = length(encoded);
mqp = zeros(1,size_en);
mqp(1) = 0;

for i=1:size_en;
    if encoded(i) == 0;
        mqp(i+1) = mqp(i) - dmin2;
    else
        mqp(i+1) = mqp(i) + dmin2;
    end
end

cf = 0.1;
n = 100;
b = fir1(n,cf);
out = conv2(mqp,b,'same');

figure(2)
subplot(211)
plot(t,mt);
hold on;
plot(t,mqp);xlabel("time(s)");ylabel("Amplitude");title("Message and decoder output signal w.r.t. dmin = 0.1");legend("m(t)","mqp(t)");
hold off;

subplot(212)
plot(t,mt);
hold on;
plot(t,out);xlabel("time(s)");ylabel("Amplitude");title("Message and LPF output signal w.r.t. dmin = 0.1");legend("m(t)","out(t)");
hold off;

mq(1) = 0;
e(1) = mt(1);
dmin3 = 0.5;
for i=1:size_mt-1;
    e(i) = mt(i+1) - mq(i);
    eq(i) = dmin3*sign(e(i));
    mq(i+1) = mq(i) + eq(i);

    if eq(i)>0
        encoded(i) = 1;
    else 
        encoded(i) = 0;
    end
end

size_en = length(encoded);
mqp = zeros(1,size_en);
mqp(1) = 0;

for i=1:size_en;
    if encoded(i) == 0;
        mqp(i+1) = mqp(i) - dmin3;
    else
        mqp(i+1) = mqp(i) + dmin3;
    end
end

cf = 0.1;
n = 100;
b = fir1(n,cf);
out = conv2(mqp,b,'same');

figure(3)
subplot(211)
plot(t,mt);
hold on;
plot(t,mqp);xlabel("time(s)");ylabel("Amplitude");title("Message and decoder output signal w.r.t. dmin = 0.5");legend("m(t)","mqp(t)");
hold off;

subplot(212)
plot(t,mt);
hold on;
plot(t,out);xlabel("time(s)");ylabel("Amplitude");title("Message and LPF output signal w.r.t. dmin = 0.5");legend("m(t)","out(t)");
hold off;












