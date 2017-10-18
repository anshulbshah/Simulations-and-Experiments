fs=44.1*10^3;
Tp = 1;
f1=20;
f2=20*10^3;
A=1;

ti = linspace(0,1,fs);
x_ss = A*sin(2*pi*f1*ti+pi*(f2-f1)./Tp*ti.^2);
x_log_ss = A*sin((2*pi*f1*Tp/log(f2/f1))*((f2/f1).^ti - 1));
subplot(2,1,1)
plot(ti,x_log_ss);
%spectrogram(x_log_ss,'yaxis');
%spectrogram(x_log_ss,,'yaxis');
subplot(2,1,2)
spectrogram(x_log_ss,256,250,256,fs,'yaxis')


figure(2)
x_rev = fliplr(x_log_ss);
alpha = (0.3*log(10)*log2(f2/f1))/Tp;
inv_filter = x_rev.*exp(-alpha*ti);
subplot(2,1,1);
plot(ti,inv_filter);
subplot(2,1,2);
spectrogram(inv_filter,256,250,256,fs,'yaxis');

figure(3)
plot([-1.*fliplr(ti) ti(2:end)],conv(x_log_ss,inv_filter));

figure(4)
b = fir2(50,[0 0.3 0.5 1],[1 1 0 0]);
plot(b)
y=filter(b,1,x_log_ss)
figure(5)
plot(y)

figure(6)
check = conv(y,inv_filter);
plot(check);


figure(7)
yn = (x_log_ss+0.5.*(x_log_ss.^3))./1.5;
non_lin = conv(yn,inv_filter);
plot(non_lin);
figure(8)
spectrogram(non_lin,256,250,256,fs,'yaxis');
