%读入原始信号 s(n) 干净语音
%[y,Fs]=audioread('p232_003.wav');
[y,Fs]=audioread('clean.wav');
% sound(y,Fs);
n=length(y);
%绘制波形
figure(1);
plot(y);
xlabel('时间');
ylabel('幅度');
title('原始信号的波形');
grid on;

%产生观测信号x(n)=s(n)+v(n) 
%v(n)噪声，选择白噪声
%v = wgn(n,1,0.2);
% sound(v,Fs);
% 产生观测信号
%x=awgn(y,20);
[x,Fs]= audioread('5dB_noisy.wav');
%sound(x,Fs);
audiowrite('snp232_003.wav',x,Fs);
figure(2);
plot(x);
xlabel('时间');
ylabel('幅度');
title('观测信号的波形');
grid on;

%维纳滤波的频域非因果实现
Rxx=xcorr(x);
Gxx=fft(Rxx,n);
Rxy=xcorr(x,y);
Gxs=fft(Rxy,n);
H=Gxs./Gxx;
Ps=fftn(y);
S=H.*Ps;
ss=real(ifft(S));
ss=ss(1:n);
sound(ss,Fs);

%audiowrite('weina.wav',ss,Fs);
figure(3);
plot(ss);
xlabel('时间');
ylabel('幅度');
title('恢复原始信号的波形');
grid on;

%结果分析
figure(4);
t=(40000:44800);
plot(t,ss(40000:44800,1 ),'-',t,y(40000:44800,1),'-.');
ylabel('幅度');
xlabel('时间');
legend('恢复的原始信号信号','原始信号');
title('信号比较');
grid on;


