%读入原始信号 s(n) 干净语音
[y,Fs]=audioread('clean.wav');
% sound(y,Fs);
n=length(y);
%产生观测信号x(n)=s(n)+v(n) 
%v(n)噪声，选择白噪声
%v = wgn(n,1,0.2);
% sound(v,Fs);
% 产生观测信号
%x=awgn(y,20);
[x,Fs]= audioread('5dB_noisy.wav');
%sound(x,Fs);
% audiowrite('snp232_003.wav',x,Fs);
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
%sound(ss,Fs);
%audiowrite('weina.wav',ss,Fs);
% %结果分析
% 
% figure(4);
% t=(40000:44800);
% plot(t,ss(40000:44800,1 ),'-',t,y(40000:44800,1),'-.');
% ylabel('幅度');
% xlabel('时间');
% legend('恢复的原始信号信号','原始信号');
% title('信号比较');
% grid on;
audiowrite('weinafaend.wav',ss,fs);
subplot(321);plot(y);xlabel('时间');ylabel('幅度');title('原始信号的波形');
subplot(323);plot(x);xlabel('时间');ylabel('幅度');title('观测信号的波形');
subplot(325);plot(ss);xlabel('时间');ylabel('幅度');title('维纳滤波恢复原始信号的波形');
subplot(322);myspectrogram(y,fs);colormap(jet);time=(0:length(x)-1)/fs;axis([0 max(time*1000) 0 8000]);title('干净信号的语谱图');xlabel('时间');ylabel('频率');
subplot(324);myspectrogram(x,fs);colormap(jet);time=(0:length(x)-1)/fs;axis([0 max(time*1000) 0 8000]);title('含噪信号的语谱图');xlabel('时间');ylabel('频率');
subplot(326);myspectrogram(ss,fs);colormap(jet);time=(0:length(x)-1)/fs;axis([0 max(time*1000) 0 8000]);title('维纳滤波增强后信号的语谱图');xlabel('时间');ylabel('频率');

