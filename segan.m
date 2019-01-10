[y,fs]=audioread('clean.wav');
[x,fs]=audioread('5dB_noisy.wav');
[s,fs]=audioread('enhanced_5dB_noisy.wav');
subplot(321);plot(y);xlabel('时间');ylabel('幅度');title('原始信号的波形');
subplot(323);plot(x);xlabel('时间');ylabel('幅度');title('观测信号的波形');
subplot(325);plot(s);xlabel('时间');ylabel('幅度');title('SEGAN增强信号的波形');
subplot(322);myspectrogram(y,fs);colormap(jet);time=(0:length(x)-1)/fs;axis([0 max(time*1000) 0 8000]);title('干净信号的语谱图');xlabel('时间');ylabel('频率');
subplot(324);myspectrogram(x,fs);colormap(jet);time=(0:length(x)-1)/fs;axis([0 max(time*1000) 0 8000]);title('含噪信号的语谱图');xlabel('时间');ylabel('频率');
subplot(326);myspectrogram(s,fs);colormap(jet);time=(0:length(x)-1)/fs;axis([0 max(time*1000) 0 8000]);title('SEGAN增强后信号的语谱图');xlabel('时间');ylabel('频率');
