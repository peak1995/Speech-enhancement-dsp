clc;
clear all;
 % 选择干净原始音频文件
[x fs]=audioread('clean.wav'); 
N=length(x);
x = x(1:N,1);     % 如果是双声道，取单通道
max_x = max(x);

% 添加含噪语音
[y fs]=audioread('5dB_noisy.wav');
%sound(y,fs);
noise_estimated = y(1:1000*fs/1000,1);  %取前1秒做为噪声进行去噪

fft_y = fft(y);
fft_n = fft(noise_estimated);
E_noise = sum(abs(fft_n)) /length(noise_estimated);
mag_y = abs(fft_y);
phase_y = angle(fft_y);   % 保留相位信息
mag_s = mag_y - E_noise;
mag_s(mag_s<0)=0;
 
% 恢复
fft_s = mag_s .* exp(1i.*phase_y);
s = ifft(fft_s);
%sound(s,fs);
audiowrite('pujian.wav',s,fs);
subplot(321);plot(x);ylim([-1.5,1.5]);title('原始干净信号');xlabel('时间');ylabel('幅度');
subplot(323);plot(y);ylim([-1.5,1.5]);title('含噪信号');xlabel('时间');ylabel('幅度');
subplot(325);plot(real(s));ylim([-1.5,1.5]);title('谱减法去噪后信号');xlabel('时间');ylabel('幅度');
subplot(322);myspectrogram(x,fs);colormap(jet);time=(0:length(x)-1)/fs;axis([0 max(time*1000) 0 8000]);title('干净信号的语谱图');xlabel('时间');ylabel('频率');
subplot(324);myspectrogram(y,fs);colormap(jet);time=(0:length(x)-1)/fs;axis([0 max(time*1000) 0 8000]);title('含噪信号的语谱图');xlabel('时间');ylabel('频率');
subplot(326);myspectrogram(s,fs);colormap(jet);time=(0:length(x)-1)/fs;axis([0 max(time*1000) 0 8000]);title('谱减法增强后信号的语谱图');xlabel('时间');ylabel('频率');
