% function [Output] = kalman(NoisyInput,Fs,Noise)
% 卡尔曼滤波函数
[Input, Fs] = audioread('clean.wav');
Input = Input(:,1);
% Noise = normrnd(0,sqrt(0.001),size(Input));
[NoisyInput Fs] = audioread('5dB_noisy.wav');
Time = (0:1/Fs:(length(Input)-1)/Fs)';
Noise=NoisyInput(1:1000*Fs/1000,1);
% 基础变量设置
Len_windowT = 0.0025; % 窗长2.5ms
Hop_Percent = 1; % 窗移占比
AR_Order = 20; % 自回归滤波器阶段
Num_Iter = 7; %迭代次数

%分帧加窗处理
Len_WinFrame = fix(Len_windowT * Fs);
Window = ones(Len_WinFrame,1);
[Frame_Signal, Num_Frame] = KFrame(NoisyInput, Len_WinFrame, Window, Hop_Percent);

%初始化
H = [zeros(1,AR_Order-1),1];   % 观测矩阵
R = var(Noise);     % 噪声方差
[FiltCoeff, Q] = lpc(Frame_Signal, AR_Order);   % LPC预测，得到滤波器的系数
P = R * eye(AR_Order,AR_Order);   % 误差协方差矩阵
Output = zeros(1,size(NoisyInput,1));   % 输出信号
Output(1:AR_Order) = NoisyInput(1:AR_Order,1)';   %初始化输出信号
OutputP = NoisyInput(1:AR_Order,1);

% 迭代器的次数.
i = AR_Order+1;
j = AR_Order+1;

%进行卡尔曼滤波
for k = 1:Num_Frame   %一次处理一帧信号
    jStart = j;     %跟踪每次迭代AR_Order+1的值.
    OutputOld = OutputP;    %为每次迭代保留第一批AROrder预估量
    
    for l = 1:Num_Iter
        A = [zeros(AR_Order-1,1) eye(AR_Order-1); fliplr(-FiltCoeff(k,2:end))];
        
        for ii = i:Len_WinFrame
            %Kalman滤波基本方程式
            OutputC = A * OutputP;
            Pc = (A * P * A') + (H' * Q(k) * H);
            K = (Pc * H')/((H * Pc * H') + R);
            OutputP = OutputC + (K * (Frame_Signal(ii,k) - (H*OutputC)));
            Output(j-AR_Order+1:j) = OutputP';
            P = (eye(AR_Order) - K * H) * Pc;
            j = j+1;
        end       
        i = 1;
        if l < Num_Iter
            j = jStart;
            OutputP = OutputOld;
        end     
        % 更新滤波后信号的lpc
        [FiltCoeff(k,:), Q(k)] = lpc(Output((k-1)*Len_WinFrame+1:k*Len_WinFrame),AR_Order);
    end
end
Output = Output';
%绘制
MAX_Am(1)=max(Input);
MAX_Am(2)=max(NoisyInput);
MAX_Am(3)=max(Output);
audiowrite('kaerman.wav',Output,Fs);
subplot(321);plot(Input);ylim([-max(MAX_Am),max(MAX_Am)]);xlabel('时间');ylabel('幅度');title('原始信号的波形');
subplot(323);plot(NoisyInput);ylim([-max(MAX_Am),max(MAX_Am)]);xlabel('时间');ylabel('幅度');title('含噪信号的波形');
subplot(325);plot(Output);ylim([-max(MAX_Am),max(MAX_Am)]);xlabel('时间');ylabel('波形');title('Kalman滤波信号的波形');
subplot(322);myspectrogram(Input,Fs);colormap(jet);time=(0:length(x)-1)/Fs;axis([0 max(time*1000) 0 8000]);title('干净信号的语谱图');xlabel('时间');ylabel('频率');
subplot(324);myspectrogram(NoisyInput,Fs);colormap(jet);time=(0:length(x)-1)/Fs;axis([0 max(time*1000) 0 8000]);title('含噪信号的语谱图');xlabel('时间');ylabel('频率');
subplot(326);myspectrogram(Output,Fs);colormap(jet);time=(0:length(x)-1)/Fs;axis([0 max(time*1000) 0 8000]);title('卡尔曼滤波增强后信号的语谱图');xlabel('时间');ylabel('频率');
