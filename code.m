clc
clear
close all
% disp('recording...');
% recObj = audiorecorder;
% recordblocking(recObj,5);
% disp('recorded');
% %%
% disp('playing recorded sound...');
% play(recObj);
% pause(7);
%%
% y = getaudiodata(recObj);
input = 'audio.wav';
[in,fs] = audioread(input);
[y,fs] = audioread(input);
y = y(:, 1);
%info = audioinfo(input);
figure,plot(y);
title('input');
y = awgn(y,40);

figure,plot(y);
title('awgn');
disp('playing added noise...');
sound(y);
pause(7)
%%
% yd = wdenoise(y,3,'Wavelet','db3','DenoisingMethod','UniversalThreshold','ThresholdRule','Soft','NoiseEstimate','LevelDependent');
disp('playing denoised sound...');
[thr,sorh,keepapp]=ddencmp( 'den' , 'wv' ,y);  
yd=wdencmp( 'gbl' ,y, 'db3' ,2,thr,sorh,keepapp);  
sound(yd);
figure,plot(yd);
title('denoise');
%% Best till now
%'Fp,Fst,Ap,Ast' (passband frequency, stopband frequency, passband ripple, stopband attenuation)
hlpf = fdesign.lowpass('Fp,Fst,Ap,Ast',3.0e3,3.5e3,0.5,50,fs);
D = design(hlpf);
x = filter(D,y);
disp('playing denoised sound');
sound(x,fs);
figure,plot(x);
title('denoise');

%% freq shaper
g = 45;
freqP = [1000, 1500, 2550, 5000];
T = 1/fs;
len = length(x);
p = log2(len);
p = ceil(p);
N = 2^p;
X = fft(x,N);
gain = linspace(0,0,N);
gain = gain.';
%for first stage
firstC = (0.3*(g-1))/freqP(1);
k = 0;
while (k<=freqP(1)*N/fs)
   gain(k+1) = firstC*k/(N*T) + 1; %paper step 1
    gain(N-k) = gain(k+1);
    k=k+1;
end
% for second stage
secondC = firstC*freqP(1)+1;
secondC2 = (freqP(2)-freqP(1))/5;
while (k<= freqP(2)*N/fs)
    gain(k+1) = ((secondC-1)*exp(-((k/(N*T))-freqP(1))/secondC2)) + 1;
    gain(N-k) = gain(k+1);
    k=k+1;
end

thirdC = ((secondC-1)*exp(-freqP(2)/secondC2)) +1;
thirdC2 = (freqP(3)-freqP(2))/5;
while (k<= freqP(3)*N/fs)
    gain(k+1) = ((thirdC-g)*exp(-((k/(N*T)-freqP(2)))/thirdC2)) +g;
    gain(N-k) = gain(k+1);
    k=k+1;
end

while (k<=freqP(4)*N/fs)
    gain(k+1) = g;
    gain(N-k) = g;
    k=k+1;
end

fifthC = g;
fifthC2 = (fs/2-freqP(4))/5;
while (k<=N/2)
    gain(k+1) = ((fifthC-1)*exp(-((k/(N*T))-freqP(4))/fifthC2)) +1;
    gain(N-k) = gain(k+1);
    k=k+1;
end
k_v = (0:N-1)/N;
figure,plot(k_v,gain);
title('Gain');

k_v = k_v*fs;
k_v = k_v(1:N/2+1);
plot(k_v,gain(1:N/2+1));
title('Frequency Shaper Transfer Function');
xlabel('Frequency (Hertz)');
ylabel('Gain');
xlim([0 10000]);
Y = X+gain; % for X refer line no.27
y = real(ifft(Y,N));
y = y(1:len);
t=(0:1/fs:(len-1)/fs);
figure;
subplot(2,1,1);
plot(t,y,'r');
title('Signal after addition of gain');
subplot(2,1,2);
plot(t,x);
title('Adjusted Signal');

%%
sound(y,fs);

%% amplitude shaper
out1=fft(out);
phse=angle(out1);
mag=abs(out1)/N;
threshold=100
for i=1:size(mag)
    if(mag((i)>threshold)
        mag(i)=threshold
    end
    out(i)=mag(i)*exp(j*phse(i));
end
outfinal=ifft(out);
