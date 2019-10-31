clc
clear
close all
disp('recording...');
recObj = audiorecorder;
recordblocking(recObj,5);
disp('recorded');
%%
disp('playing recorded sound...');
play(recObj);
pause(7);
%%
Fs = 8000;
fs = 8000;
y = getaudiodata(recObj);
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
Fs = 8000;
hlpf = fdesign.lowpass('Fp,Fst,Ap,Ast',3.0e3,3.5e3,0.5,50,Fs);
D = design(hlpf);
yd = filter(D,y);
disp('playing denoised sound');
sound(yd,Fs);
figure,plot(yd);
title('denoise');