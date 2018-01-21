function ANI = ANIindex(final_rr,fs) 
% fs=16;
% ----------------------------------- Filtering -------------------------

% PLANO A: Wavelets
% n=5;
% %required frequency levels
% levels=[5 4 3];
% new_sig=wavelets(final_rr,n,fs,levels);%M
% USAR SO DETAIL COEFS DO NIVEL 5 E 4 tambem para testar
% APPCOEFS É PARTE LOW FREQ. NAO INTERESSA EM WAVELETS
% Wavelet Band pass filtering between 0.15Hz and 0.4_0.5Hz (variações parasimpáticas muito influenciadas pelo ciclo respiratorio)
% To do so, we'll use a daubechie wavelet with 4 coefficients and we'll reconstruct the levels 3 to 5


% %% PLANO B: lomb scargle 
% x=1:length(data_filtered);
% [Pr,fr,A,z0,A0,ofac]=lomb_scargle(x,data_filtered,[0.15 0.5]);
% %%
% [Pxx,F]=lomb([x' data_filtered']);
% %%
% figure('Name', 'Lomb-scargle periodogram');
% plot(fr{1,1},Pr{1,1},'r-',F,Pxx,'b-');
% %%
% figure;
% Fs=200;
% L=144;
% Y = fft(data_filtered);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
% plot(f,P1)
% title('Single-Sided Amplitude Spectrum of S(t)')
% xlabel('f (Hz)')
% ylabel('|PPG(f)|')

%% PLANO C
ordem=4;
wc= 0.5;
fc= wc/(0.5*fs); 
[b,a]=butter(ordem,fc); 
e1= filter(b, a, final_rr);

%Passo 2
wc= 0.15;
fc= wc/(0.5*fs); 
[b,a]=butter(ordem,fc,'high'); 
new_sig= filter(b, a, e1);
% filtfilt
%% normalization in 0 to 0.2 normalized units
norm_sig=new_sig*0.1/max(abs(new_sig));
%% UPDATE: TRY TO CATCH MAX-MIN-MAX!!!!!!!!!!!!!!
% Detection of local maxima and minima
rr_wave=norm_sig; % .^2; highlights more the differences
[max_val,max_pos]=findpeaks(rr_wave,'MinPeakHeight',0);
[min_val,min_pos]=findpeaks(-rr_wave,'MinPeakHeight',0);
%%
% Interpolation among maximum and minimum points
max_interp = interp1([1 max_pos length(rr_wave)],[0 max_val 0],1:length(rr_wave),'linear');
min_interp = interp1([1 min_pos length(rr_wave)],[0 rr_wave(min_pos) 0],1:length(rr_wave),'linear');
%% area calculation (ou apenas sum). Division into four subareas (A1,A2,A3,A4), measured between the lower and upper envelope

A1_max=sum(max_interp(1:fs:16*fs),'omitnan');
A1_min=abs(sum(min_interp(1:fs:16*fs),'omitnan'));
A1=A1_max+A1_min;

A2_max=sum(max_interp(16*fs+1:fs:32*fs),'omitnan');
A2_min=abs(sum(min_interp(16*fs+1:fs:32*fs),'omitnan'));
A2=A2_max+A2_min;

A3_max=sum(max_interp(32*fs+1:fs:48*fs),'omitnan');
A3_min=abs(sum(min_interp(32*fs+1:fs:48*fs),'omitnan'));
A3=A3_max+A3_min;

A4_max=sum(max_interp(48*fs+1:fs:end),'omitnan');
A4_min=abs(sum(min_interp(48*fs+1:fs:end),'omitnan'));
A4=A4_max+A4_min;

A=[A1,A2,A3,A4];

% AUCmin (minimum of the 4 subareas)
AUCmin=min(A);

a=5.1; b=1.2; % determined on a 200 patients dataset in order to keep coherence between
% respiratory influence on RR series and ANI

% Calculation of ANI
ANI=100*(a*AUCmin+b)/(12.8);
% 0: high ANS response to stress; (low anti-nociception).is often but not always related to an inadequate analgesia/nociception balance.
% 100: low ANS response to stress
figure('Name','rr_wave');
plot((1:length(rr_wave))/fs,rr_wave,'r-',(1:length(max_interp))/fs,max_interp,'b-',(1:length(min_interp))/fs,min_interp,'b-');
str = sprintf('ANI is %f', ANI);
set(gca,'fontsize',13)
title(str)
xlabel('Time (seconds)')
ylabel('normalized units')
legend('RR series','interpolation');
