close all; clear all;
% ----------------------------- development ------------------------------
%%
%  [siginfo,Fs]=wfdbdesc('mimicdb/224/');
% intensive care unit patient
db=1;
if db==1
    %------------- mimic
fs=125;
[tm, signal]=rdsamp('mimicdb/218/',[],50000);
ecg=signal(:,2); % lead II

elseif db==2
    %------------ philips
fs=250;
ecg=measure(1).data(2.03001*10^5:7.03*10^5)';
end
%%
% Pan tompkins
rr=Pan_Tompkins_rt(ecg,fs);
%%
% Normalization (mean centered and resampling at 8Hz: valor arbitrario)
rr_centered=rr-mean(rr);
n_fs=8; %nova frequencia de sampling
rrc_resampled=resample(rr_centered,n_fs,1); %experimentar com 4

%% artificial test
% n_fs=16; %nova frequencia de sampling
% t=0:530;
% rr_centered=sin(2*pi*n_fs*500*t);
% rrc_resampled=resample(rr_centered,n_fs,1); %experimentar com 4
%% NORMALIZATION
final_rr=[];
sr_rr=0;
contador=0;
init=1;
ANI_values=[];

% Sem sobreposiçao da janela (com sobreposiçao, ciclo for extende-se a todo
% o sinal).The PhysioDoloris monitor continuously displays an average measurement of ANI recorded during the previous 60 seconds.
%(64 IN THIS CASE). 4S MOVING PERIOD É SUFICIENTE PARA UMA BOA ANALISE

for j=1:4*n_fs:length(rrc_resampled)-64*n_fs
    sr_rr=rrc_resampled(j:j+64*n_fs).^2;
    S=sqrt(sum(sr_rr));
    final_rr=rrc_resampled(j:j+64*n_fs)./S;
    ANI=ANIindex(final_rr,n_fs);
%     fprintf ('ANI= %4.2f ;\n',ANI);
    ANI_values=cat(2,ANI_values,ANI);
    %ANI=100 bom
    %ANI=0 stress
    pause
end
%%
temp=(1:4*n_fs:length(rrc_resampled)-64*n_fs)/n_fs;
figure('Name','RR series and respective ANI index')

% subplot(312)
plot((1:length(rrc_resampled))/n_fs,rrc_resampled*100,'r-',temp,ANI_values,'b-');
xlabel('Time (seconds)')
ylabel('RR (seconds)')
legend('RR series (resampled)', 'ANI')
% 
% subplot(311)
% plot((1:length(rr_centered))/n_fs,rr,'r-');
% ylabel('Time (seconds)')
% xlabel('Heart Beats')
% legend('RR series')
% 
% subplot(313)
% plot(temp,ANI_values,'b-');
% xlabel('Time (seconds)')
% ylabel('ANI')
% legend('ANI values')

