% clear all;
% clc;
% INTRODUZIR OUTROS PACIENTES TODOS!!!

% % ------------------------------- SPI index ------------------------------

%% Virtual test
% n_fs=8; %new sampling frequency
% ppga_temp=[interp1([0 7],[0 10],0:7/(n_fs*60*7):7) interp1([7 7.5],[10 3],7:0.5/(n_fs*60*0.5):7.5) interp1([7.5 8],[3 10],7.5:0.5/(n_fs*60*0.5):8) interp1([8 8.5],[10 0],8:0.5/(n_fs*60*0.5):8.5) interp1([8.5 8.6],[0 5],8.5:0.1/(n_fs*60*0.1):8.6) interp1([8.6 9.3],[5 7],8.6:0.7/(n_fs*60*0.7):9.3) interp1([9.3 9.4],[7 2],9.3:0.1/(n_fs*60*0.1):9.4) interp1([9.4 9.5],[2 10],9.4:0.1/(n_fs*60*0.1):9.5) interp1([9.5 11],[10 15],9.5:1.5/(n_fs*60*1.5):11) interp1([11 12],[15 15],11:1/(n_fs*60*1):12) interp1([12 14],[15 10],12:2/(n_fs*60*2):14) randi([10 13],1,(n_fs*60*5)) interp1([19 20],[13 2],19:1/(n_fs*60*1):20) randi([2 15],1,(n_fs*60*9)) interp1([29 30],[10 15],29:1/(n_fs*60*1):30) interp1([30 30.5],[15 2],30:0.5/(n_fs*60*0.5):30.5) interp1([30.5 31],[2 14],30.5:0.5/(n_fs*60*0.5):31)];
% rr_t=[interp1([0 7],[85 70],0:7/(n_fs*60*7):7) randi([70 77],1,(n_fs*60*2)) interp1([9 9.1],[70 110],9:0.1/(n_fs*60*0.1):9.1) interp1([9.1 12],[110 80],9.1:2.9/(n_fs*60*2.9):12) interp1([12 12.1],[80 60],12:0.1/(n_fs*60*0.1):12.1) interp1([12 13],[60 100],12:1/(n_fs*60*1):13) interp1([13 15],[100 83],13:2/(n_fs*60*2):15) interp1([15 15.1],[83 75],15:0.1/(n_fs*60*0.1):15.1) interp1([15.1 15.2],[75 80],15.1:0.1/(n_fs*60*0.1):15.2) randi([79 81],1,(n_fs*60*3.8)) interp1([19 20],[79 110],19:1/(n_fs*60*1):20) randi([105 110],1,(n_fs*60*3)) interp1([23 23.1],[100 80],23:0.1/(n_fs*60*0.1):23.1) randi([80 90],1,(n_fs*60*3.9)) interp1([27 28],[80 70],27:1/(n_fs*60*1):28) interp1([28 31],[70 67],28:3/(n_fs*60*3):31)];
% rr_temp=diff(rr_t);
%
% figure('Name', 'Artificial signal');
% subplot 211
% tempo=0:31/(length(ppga_temp)):31;
% plot(tempo(2:end),ppga_temp./100);
% legend('PPGA')
% subplot 212
% tempo=0:31/(length(rr_temp)-1):31;
% plot(tempo,rr_temp);
% legend('HR')

%% teste de qualidade do sinal!
[siginfo,Fs]=wfdbdesc('mimicdb/408/');
%%
[tm, signal]=rdsamp('mimicdb/401/',[],50000); %RETIRAR
%%
figure;
subplot 211
plot(1:length(signal(:,1)),signal(:,4),'b-');
subplot 212
plot(1:length(signal(:,1)),signal(:,3),'y-');
%% load train database 54 people (or 11)
abp_train=[];
ppg_train=[];

%%

[tm, signal]=rdsamp('mimicdb/410/',[],50000);
ppg_train=cat(2,ppg_train,signal(:,5));
abp_train=cat(2,abp_train,signal(:,3));

[tm, signal]=rdsamp('mimicdb/438/',[],50000);
ppg_train=cat(2,ppg_train,signal(:,5));
abp_train=cat(2,abp_train,signal(:,4));

[tm, signal]=rdsamp('mimicdb/437/',[],50000);
ppg_train=cat(2,ppg_train,signal(:,5));
abp_train=cat(2,abp_train,signal(:,4));

[tm, signal]=rdsamp('mimicdb/276/',[],50000);
ppg_train=cat(2,ppg_train,signal(:,6));
abp_train=cat(2,abp_train,signal(:,3));

[tm, signal]=rdsamp('mimicdb/281/',[],50000);
ppg_train=cat(2,ppg_train,signal(:,7));
abp_train=cat(2,abp_train,signal(:,3));

[tm, signal]=rdsamp('mimicdb/408/',[],50000);
ppg_train=cat(2,ppg_train,signal(:,5));
abp_train=cat(2,abp_train,signal(:,4));

[tm, signal]=rdsamp('mimicdb/401/',[],50000);
ppg_train=cat(2,ppg_train,signal(:,4));
abp_train=cat(2,abp_train,signal(:,3));

[tm, signal]=rdsamp('mimicdb/284/',[],50000);
ppg_train=cat(2,ppg_train,signal(:,5));
abp_train=cat(2,abp_train,signal(:,4));

[tm, signal]=rdsamp('mimicdb/240/',[],50000);
ppg_train=cat(2,ppg_train,signal(:,5));
abp_train=cat(2,abp_train,signal(:,3));

[tm, signal]=rdsamp('mimicdb/252/',[],50000);
ppg_train=cat(2,ppg_train,signal(:,5));
abp_train=cat(2,abp_train,signal(:,4));

[tm, signal]=rdsamp('mimicdb/230/',[],50000);
ppg_train=cat(2,ppg_train,signal(:,5));
abp_train=cat(2,abp_train,signal(:,3));

[tm, signal]=rdsamp('mimicdb/237/',[],50000);
ppg_train=cat(2,ppg_train,signal(:,5));
abp_train=cat(2,abp_train,signal(:,3));

[tm, signal]=rdsamp('mimicdb/221/',[],50000);
ppg_train=cat(2,ppg_train,signal(:,4));
abp_train=cat(2,abp_train,signal(:,3));

[tm, signal]=rdsamp('mimicdb/253/',[],50000);
ppg_train=cat(2,ppg_train,signal(:,6));
abp_train=cat(2,abp_train,signal(:,4));

% [tm, signal]=rdsamp('mimicdb/226/',[],50000);
% ppg_train=cat(2,ppg_train,signal(:,4));
% abp_train=cat(2,abp_train,signal(:,3));

[tm, signal]=rdsamp('mimicdb/224/',[],50000);
ppg_train=cat(2,ppg_train,signal(:,4));
abp_train=cat(2,abp_train,signal(:,3));

[tm, signal]=rdsamp('mimicdb/211/',[],50000);
ppg_train=cat(2,ppg_train,signal(:,4));
abp_train=cat(2,abp_train,signal(:,3));

[tm, signal]=rdsamp('mimicdb/225/',[],50000);
ppg_train=cat(2,ppg_train,signal(:,5));
abp_train=cat(2,abp_train,signal(:,4));

%% PHILIPS DATA


ppg_train=cat(2,ppg_train,signal(:,5));
abp_train=cat(2,abp_train,signal(:,3));

%% ------------------------------PRE PROCESSING ---------------------------------
%powerline interference:freq. fundamental e harmonicos (notch filter)



%baseline drift or motion artifact

%% SAVING DATA
save('DADOS.mat','abp_train','ppg_train')
%% LOADING DATA
load('DADOS.mat')
%% add philips data
load('#1.mat')

abp_train_new=measure(2).data(6.50001*10^5:7*10^5)';
ppg_train_new=measure(3).data(6.50001*10^5:7*10^5)';
%%
ppg_train=cat(2,ppg_train,ppg_train_new);
abp_train=cat(2,abp_train,abp_train_new);
%% -------------------------------- TRAIN GROUP --------------------------------
n_fs=8; %new sampling frequency
tam=18;
abp_group=cell(tam,1);
ppga_group=cell(tam,1);
rr_group=cell(tam,1);
sbp_fin=cell(tam,1);
ppga_fin=cell(tam,1);
diastolic_fin=cell(tam,1);
rr_fin=cell(tam,1);
locaisppg=cell(tam,1);
locaisabp=cell(tam,1);
locaisdiastolic=cell(tam,1);
rrsbp_fin=cell(tam,1);
b_a=cell(tam,1);
c_a=cell(tam,1);
RI=cell(tam,1);
as=cell(tam,1);
bs=cell(tam,1);
cs=cell(tam,1);
sysppg=cell(tam,1);
CT=cell(tam,1);
DELTAT=cell(tam,1);

for j=1:size(ppg_train,2)
    sys_temp=[];
    diast_temp=[];
    ppga_temp=[];
    abp_temp=[];
    rr_temp=[];
    [picoppg,localppg]=findpeaks(ppg_train(:,j),'MinPeakDistance',60,'MinPeakWidth',15); % distance of 50 in order to avoid capturing small peaks
    sysppg{j,1}=cat(1,sysppg{j,1}, picoppg(3:end));

    ppg_train(:,j)=fillmissing(ppg_train(:,j));
    locaisppg{j,1}=localppg(2:end,1);
%     figure;
    for i=3:length(localppg) %starts in 3 in order to avoid trespassing the signal 50000 samples, given the considerated intervals
        interval=localppg(i-1):localppg(i);
        
        %min global
        [minimo,local_min]= min(ppg_train(interval,j));
        
        %RR
        rr_i=localppg(i)-localppg(i-1);
        rr_temp=[rr_temp rr_i];
        
        %AMP PPG
        amp_ppg=ppg_train(localppg(i),j)-minimo;
        ppga_temp=[ppga_temp amp_ppg];
        
        
        %AMP SBP
        new_int=[localppg(i-1)-30:localppg(i-1)+30];
        [maximo,local_max]= findpeaks(abp_train(new_int,j),'MinPeakWidth',2,'SortStr','descend','NPeaks',1);
        if isempty(maximo)
            [maximo,local_max]=max(abp_train(new_int,j));
        end
        
        [minimo,local_min]= min(abp_train(interval,j)); %new_int
        locaisabp{j,1}=cat(1,locaisabp{j,1},localppg(i-1)-30+local_max);
        
        amp=maximo-minimo;
        abp_temp=[abp_temp amp];
        
        %2ND LEVEL OF PPG FEATURES
        interval=[localppg(i-1)-25:localppg(i-1)+50];
        
        %PPG PULSE WIDTH
        %widthloc=find(ppg_train(interval,j)>0.5*amp_ppg);
%         
%         plot(interval,ppg_train(interval,j),'b-');
%        
%         hold on;
        
        %derivatives of PPG event
        second=(diff(diff(ppg_train(interval,j))));
        
        %clear out 2nd derivative
        wc2=17;
        fs=125;
        ordem=5;
        second=process(ordem,fs,1,wc2,second);
        delay=5;
        for k=1:length(second)-delay
           second(k)=second(k+delay); 
        end
%         plot(interval(1:end-2),second*100,'r-');
%         hold on;
        
        %first=diff(ppg_train(interval,j));
        
        % %%
%         figure
%         % pwelch(second);
%         pwelch(fil,[],[],[],125,'twosided');
%         %%
%         figure
%         wc2=17;
%         ecg=second;
%         fs=125;
%         ordem=4;
%         fil=process(ordem,fs,1,wc2,ecg);
%         plot(1:length(second),fil);
        % %%
        [picosec,localsec]=findpeaks(second,'MinPeakWidth',3);
        if (length(picosec)<2)
            [picosec,localsec]=findpeaks(second);
            if (length(picosec)<2)
                ai=interval(1)+aip; ci=interval(1)+cip;
            else
                ai=interval(1)+localsec(1);
                ci=interval(1)+localsec(2);
                aip=localsec(1); cip=localsec(2);
            end
        else
            ai=interval(1)+localsec(1);
            ci=interval(1)+localsec(2);
            aip=localsec(1); cip=localsec(2);
        end
        as{j,1}=cat(1,as{j,1},ai);
        cs{j,1}=cat(1,cs{j,1},ci);
%         
%         plot(ai,ppg_train(ai,j),'y*');
%         hold on;
% 
%         plot(ci,ppg_train(ci,j),'r*');
%         hold on;
        
        [picosec1,localsec1]=findpeaks(-second,'MinPeakWidth',5); % distance of 50 in order to avoid capturing small peaks
        if (length(picosec1)<1)
            [picosec1,localsec1]=findpeaks(-second,'MinPeakWidth',3);
            if (length(picosec1)<1)
                bi=interval(1)+bip;
            else
                bi=interval(1)+localsec1(1);
                bip=localsec1(1);
            end
        else
            bi=interval(1)+localsec1(1);
            bip=localsec1(1);
        end
        bs{j,1}=cat(1,bs{j,1},bi);
        
%         plot(bi,ppg_train(bi,j),'g*');
%         hold on;
%         legend('PPG','2nd PPG','a','c','b');
%         ylabel('Amplitude')
%         xlabel('time (samples)')
%         
%         pause;
%         hold off;
        %if (ai
        %SI related indexes
        b_a{j,1}=cat(1,b_a{j,1},ppg_train(bi,j)/ppg_train(ai,j));
        c_a{j,1}=cat(1,c_a{j,1},ppg_train(ci,j)/ppg_train(ai,j));
        
        %AMP DBP
        interval=localppg(i-1)-30+local_max+5:localppg(i-1)-30+local_max+40; %INTERVALO PARA ENCONTRAR A POSIÇAO DO PICO DIASTOLICO DO BP
        
%         plot(1:length(interval),100*ppg_train(interval,j),'b-');
%         legend('PPG');
%         hold on;
%         plot(1:length(interval),abp_train(interval,j),'g-',1:length(second),100*second,'y-');
%         legend('ABP','2nd der');

        % 2nd derivative ABP --------------------------------------------
        second=(diff(diff(abp_train(interval,j))));
        wc2=17;
        fs=125;
        ordem=5;
        second=process(ordem,fs,1,wc2,second);
        delay=5;
        for k=1:length(second)-delay
           second(k)=second(k+delay); 
        end
        
        % 1st derivative ABP --------------------------------------------

        first=diff(abp_train(interval,j));
        wc2=17;
        fs=125;
        ordem=5;
        first=process(ordem,fs,1,wc2,first);
        delay=5;
        for k=1:length(first)-delay
           first(k)=first(k+delay); 
        end
        
        compara= find((second>first(1:end-1))) ;
        compara2= find(((second<first(1:end-1)))) ;
        
        loc1=find(first<0);
        loc2=find(diff(compara)>1);
        loc3=find(diff(compara2)>1);
        local_diastolic=[];
        diastolic=[];
        [minimo,local_min]= min(abp_train(interval,j)); %new_int

        if ( isempty(loc2)==0 )
            if( first(compara(loc2(1)+1))<=0)
                %                 1
                local_diastolic=interval(1)+compara(loc2(1)+1);
                diastolic=abp_train(local_diastolic);
%                                 plot(compara(loc2(1)+1),abp_train(local_diastolic,j),'k*');
%                                 hold on;
            end
        end
        if ( isempty(loc3)==0 & isempty(local_diastolic) )
            if( first(compara2(loc3(1)+1))<0)
                %                 2
                local_diastolic=interval(1)+compara2(loc3(1)+1);
                diastolic=abp_train(local_diastolic,j);
%                 plot(compara2(loc3(1)+1),abp_train(local_diastolic,j),'b*');
%                 hold on;
            end
        end
        if (length(loc1)>1 & isempty(local_diastolic))
            %                   3
            diferenca=abs(first(loc1(1:end-1))-second(loc1(1:end-1)));
            [diastolic,temp_pos]=min(diferenca);
            local_diastolic=interval(1)+temp_pos(1);
%             plot(temp_pos,abp_train(local_diastolic,j),'r*');
%             hold on;
        end
        if ( length(loc1)<2 & isempty(loc2) & isempty(loc3) & isempty(local_diastolic))
            
            [diastolic,local_diastolic_t]=min(abp_train(interval,j));
            local_diastolic=interval(1)+local_diastolic_t;
%             plot(local_diastolic,abp_train(local_diastolic,j),'y*');
%             hold on;
        end
        %         5
%         plot(1:length(interval),abp_train(interval,j),'g-',1:length(second),100*second,'y-',1:length(first),100*first,'r-');
%         legend('picos diastolicos','ABP','2nd der','1st der');
%         xlabel('Samples');ylabel('Blood pressure (mmHg)');
%         pause
%         hold off;
        
        if (isempty(local_diastolic) & isempty(diastolic))
            1
            i
        end
        diast_temp=[diast_temp diastolic];
        locaisdiastolic{j,1}=cat(1,locaisdiastolic{j,1},local_diastolic);
        
    end
    
    abp_resampled=resample(abp_temp,n_fs,1); %experimentar com 4 ou UPSAMPLE!!!!!!!!!!!!!!!!!!!!!!
    abp_group{j,:}=abp_resampled;
    sbp_fin{j,:}=abp_temp;
    
    %         Normalization (mean centered and resampling at 16Hz: valor arbitrario)
    %         ppga_centered=ppga_temp-mean(ppga_temp);
    
    ppga_resampled=resample(ppga_temp./max(ppga_temp),n_fs,1); %experimentar com 4
    ppga_group{j,:}=ppga_resampled;
    ppga_fin{j,:}=ppga_temp./max(ppga_temp); %normalization according to max value
    
    %     % Normalization (mean centered and resampling at 16Hz: valor arbitrario)
    %     %     rr_centered=rr_temp-mean(rr_temp);
    
    rrc_resampled=resample(rr_temp,n_fs,1); %experimentar com 4
    rr_group{j,:}=rrc_resampled;
    rr_fin{j,:}=rr_temp;
    
    %RR series of BP
    rrsbp_fin{j,:}=diff(locaisabp{j,1});
    
    %Diastolic peaks of BP
    diastolic_fin{j,:}=diast_temp;
    
    %RI
    RI{j,1}=ppg_train(cs{num,1},num)'./ppga_fin{j,:}';
%     
%     %1st derivative features
%     %CT
%     CT{j,1}=locaisppg{j,1}(1:end-1)-ppg_train(as{num,1},num);
%     
%     %DELTAT
%     DELTAT{j,1}=ppg_train(cs{num,1},num)-locaisppg{j,1}(1:end-1);
%     %num=12; EXEMPLO DE SINAL PESSIMO CHEIO DE SATURAÇAO

 end
%%
for num=17:tam
%     num=1;
    second=(diff(diff(ppg_train(:,num))));
    %clear out 2nd derivative
    wc2=17;
    fs=125;
    ordem=5;
    second=process(ordem,fs,1,wc2,second);
    delay=5;
    for k=1:length(second)-delay
       second(k)=second(k+delay); 
    end
    
    figure;
    subplot 211
    plot(locaisdiastolic{num,1},abp_train(locaisdiastolic{num,1},num),'r*',locaisabp{num,1},abp_train(locaisabp{num,1},num),'b*',1:length(abp_train(:,num)),abp_train(:,num),'g-');
    ylabel('mmHg')
    xlabel('Samples')
    legend('Diastolic peaks', 'Systolic peaks', 'BP')
    subplot 212
    plot(1:length(second),100*second,'r-',as{num,1},ppg_train(as{num,1},num),'r*',bs{num,1},ppg_train(bs{num,1},num),'y*',cs{num,1},ppg_train(cs{num,1},num),'b*',locaisppg{num,1},ppg_train(locaisppg{num,1},num),'g*',1:length(ppg_train(:,num)),ppg_train(:,num),'b-');
    ylabel('Amplitude')
    xlabel('Samples')
    legend('2nd der','a','b','c','sistolicos','PPG')
    pause
    
end
%% 1) 1st LEVEL of features (inspired in SPI, a RR and PPGA based algorithm: )

% LEAVE ONE OUT CROSS VALIDATION
level=1;

if level==1
    
    fs=n_fs;
    figure('Name','1st level')
    for l=1:size(ppga_group,1)
        
        ppga_train=cell(tam,1);
        rr_train=cell(tam,1);
        
        for k=1:size(ppga_group,1)
            if(k==l)
                ppga=ppga_group{k,:};
                rr=rr_group{k,:};
            else
                ppga_train{k,1}=ppga_fin{k,1};%ppga_group
                rr_train{k,1}=rr_fin{k,1};%ppga_group
            end
        end
        % ------------------------------ PPG and HBI normalization --------------------
        
        task='train';
        [comb_ppga_train,ppga_mean,ppga_std] = normalization(ppga_train,' ',' ','',fs,task);
        [comb_rr_train,rr_mean,rr_std] = normalization(rr_train,' ',' ','',fs,task);
        
        % artificial train
        % comb_ppga_train=[0 1];
        % ppga_mean=mean(ppga_temp);
        % ppga_std=std(ppga_temp);
        % comb_rr_train=[0 1];
        % rr_mean=mean(rr_temp);
        % rr_std=std(rr_temp);)
        % ppga=ppga_temp;
        % rr=rr_temp;
        % fs=n_fs;
        %
        task='test';
        fs=8;
        
        %funçao para reavaliaçao dos pesos individuais/grupo
        V=[1,0.3];
        X=[0,0.7];
        Xq=0:0.7/(5*60*fs):0.7;
        gr_val_fun=[interp1(X,V,Xq,'spline') 0.3*ones(1,15000)];
        
        SPI_values=[];
        
        % simulation of incision and intubation (decrease in PPGA 1 min before the ending of analysis)
        
        % ppga(1,end-n_fs*60.2:end-n_fs*60)=0.7;
        % ppga(1,n_fs*50:n_fs*50.2)=0.7;
        
        for i=1:n_fs:length(ppga)-1% signal size
            % pacientes do validation set. amostra i de todos os pacientes
            [comb_ppga_test,ppga_media_ind] = normalization(ppga(1,i:i+n_fs-1),round(i/n_fs)+1,comb_ppga_train,ppga_mean,fs,task);
            [comb_rr_test,rr_media_ind] = normalization(rr(1,i:i+n_fs-1),round(i/n_fs)+1,comb_rr_train,rr_mean,fs,task);
            
            ppga_mean=ppga_media_ind;
            rr_mean=rr_media_ind;
            comb_rr_train=[1-gr_val_fun(round(i/n_fs)+1) gr_val_fun(round(i/n_fs)+1)];
            comb_ppga_train=[1-gr_val_fun(round(i/n_fs)+1) gr_val_fun(round(i/n_fs)+1)];
            
            %valores normalizados de 0 a 100
            spi_ppga = normcdf(ppga(1,i),ppga_media_ind,ppga_std)*100;
            spi_rr= normcdf(rr(1,i),rr_media_ind,rr_std)*100;
            
            coef_HBI=0.3; %test with +-0.22
            coef_ppga=0.7;
            spi=100-(coef_ppga*spi_ppga + coef_HBI*spi_rr); % impressao do spi para cada um dos pacientes
            % fprintf ('SPI= %4.2f \n',spi);
            if i==1
                SPI_values=cat(2,SPI_values,spi);
            end
            
            if i>1
                SPI_values=cat(2,SPI_values,spi*0.5+SPI_values(round(i/n_fs))*0.5);
            end
        end
        
        %100: stress
        %0: sem stress
        
        % ------------------------ ABP FEATURES
        sdbp_fin{l,1}=sbp_fin{l,1}-diastolic_fin{l,1};
        
        fs=125;
        %     subplot(711)
        %     plot((1:50000)/fs,abp_train(:,l),'b-');
        %     xlabel('Time (seconds)')
        %     ylabel('ABP (mmHg)')
        %     legend('ABP')
        
        subplot(511)
        plot((locaisabp{l,1})/fs,sdbp_fin{l,1},'b-');
        xlabel('Time (seconds)')
        ylabel('SBP (mmHg)')
        legend('SBP')
        
        subplot(512)
        plot((locaisppg{l,1}(2:end))/fs,SPI_values,'r-');
        xlabel('Time (seconds)')
        ylabel('SPI')
        legend('SPI')
        
        subplot(513)
        plot((locaisabp{l,1}(1:end-1))/fs,rrsbp_fin{l,1},'b-');
        xlabel('Time (seconds)')
        ylabel('SBP-DBP (mmHg)')
        legend('SBP-DBP')
        
        subplot(514)
        plot(sdbp_fin{l,1},SPI_values,'b*');
        xlabel('SBP (mmHg)')
        ylabel('SPI')
        
        subplot(515)
        plot(rrsbp_fin{l,1},SPI_values(1:end-1),'b*');
        xlabel('SBP-DBP (mmHg)')
        ylabel('SPI')
        
        fprintf ('sample= %2.0f \n',l);
        fprintf ('SBP vs SPI');
        [R,P]=corrcoef(sdbp_fin{l,1},SPI_values) %PPGA
        fprintf ('RR vs SPI');
        [R1,P1]=corrcoef(rrsbp_fin{l,1},SPI_values(1:end-1)) %ABP
        
        % ------------------------ PPG FEATURES
        pause
        %
        %     subplot(511)
        %     plot((locaisppg{l,1}(2:end))/fs,ppga_fin{l,1},'g-');
        %     xlabel('Time (seconds)')
        %     ylabel('PPGA (amplitude)')
        %     legend('PPGA')
        %
        %     subplot(512)
        %     plot((locaisppg{l,1}(2:end))/fs,SPI_values,'r-');
        %     xlabel('Time (seconds)')
        %     ylabel('SPI')
        %     legend('SPI')
        %
        %     subplot(513)
        %     plot((locaisppg{l,1}(2:end))/fs,rr_fin{l,1}/fs,'b-');
        %     xlabel('Time (seconds)')
        %     ylabel('RR time series (seconds)')
        %     legend('RR')
        %
        %     subplot(514)
        %     plot(ppga_fin{l,1},SPI_values,'b*');
        %     xlabel('PPGA (amplitude)')
        %     ylabel('SPI')
        %
        %     subplot(515)
        %     plot(rr_fin{l,1},SPI_values,'b*');
        %     xlabel('RR time series (seconds)')
        %     ylabel('SPI')
        %
        %
        % %     subplot(717)
        % %     plot(ppga_fin{l,1},abp_fin{l,1},'b*');
        % %     xlabel('PPGA (amplitude)')
        % %     ylabel('SBP (mmHg)')
        %
        %     fprintf ('sample= %2.0f \n',l);
        %     fprintf ('PPGA');
        %     [R,P]=corrcoef(ppga_fin{l,1},SPI_values) %PPGA
        %     fprintf ('RR');
        %     [R1,P1]=corrcoef(rr_fin{l,1},SPI_values) %RR
        %     pause
        %checar se P1 e P possuem valores menores que 0.0005, validando
        %correlaçao entre variaveis
        
    end
    
    % ------------------------------------------- 2) 2nd LEVEL of features (inspired in lower level features in PPG)
    
    % LEAVE ONE OUT CROSS VALIDATION
elseif level==2
    corrs=[];
    for l=1:17
%         figure('Name','RI')
        
        %100: stress
        %0: sem stress
        
        % ------------------------ ABP FEATURES
        sdbp_fin{l,1}=sbp_fin{l,1}-diastolic_fin{l,1};
        
%         fs=125;
%         
%         subplot(311)
%         plot(RI{l,1},diastolic_fin{l,1},'r*');
%         xlabel('RI')
%         ylabel('Diastolic')
%         
%         subplot(312)
%         plot(RI{l,1},sdbp_fin{l,1},'b*');
%         xlabel('RI')
%         ylabel('SIS-DIAST')
%         
%         subplot(313)
%         plot(RI{l,1},sbp_fin{l,1},'b*');
%         ylabel('SBP (mmHg)')
%         xlabel('RI')
        
%         fprintf ('sample= %2.0f \n',l);     
        
%         fprintf ('DIAST vs RI');
        [R1]=corrcoef(DELTAT{l,1},sbp_fin{l,1}); %DBP
        
%         fprintf ('SDBP vs RI');
        [R2]=corrcoef(DELTAT{l,1},diastolic_fin{l,1}); %SBP
        
%         fprintf ('SBP vs RI');
        [R3]=corrcoef(DELTAT{l,1},sdbp_fin{l,1}); %SBP
        
        %         subplot(511)
        %         plot((locaisabp{l,1})./fs,sdbp_fin{l,1},'b-');
        %         xlabel('Time (seconds)')
        %         ylabel('SBP (mmHg)')
        %         legend('SBP')
        %
        %         subplot(512)
        %         plot((locaisppg{l,1}(2:end-1))./fs,SPI_values,'r-');
        %         xlabel('Time (seconds)')
        %         ylabel('SPI')
        %         legend('SPI')
        %
        %         subplot(513)
        %         plot((locaisdiastolic{l,1})./fs,diastolic_fin{l,1},'b-');
        %         xlabel('Time (seconds)')
        %         ylabel('DBP (mmHg)')
        %         legend('DBP')
        %
        %         subplot(514)
        %         plot(sdbp_fin{l,1}(2:end),SPI_values,'b*');
        %         xlabel('SBP (mmHg)')
        %         ylabel('SPI')
        %
        %         subplot(515)
        %         plot(diastolic_fin{l,1}(2:end),SPI_values,'b*');
        %         xlabel('DBP (mmHg)')
        %         ylabel('SPI')
        %
        %         fprintf ('sample= %2.0f \n',l);
        %         fprintf ('SBP');
        %         [R]=corrcoef(sdbp_fin{l,1}(2:end),SPI_values) %SBP
        %         fprintf ('DBP');
        %         [R1]=corrcoef(diastolic_fin{l,1}(2:end),SPI_values) %DBP
        %
        % ------------------------ SBP comparison
%         pause
        
%         figure('Name', 'SBP AND DBP');
%         subplot(311)
%         plot(b_a{l,1},c_a{l,1},'b*');
%         xlabel('B/A')
%         ylabel('C/A')
        
%         fprintf ('C/A vs B/A');
        [R4]=corrcoef(b_a{l,1},c_a{l,1}); %C/A
%         
%         subplot(312)
%         plot(b_a{l,1},sbp_fin{l,1},'b*',b_a{l,1},diastolic_fin{l,1},'g*');
%         xlabel('B/A')
%         ylabel('BP (mmHg)')
%         legend('SBP','DBP')
         
%         subplot(313)
%         plot(c_a{l,1},sbp_fin{l,1},'b*',c_a{l,1},diastolic_fin{l,1},'g*');
%         xlabel('C/A')
%         ylabel('BP (mmHg)')
%         legend('SBP','DBP')
          
%         fprintf ('sample= %2.0f \n',l);
%         fprintf('SISTOLIC ');
%         fprintf ('C/A');
%         [R5]=corrcoef(sbp_fin{l,1},c_a{l,1}); %C/A
%         fprintf ('B/A');
%         [R6]=corrcoef(sbp_fin{l,1},b_a{l,1}); %B/A
        
        % ------------------------ Diastolic comparison
         
%         fprintf ('sample= %2.0f \n',l);
%         fprintf('DIASTOLIC ');
%         fprintf ('C/A');
%         [R7]=corrcoef(diastolic_fin{l,1},c_a{l,1}); %C/A
%         fprintf ('B/A');
%         [R8]=corrcoef(diastolic_fin{l,1},b_a{l,1}); %B/A
%         pause
        corrs=cat(2,corrs,[R1(1,2);R2(1,2);R3(1,2);R4(1,2)]);%;R5(1,2);R6(1,2);R7(1,2);R8(1,2)]);
    end
    T=array2table(corrs,'RowNames',{'R1'; 'R2'; 'R3' ;'R4'})% ;'R5'; 'R6' ;'R7'; 'R8'})
end