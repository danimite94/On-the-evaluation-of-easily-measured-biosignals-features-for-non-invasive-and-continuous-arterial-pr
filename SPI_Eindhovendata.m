%% load train database 3 people
abp_group=cell(3,1);
ppg_group=cell(3,1);

%% PHILIPS DATA

addpath(genpath('C:\Users\310277823\Desktop\tese\matlab\'));

signal=1;
% Philips DATA
load('#1.mat')

% SPECIFIC REGION IN ANALYSIS
surgery=measure(7).time(measure(7).data()>0);
% example=measure(7).time(measure(7).time()>6000);
begin=surgery(1);
fim=surgery(end);

ex_begin=begin;
ex_end=fim;

%MANUAL DETECTION OF NOISY REGIONS
restrict=[begin 3005; 5890 5905; 10263 10267; 10606 10608; 9721 9723; 9853 9855; 9637 9639; 9643 9644; 9661 9666; 3040 3070; 3210 3281; 3532 3550; 3671 3700; 7003 7004; 7125 7127; 7681 7683; 7933 7935; 8400 8430; 8490 8517; 8631 8632; 8694 8771; 9041 9043; 9132 9133; 9156 9157; 9190 9191; 9322 9323; 9499 9500; 9516 9519; 9746 9748; 9887 9888; 9908 9911; 10587 10593; 10599 10601; 10625 10637; 10735 10743; 10995 11025; 11042 11044; 11120 11130; 11355 11360; 11207 11228; 10860 10875; 6333 6335; 6238 6241; 6120 6124; 5652 5656; 5580 5585; 5205 5280; 4612 4615; 4450 4454; 4654 4657; 4733 4737; 6817 6837; 6317 6319; 6678 6680; 6455 6460];

ex_clean_abp=measure(2).data(find(measure(2).time()>ex_begin & measure(2).time()<ex_end))';
t1=measure(2).time(find(measure(2).time()>ex_begin & measure(2).time()<ex_end));
ex_clean_ppg_DC=measure(3).data(find(measure(3).time()>ex_begin & measure(3).time()<ex_end))'; %DC+AC
t2=measure(3).time(find(measure(3).time()>ex_begin & measure(3).time()<ex_end));
ex_clean_ECG=measure(1).data(find(measure(1).time()>ex_begin & measure(1).time()<ex_end))';
t3=measure(1).time(find(measure(1).time()>ex_begin & measure(1).time()<ex_end));

% DC removal and y-axis flip
fs=125;
wc1=0.05;
ordem=5;
ppg_train_new=process(ordem,fs,wc1,1,ex_clean_ppg_DC,2);
ppg_train_minus=-ppg_train_new;
ppg_train=ppg_train_minus;

ppg_group{1,1}=ppg_train;
abp_group{1,1}=ex_clean_abp;

signal=3;
% Philips DATA
load('#3.mat')

% SPECIFIC REGION IN ANALYSIS
surgery=measure(7).time(measure(7).data()>0);
% example=measure(7).time(measure(7).time()>6000);
begin=surgery(1);
fim=surgery(end);

ex_begin=begin;
ex_end=fim;

%MANUAL DETECTION OF NOISY REGIONS
restrict=[begin 2140; 2151 2160; 2174 2180; 2194 2199; 2232 2234; 2242 2250; 2284 2290; 2298 2322; 2410 2423; 2436 2480; 2581 2585; 3187 3189; 3242 3244; 3344 3401; 3536 3591; 3613 3690; 3824 3840; 3939 3943; 4171 4172; 4783 4787; 4806 4808; 4988 4990; 5146 5181; 5940 5943; 6320 6328; 6405 6408; 6456 6467; 6820 6823; 6865 6868; 6942 7035; 7443 7448; 7698 7704; 8464 8496; 8693 8791; 9527 9555; 10540 10580; 10650 10680; 10800 10810; 10920 10950; 11034 11037; 11553 11557; 11654 11660; 12050 12055; 12190 12195; 12337 12375; 12520 12524; 12695 12742; 12867 12890; 12910 12920; 12948 12985; 13057 13060; 13110 13440; 13478 13715; 13750 13762]; %#3

ex_clean_abp=measure(2).data(find(measure(2).time()>ex_begin & measure(2).time()<ex_end))';
t1=measure(2).time(find(measure(2).time()>ex_begin & measure(2).time()<ex_end));
ex_clean_ppg_DC=measure(3).data(find(measure(3).time()>ex_begin & measure(3).time()<ex_end))'; %DC+AC
t2=measure(3).time(find(measure(3).time()>ex_begin & measure(3).time()<ex_end));
ex_clean_ECG=measure(1).data(find(measure(1).time()>ex_begin & measure(1).time()<ex_end))';
t3=measure(1).time(find(measure(1).time()>ex_begin & measure(1).time()<ex_end));

% DC removal and y-axis flip
fs=125;
wc1=0.05;
ordem=5;
ppg_train_new=process(ordem,fs,wc1,1,ex_clean_ppg_DC,2);
ppg_train_minus=-ppg_train_new;
ppg_train=ppg_train_minus;

ppg_group{2,1}=ppg_train;
abp_group{2,1}=ex_clean_abp;

signal=5;
% Philips DATA
load('#5.mat')

% SPECIFIC REGION IN ANALYSIS
surgery=measure(7).time(measure(7).data()>0);
% example=measure(7).time(measure(7).time()>6000);
begin=surgery(1);
fim=surgery(end);

ex_begin=begin;
ex_end=fim;

%MANUAL DETECTION OF NOISY REGIONS
restrict=[2579 2600; 9096 9450; 7835 7840; 1245 1430; 1693 1695; 1929 1946; 2437 2444; 2677 2686; 3330 3333; 3700 4603; 4940 4945; 4989 4992; 5396 5400; 5848 5998; 6666 6669; 6927 7045; 8362 8397; 8505 8640; 8865 8943; 9538 9540; 9779 9781; 10455 10530; 10795 10895; 11010 11045; 11275 11320; 11390 11400; 11430 11450; 11498 11503; 11547 ex_end]; %; #5

ex_clean_abp=measure(2).data(find(measure(2).time()>ex_begin & measure(2).time()<ex_end))';
t1=measure(2).time(find(measure(2).time()>ex_begin & measure(2).time()<ex_end));
ex_clean_ppg_DC=measure(3).data(find(measure(3).time()>ex_begin & measure(3).time()<ex_end))'; %DC+AC
t2=measure(3).time(find(measure(3).time()>ex_begin & measure(3).time()<ex_end));
ex_clean_ECG=measure(1).data(find(measure(1).time()>ex_begin & measure(1).time()<ex_end))';
t3=measure(1).time(find(measure(1).time()>ex_begin & measure(1).time()<ex_end));

% DC removal and y-axis flip
fs=125;
wc1=0.05;
ordem=5;
ppg_train_new=process(ordem,fs,wc1,1,ex_clean_ppg_DC,2);
ppg_train_minus=-ppg_train_new;
ppg_train=ppg_train_minus;

ppg_group{3,1}=ppg_train;
abp_group{3,1}=ex_clean_abp;
%% --------------------------------------------------- PPG and BP features --------------------------------------------------
tot_tam=3;
n_fs=8;
rr_group=cell(tot_tam,1);
sbp_fin=cell(tot_tam,1);
SIG1_amp_fin=cell(tot_tam,1);
SIG2_amp_fin=cell(tot_tam,1);
DC_ppga_fin=cell(tot_tam,1);
locaisDC_ppga_fin=cell(tot_tam,1);
DC_max=cell(tot_tam,1);
diastolic_fin=cell(tot_tam,1);
rr_fin=cell(tot_tam,1);
locaisspSIG1=cell(tot_tam,1);
locaisspSIG2=cell(tot_tam,1);
locaisdiastSIG2=cell(tot_tam,1);
rrsbp_fin=cell(tot_tam,1);
b_a=cell(tot_tam,1);
c_val=cell(tot_tam,1);
b_val=cell(tot_tam,1);
a_val=cell(tot_tam,1);
e_val=cell(tot_tam,1);
RI_temp=cell(tot_tam,1);
SI_temp=cell(tot_tam,1);
es=cell(tot_tam,1);
as=cell(tot_tam,1);
bs=cell(tot_tam,1);
cs=cell(tot_tam,1);
dicnotch_temp=cell(tot_tam,1);
locaisdicnotch=cell(tot_tam,1);
diast_temp=cell(tot_tam,1);
locaisdiast=cell(tot_tam,1);
locaisdiastSIG1=cell(tot_tam,1);
TEMPPAT=cell(tot_tam,1);
PAT_time=cell(tot_tam,1);
AC_DC=cell(tot_tam,1);
RR_fin=cell(tot_tam,1);

CORR_SHOW=0;

for tam=1:tot_tam
    
    [picoppg,localppg]=findpeaks(ppg_group{tam,1},'MinPeakDistance',60,'MinPeakWidth',15); % distance of 60 in order to avoid capturing small peaks
    [picoabp,localabp]=findpeaks(abp_group{tam,1},'MinPeakDistance',60,'MinPeakWidth',15,'MinPeakHeight',70); % distance of 60 in order to avoid capturing small peaks
    
    % ppg_train(:,j)=fillmissing(ppg_train(:,j));
    
    if CORR_SHOW==1
        figure;
    end
    
    SIG1=ex_clean_abp; %ex_clean_ECG
    SIG2=ppg_train;
    
    % -------------------------------------------------------------- Peak detection
    for i=2:length(localabp)-2
        interval1=localabp(i-1):localabp(i); %local1(i-1):local1(i)
        [minimo1,local_min1]= min(SIG1(interval1,tam));         %min global
        %-------------------------------------------------------------- RR
        %     rr_i=local1(i)-local1(i-1);
        
        %------------------------------------------------------------ AMP AC BP
        amp_SIG1=SIG1(localabp(i),tam)-minimo1;
        
        SIG1_amp_fin{tam,1}=cat(1,SIG1_amp_fin{tam,1},amp_SIG1);
        locaisspSIG1{tam,1}=cat(1,locaisspSIG1{tam,1},localabp(i));
        locaisdiastSIG1{tam,1}=cat(1,locaisdiastSIG1{tam,1},interval1(1)+local_min1);
        
        %     rr_temp=[rr_temp rr_i];
    end
    
    for i=2:length(localppg)-1
        %---------------------------------------------------------- PPG features
        %starts in 3 in order to avoid trespassing the signal 50000 samples, given the considerated intervals
        interval2=localppg(i-1):localppg(i); %local2(i-1):local2(i)
        [minimo2,local_min2]= min(SIG2(interval2,tam));
        amp_SIG2=SIG2(localppg(i),tam)-minimo2; %AC amplitude
        
        %----------------------------------------------------- PPG Diastolic value
        locaisdiastSIG2{tam,1}=cat(1,locaisdiastSIG2{tam,1},local_min2+interval2(1));
        locaisspSIG2{tam,1}=cat(1,locaisspSIG2{tam,1},localppg(i));
        SIG2_amp_fin{tam,1}=cat(1,SIG2_amp_fin{tam,1},amp_SIG2);
        
        % ------------------------------------------------------------ AMP DC PPG
        [picotemp,localtemp]=findpeaks(ex_clean_ppg_DC(interval2,1),'MinPeakDistance',50,'MinPeakWidth',15); % distance of 50 in order to avoid capturing small peaks
        if isempty(localtemp)
            [picotemp,localtemp]=max(ex_clean_ppg_DC(interval2,1));
        end
        [minimo2,local_min2]= min(ex_clean_ppg_DC(interval2,1));         %min global
        amp_ppg_DC=(picotemp(1)+minimo2(1))/2;
        
        DC_max{tam,1}=cat(1,DC_max{tam,1},picotemp);
        DC_ppga_fin{tam,1}=cat(1,DC_ppga_fin{tam,1},amp_ppg_DC);
        locaisDC_ppga_fin{tam,1}=cat(1,locaisDC_ppga_fin{tam,1},interval2(1)+localtemp);
        
        %------------------------------------------------------------------ 2ND LEVEL OF PPG FEATURES
        interval=localppg(i)-50:localppg(i)+50;
        
        if CORR_SHOW==1
            temp_inter=localppg(i)-50:localppg(i)+50;
            plot(temp_inter,ppg_train(temp_inter,tam),'b-');
            hold on;
        end
        
        %-------------------------------------------------------- a,b,c and e/DN of PPG
        second_ppg_spe=(diff(diff(ppg_train(interval,tam))));
        wc2=17;
        fs=125;
        ordem=5;
        second_ppg_spe=process(ordem,fs,1,wc2,second_ppg_spe,1);
        
        [picosec,localsec]=findpeaks(second_ppg_spe,'MinPeakWidth',3,'SortStr','descend');
        if (length(picosec)<1)
            [picosec,localsec]=findpeaks(second_ppg_spe,'MinPeakWidth',2,'SortStr','descend');
            if (length(picosec)<1) %if clause for the manual removed areas (yet to be removed)
                ai=interval(1)+aip;
                picosec=a_val{tam,1}(end);
            else
                ai=interval(1)+localsec(1);
                aip=localsec(1);
            end
        else
            ai=interval(1)+localsec(1);
            aip=localsec(1);
        end
        as{tam,1}=cat(1,as{tam,1},ai);
        a_val{tam,1}=cat(1,a_val{tam,1},picosec(1));
        
        [picosec1,localsec1]=findpeaks(-second_ppg_spe(localsec(1):end),'MinPeakWidth',3); % distance of 50 in order to avoid capturing small peaks
        if (length(picosec1)<1)
            [picosec1,localsec1]=findpeaks(-second_ppg_spe(localsec(1):end),'MinPeakWidth',1);
            if (length(picosec1)<1) %if clause for the manual removed areas (yet to be removed)
                bi=ai+bip;
                picosec1=-b_val{tam,1}(end);
                localsec1=-5; %noisy regions purposes (no implications in the results)
            else
                bi=ai+localsec1(1)-1;
                bip=localsec1(1);
            end
        else
            bi=ai+localsec1(1)-1;
            bip=localsec1(1);
        end
        bs{tam,1}=cat(1,bs{tam,1},bi);
        b_val{tam,1}=cat(1,b_val{tam,1},-picosec1(1));
        
        int=localsec(1)+localsec1(1);
        %----------------------------------------------------------diastolic peak in elgendi article
        interval=localppg(i)+10:localppg(i)+50; % INTERVAL TO FIND diastolic POSITION ("local2(i)+10" is to help finding the nearest point to zero (which will be the diastolic peak) )
        
        first_ppg=diff(ppg_train(interval,tam));
        wc2=17;
        fs=125;
        ordem=5;
        first_ppg=process(ordem,fs,1,wc2,first_ppg,1);
        
        loc1=find(first_ppg>0);
        compara=any(diff(loc1)>2); %when there is a gap in this vector, it means that f' went from positive to negative
        [vv,loc2]=min(abs(first_ppg));
        
        if ( compara==1 ) %when 1st derivative goes from positive to negative values
            pos=find(diff(loc1)>2);
            local_diast=interval(1)+loc1(pos);
            diast=ppg_train(local_diast,tam);
            
            if CORR_SHOW==1
                plot(local_diast,diast,'k*');
                hold on;
            end
        elseif ( compara==0 & isempty(loc2)==0 ) %when 1st derivative is the closest to the zero
            
            local_diast=interval(1)+loc2(1);
            diast=ppg_train(local_diast,tam);
            
            if CORR_SHOW==1
                plot(local_diast,diast,'r*');
                hold on;
            end
        end
        
        diast_temp{tam,1}=cat(1,diast_temp{tam,1},diast(1));
        locaisdiast{tam,1}=cat(1,locaisdiast{tam,1},local_diast(1));
        %--------------------------------------------------------------c,e/DICROTIC NOTCH (same position as "e" feature in elgendi)
        int2=local_diast(1)-interval(1)+(interval(1)-(localppg(i)-50))-1;%posicao do diastolico no intervalo para detecao do a,b
        
        if int2-int<4 %noisy regions purposes (no implications in the results)
            int=int-50;
        end
        
        [picosec2,localsec2]=findpeaks(second_ppg_spe(int:int2),'MinPeakWidth',3);
        if(length(picosec2)<2)
            [picosec2,localsec2]=findpeaks(second_ppg_spe(int:int2),'MinPeakWidth',1);
            if (length(picosec2)==1) %if clause for the manual removed areas (yet to be removed)
                ci=bi+localsec2(1);
                ei=bi+localsec2(1);
                e_val{tam,1}=cat(1,e_val{tam,1},picosec2(1));
            elseif (isempty(picosec2))
                ci=interval(1)+int;
                ei=interval(1)+int;
                e_val{tam,1}=cat(1,e_val{tam,1},e_val{tam,1}(end));
                picosec2(1)=c_val{tam,1}(end);
            else
                ci=bi+localsec2(1);
                ei=bi+localsec2(2);
                e_val{tam,1}=cat(1,e_val{tam,1},picosec2(2));
            end
        else
            ci=bi+localsec2(1);
            ei=bi+localsec2(2);
            e_val{tam,1}=cat(1,e_val{tam,1},picosec2(2));
        end
        
        if CORR_SHOW==1
            plot(ai,a_val{tam,1}(end)*100,'y*',bi,b_val{tam,1}(end)*100,'g*');
            hold on;
        end
        
        if length(picosec2)>2
            local_dicnotch=bi+localsec2(3);
            dicnotch=ppg_train(local_dicnotch,tam);
        else
            local_dicnotch=ei;
            dicnotch=ppg_train(local_dicnotch,tam);
        end
        
        es{tam,1}=cat(1,es{tam,1},ei);
        cs{tam,1}=cat(1,cs{tam,1},ci);
        c_val{tam,1}=cat(1,c_val{tam,1},picosec2(1));
        dicnotch_temp{tam,1}=cat(1,dicnotch_temp{tam,1},dicnotch);
        locaisdicnotch{tam,1}=cat(1,locaisdicnotch{tam,1},local_dicnotch);
        
        if CORR_SHOW==1
            interval=localppg(i)-50:localppg(i)+50; %INTERVAL TO FIND DIC NOTCH POSITION
            second_ppg=(diff(diff(ppg_train(interval,tam))));
            
            wc2=17;
            fs=125;
            ordem=5;
            second_ppg=process(ordem,fs,1,wc2,second_ppg,1);
            
            first_ppg=diff(ppg_train(interval,tam));
            wc2=17;
            fs=125;
            ordem=5;
            first_ppg=process(ordem,fs,1,wc2,first_ppg,1);
            
            plot(interval(2:end-1),second_ppg*100,'k-',ci,c_val{tam,1}(end)*100,'go',local_dicnotch,dicnotch,'k*'); %local_dicnotch,dicnotch,'go',
            legend('PPG','diast','a','b','2nd*100','c','DN(e)');
            xlabel('Samples');ylabel('Amplitude');
            pause
        end
        
    end
    
    %----------------------------------------------------------DN
    dicnotch_fin{tam,1}=dicnotch_temp{tam,1}-ppg_train(locaisdiastSIG2{tam,1},tam);
    
    %----------------------------------------------------------RI
    RI_temp{tam,1}=(diast_temp{tam,1}-ppg_train(locaisdiastSIG2{tam,1},tam))./SIG2_amp_fin{tam,1};
    
    %----------------------------------------------------------SI
    SI_temp{tam,1}=b_val{tam,1}./a_val{tam,1};
    
    %--------------------------------------------------------RR series BP
    rrsbp_fin{tam,1}=diff(locaisspSIG1{tam,1});
    
    %------------------------------------------------------diastolic BP peak
    diastolic_fin{tam,1}=SIG1(locaisdiastSIG1{tam,1},tam);
    
    %------------------------------------------------------systolic BP peak
    sbp_fin{tam,1}=SIG1(locaisspSIG1{tam,1},tam);
    
    %--------------------------------------------- Mean arterial pressure (MAP)
    MAP_temp{tam,1}=(SIG1(locaisspSIG1{tam,1},tam)+2.*SIG1(locaisdiastSIG1{tam,1},tam))./3;
    
    %--------------------------------------------- Pulse pressure (PP)
    PP_fin{tam,1}=sbp_fin{tam,1}-diastolic_fin{tam,1};
    
    %---------------------------------------------- Pulse transit time (PTT)
    % PTT{j,1}=locaisppg{j,1}-locaissbp{j,1};
    
    %---------------------------------------------- DELTA T
    DT{tam,1}=abs(t2(locaisspSIG2{tam,1})-t2(locaisdiast{tam,1}));
    
    % Normalized volume
    AC_DC{tam,1}=SIG2_amp_fin{tam,1}./DC_ppga_fin{tam,1};
    
    %---------------------------------------------- Pulse Arrival time (PAT)
    
    %Pan tompkins (R peak detection)
    [b,c,ECGrpeak_pos,d,e] = pan_tompkins_qrs(ex_clean_ECG,250,0);
    
    exemp_ECG=ECGrpeak_pos;
    % exemp_ECG_resamp=resample(exemp_ECG,8,1);
    
    RR_fin{tam,1}=diff(exemp_ECG)./250;
    
    ecg_PAT=t3(exemp_ECG);
    [ppg80_values,ppg80_time] = ppgp(SIG2,t2,locaisspSIG2{tam,1},locaisdiastSIG2{tam,1},0.2); % increases
    % precision concerning the time features
    
    ppg_PAT=ppg80_time; %t2(locaisdiastSIG2) ou ppg80_time
    [patt,locsr,locsppg]=PATcalculator(ecg_PAT,ppg_PAT); % returns patt (in deltasec) and locs (in sec)
    
    TEMPPAT{tam,1}=cat(1,TEMPPAT{tam,1},patt);
    PAT_time{tam,1}=cat(1,PAT_time{tam,1},locsr);
    
end

%% 1) 1st LEVEL of features (inspired in SPI, a RR and PPGA based algorithm)

% LEAVE ONE OUT CROSS VALIDATION
level=1;

if level==1
    
    fs=n_fs;
    figure('Name','1st level')
    for l=1:size(ppga_group,1)
        
        ppga_train=cell(tot_tam,1);
        rr_train=cell(tot_tam,1);
        
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
%%
