% sem erros ainda
clear all;

addpath(genpath('C:\Users\danis\Desktop\fff\tese\matlab'));

%signal=1;
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
% restrict=[begin 3005; 5890 5905; 10263 10267; 10606 10608; 9721 9723; 9853 9855; 9637 9639; 9643 9644; 9661 9666; 3040 3070; 3210 3281; 3532 3550; 3671 3700; 7003 7004; 7125 7127; 7681 7683; 7933 7935; 8400 8430; 8490 8517; 8631 8632; 8694 8771; 9041 9043; 9132 9133; 9156 9157; 9190 9191; 9322 9323; 9499 9500; 9516 9519; 9746 9748; 9887 9888; 9908 9911; 10587 10593; 10599 10601; 10625 10637; 10735 10743; 10995 11025; 11042 11044; 11120 11130; 11355 11360; 11207 11228; 10860 10875; 6333 6335; 6238 6241; 6120 6124; 5652 5656; 5580 5585; 5205 5280; 4612 4615; 4450 4454; 4654 4657; 4733 4737; 6817 6837; 6317 6319; 6678 6680; 6455 6460];
% restrict=[begin 2140; 2151 2160; 2174 2180; 2194 2199; 2232 2234; 2242 2250; 2284 2290; 2298 2322; 2410 2423; 2436 2480; 2581 2585; 3187 3189; 3242 3244; 3344 3401; 3536 3591; 3613 3690; 3824 3840; 3939 3943; 4171 4172; 4783 4787; 4806 4808; 4988 4990; 5146 5181; 5940 5943; 6320 6328; 6405 6408; 6456 6467; 6820 6823; 6865 6868; 6942 7035; 7443 7448; 7698 7704; 8464 8496; 8693 8791; 9527 9555; 10540 10580; 10650 10680; 10800 10810; 10920 10950; 11034 11037; 11553 11557; 11654 11660; 12050 12055; 12190 12195; 12337 12375; 12520 12524; 12695 12742; 12867 12890; 12910 12920; 12948 12985; 13057 13060; 13110 13440; 13478 13715; 13750 13762]; %#3
restrict=[2579 2600; 9096 9450; 7835 7840; 1245 1430; 1693 1695; 1929 1946; 2437 2444; 2677 2686; 3330 3333; 3700 4603; 4940 4945; 4989 4992; 5396 5400; 5848 5998; 6666 6669; 6927 7045; 8362 8397; 8505 8640; 8865 8943; 9538 9540; 9779 9781; 10455 10530; 10795 10895; 11010 11045; 11275 11320; 11390 11400; 11434 ex_end]; %; #5

% ppg_train=measure(2).data(find(measure(2).time()>begin & measure(2).time()<fim))';
% ex_clean_ppg_DC=measure(3).data(find(measure(3).time()>begin & measure(3).time()<fim))'; %DC+AC
ex_clean_abp=measure(2).data(find(measure(2).time()>ex_begin & measure(2).time()<ex_end))';
t1=measure(2).time(find(measure(2).time()>ex_begin & measure(2).time()<ex_end));
ex_clean_ppg_DC=measure(3).data(find(measure(3).time()>ex_begin & measure(3).time()<ex_end))'; %DC+AC
t2=measure(3).time(find(measure(3).time()>ex_begin & measure(3).time()<ex_end));
ex_clean_ECG=measure(1).data(find(measure(1).time()>ex_begin & measure(1).time()<ex_end))';
t3=measure(1).time(find(measure(1).time()>ex_begin & measure(1).time()<ex_end));
% %%
% PLOT OF THE 3 SIGNALS
figure
ax(1)=subplot(2,1,1);
plot(measure(2).time()./60,measure(2).data(),'g-',measure(2).time(find(measure(2).time()>begin & measure(2).time()<fim))./60,measure(2).data(find(measure(2).time()>begin & measure(2).time()<fim)),'b-',t1./60,ex_clean_abp,'r-');
ylabel('mmHg');
xlabel('time (min)');
ax(2)=subplot(2,1,2);
% plot(measure(1).time(),measure(1).data(),'g-',measure(1).time(find(measure(1).time()>begin & measure(1).time()<fim)),measure(1).data(find(measure(1).time()>begin & measure(1).time()<fim)),'b-',t3,ex_clean_ECG,'r-');
% ylabel('Amplitude (mV)');
% xlabel('time (sec)');
% legend('full ECG','surgery ECG','analyzed ECG');
% ax(3)=subplot(3,1,3);
plot(measure(3).time()./60,measure(3).data(),'g-',measure(3).time(find(measure(3).time()>begin & measure(3).time()<fim))./60,measure(3).data(find(measure(3).time()>begin & measure(3).time()<fim)),'b-',t2./60,ex_clean_ppg_DC,'r-');
ylabel('amplitude (mNp)');
xlabel('time (min)');
% % ax(3)=subplot(3,1,3);
% % plot(measure(7).time(),measure(7).data(),'b-');
% % legend('Surgery period');
linkaxes(ax,'x');
%%
% DC removal and y-axis flip
fs=125;
wc1=0.5;
wc2=10;
ordem=5;
ppg_train_new=process(ordem,fs,wc1,wc2,ex_clean_abp,1);
abp_2=ppg_train_new;

fs=125;
wc1=0.5;
wc2=1;
ordem=5;
ppg_train_new=process(ordem,fs,wc1,wc2,ex_clean_ppg_DC,2);
ppg_train=-ppg_train_new;

fs=250;
wc1=4;
wc2=20;
ordem=5;
ppg_train_new=process(ordem,fs,wc1,wc2,ex_clean_ECG,3);
ECG_2=ppg_train_new;

%
% figure
% ax(1)=subplot (2,1,1);
% plot(ppg_train)
% legend('PPG AC component');
% ylabel('Np')
% xlabel('Time(sec)')
% ax(2)=subplot (2,1,2);
% plot(ex_clean_ppg_DC)
% legend('PPG DC component');
% ylabel('Np')
% xlabel('Time(sec)')
% linkaxes(ax,'x');

% --------------------------------------------------- PPG and BP features --------------------------------------------------
tam=1;
n_fs=8;
rr_group=cell(tam,1);
sbp_fin=cell(tam,1);
SIG2_amp_fin=cell(tam,1);
DC_ppga_fin=cell(tam,1);
locaisDC_ppga_fin=cell(tam,1);
DC_max=cell(tam,1);
locaisDC_min=cell(tam,1);
DC_min=cell(tam,1);
diastolic_fin=cell(tam,1);
rr_fin=cell(tam,1);
locaisspSIG1=cell(tam,1);
locaisspSIG2=cell(tam,1);
locaisdiastSIG2=cell(tam,1);
b_a=cell(tam,1);
c_val=cell(tam,1);
b_val=cell(tam,1);
a_val=cell(tam,1);
e_val=cell(tam,1);
RI_temp=cell(tam,1);
SI_temp=cell(tam,1);
es=cell(tam,1);
as=cell(tam,1);
bs=cell(tam,1);
cs=cell(tam,1);
dicnotch_temp=cell(tam,1);
locaisdicnotch=cell(tam,1);
mm=1;
diast_temp=cell(tam,1);
locaisdiast=cell(tam,1);
locaisdiastSIG1=cell(tam,1);
TEMPPAT=cell(tam,1);
PAT_time=cell(tam,1);
TEMPPAT2=cell(tam,1);
PAT_time2=cell(tam,1);

CORR_SHOW=0;

[picoppg,localppg]=findpeaks(ppg_train(:,mm),'MinPeakDistance',60,'MinPeakWidth',15); % distance of 60 in order to avoid capturing small peaks
[picoabp,localabp]=findpeaks(ex_clean_abp(:,mm),'MinPeakDistance',60,'MinPeakWidth',15,'MinPeakHeight',70); % distance of 60 in order to avoid capturing small peaks

if CORR_SHOW==1
    figure;
end

SIG1=ex_clean_abp; %ex_clean_ECG
SIG2=ppg_train;

% -------------------------------------------------------------- BP features
for i=2:length(localabp)-2
    interval1=localabp(i-1):localabp(i);
    [minimo1,local_min1]= min(SIG1(interval1,mm));         %min global
    
    locaisspSIG1{mm,1}=cat(1,locaisspSIG1{mm,1},localabp(i));
    locaisdiastSIG1{mm,1}=cat(1,locaisdiastSIG1{mm,1},interval1(1)+local_min1);
end

%------------------------------------------------------diastolic BP peak
diastolic_fin{mm,:}=SIG1(locaisdiastSIG1{mm,1},mm);

%------------------------------------------------------systolic BP peak
sbp_fin{mm,:}=SIG1(locaisspSIG1{mm,1},mm);

%--------------------------------------------- Mean arterial pressure (MAP)
MAP_temp{mm,1}=(SIG1(locaisspSIG1{mm,1},mm)+2.*SIG1(locaisdiastSIG1{mm,1},mm))./3;

%--------------------------------------------- Pulse pressure (PP)
PP_fin{mm,1}=sbp_fin{mm,1}-diastolic_fin{mm,1};


%---------------------------------------------------------- PPG features
for i=2:length(localppg)-1
    interval2=localppg(i-1):localppg(i);
    
    %---------------------------------------------------------- AC amplitude
    [minimo2,local_min2]= min(SIG2(interval2,mm));
    amp_SIG2=SIG2(localppg(i),mm)-minimo2;
    SIG2_amp_fin{mm,1}=cat(1,SIG2_amp_fin{mm,1},amp_SIG2);
    
    %---------------------------------------------------------- PPG Diastolic and systolic locations
    locaisdiastSIG2{mm,1}=cat(1,locaisdiastSIG2{mm,1},local_min2+interval2(1));
    locaisspSIG2{mm,1}=cat(1,locaisspSIG2{mm,1},localppg(i));
    
    % ------------------------------------------------------------ AMP DC PPG
    [picotemp,localtemp]=findpeaks(ex_clean_ppg_DC(interval2,1),'MinPeakDistance',50,'MinPeakWidth',15); % distance of 50 in order to avoid capturing small peaks
    if isempty(localtemp)
        [picotemp,localtemp]=max(ex_clean_ppg_DC(interval2,1));
    end
    [minimo2,local_min2]= min(ex_clean_ppg_DC(interval2,1));         %min global
    amp_ppg_DC=(picotemp(1)+minimo2(1))/2;
    
    DC_max{mm,1}=cat(1,DC_max{mm,1},picotemp);
    DC_ppga_fin{mm,1}=cat(1,DC_ppga_fin{mm,1},amp_ppg_DC);
    locaisDC_ppga_fin{mm,1}=cat(1,locaisDC_ppga_fin{mm,1},interval2(1)+localtemp);
    
    DC_min{mm,1}=cat(1,DC_min{mm,1},minimo2);
    locaisDC_min{mm,1}=cat(1,locaisDC_min{mm,1},interval2(1)+local_min2);
    
    
    %------------------------------------------------------------------ 2ND LEVEL OF PPG FEATURES
    interval=localppg(i)-50:localppg(i)+50;
    
    if CORR_SHOW==1
        temp_inter=localppg(i)-50:localppg(i)+50;
        plot(temp_inter,ppg_train(temp_inter,mm),'b-');
        hold on;
    end
    
    %-------------------------------------------------------- a,b,c and e/DN of PPG
    second_ppg_spe=(diff(diff(ppg_train(interval,mm))));
    wc2=17;
    fs=125;
    ordem=5;
    second_ppg_spe=process(ordem,fs,1,wc2,second_ppg_spe,1);
    
    [picosec,localsec]=findpeaks(second_ppg_spe,'MinPeakWidth',3,'SortStr','descend');
    if (length(picosec)<1)
        [picosec,localsec]=findpeaks(second_ppg_spe,'MinPeakWidth',2,'SortStr','descend');
        if (length(picosec)<1) %if clause for the manual removed areas (yet to be removed)
            ai=interval(1)+aip;
            picosec=a_val{mm,1}(end);
            localsec=20; %random number just for absent regions
        else
            ai=interval(1)+localsec(1);
            aip=localsec(1);
        end
    else
        ai=interval(1)+localsec(1);
        aip=localsec(1);
    end
    as{mm,1}=cat(1,as{mm,1},ai);
    a_val{mm,1}=cat(1,a_val{mm,1},picosec(1));
    
    [picosec1,localsec1]=findpeaks(-second_ppg_spe(localsec(1):end),'MinPeakWidth',3); % distance of 50 in order to avoid capturing small peaks
    if (length(picosec1)<1)
        [picosec1,localsec1]=findpeaks(-second_ppg_spe(localsec(1):end),'MinPeakWidth',1);
        if (length(picosec1)<1) %if clause for the manual removed areas (yet to be removed)
            bi=ai+bip;
            picosec1=-b_val{mm,1}(end);
            localsec1=-5; %noisy regions purposes (no implications in the results)
        else
            bi=ai+localsec1(1)-1;
            bip=localsec1(1);
        end
    else
        bi=ai+localsec1(1)-1;
        bip=localsec1(1);
    end
    bs{mm,1}=cat(1,bs{mm,1},bi);
    b_val{mm,1}=cat(1,b_val{mm,1},-picosec1(1));
    
    int=localsec(1)+localsec1(1);
    %----------------------------------------------------------diastolic peak in elgendi article
    interval=localppg(i)+10:localppg(i)+50; % INTERVAL TO FIND diastolic POSITION ("local2(i)+10" is to help finding the nearest point to zero (which will be the diastolic peak) )
    
    first_ppg=diff(ppg_train(interval,mm));
    wc2=17;
    fs=125;
    ordem=5;
    first_ppg=process(ordem,fs,1,wc2,first_ppg,1);
    
    second_ppg_dias=(diff(diff(ppg_train(interval,mm))));
    wc2=17;
    fs=125;
    ordem=5;
    second_ppg_dias=process(ordem,fs,1,wc2,second_ppg_dias,1);
    [~,pos] = findpeaks(-second_ppg_dias);
    
    loc1=find(first_ppg>0);
    compara=any(diff(loc1)>2); %when there is a gap in this vector, it means that f' went from positive to negative
    [vv,loc2]=min(abs(first_ppg));
    
    if isempty(pos)==0
        local_diast=interval(1)+pos(1);
        diast=ppg_train(local_diast,mm);
        
    elseif ( compara==1 ) %when 1st derivative goes from positive to negative values
        pos=find(diff(loc1)>2);
        local_diast=interval(1)+loc1(pos);
        diast=ppg_train(local_diast,mm);
        
        if CORR_SHOW==1
            plot(local_diast,diast,'k*');
            hold on;
        end
    elseif ( compara==0 & isempty(loc2)==0 ) %when 1st derivative is the closest to the zero
        
        local_diast=interval(1)+loc2(1);
        diast=ppg_train(local_diast,mm);
        
        if CORR_SHOW==1
            plot(local_diast,diast,'r*');
            hold on;
        end
    end
    
    diast_temp{mm,1}=cat(1,diast_temp{mm,1},diast(1));
    locaisdiast{mm,1}=cat(1,locaisdiast{mm,1},local_diast(1));
    %--------------------------------------------------------------c,e/DICROTIC NOTCH (same position as "e" feature in elgendi)
    int2=local_diast(1)-interval(1)+(interval(1)-(localppg(i)-50))-1;%detecao do a,b,e ate posiçao do diastolico de elgendi
    
    if int2-int<4 %noisy regions purposes (no implications in the results)
        int=int-50;
    end
    
    [picosec2,localsec2]=findpeaks(second_ppg_spe(int:int2),'MinPeakWidth',3);
    if(length(picosec2)<2)
        [picosec2,localsec2]=findpeaks(second_ppg_spe(int:int2),'MinPeakWidth',1);
        if (length(picosec2)==1) %if clause for the manual removed areas (yet to be removed)
            ci=bi+localsec2(1);
            ei=bi+localsec2(1);
            e_val{mm,1}=cat(1,e_val{mm,1},picosec2(1));
        elseif (isempty(picosec2))
            ci=interval(1)+int;
            ei=interval(1)+int;
            e_val{mm,1}=cat(1,e_val{mm,1},e_val{mm,1}(end));
            picosec2(1)=c_val{mm,1}(end);
        else
            ci=bi+localsec2(1);
            ei=bi+localsec2(2);
            e_val{mm,1}=cat(1,e_val{mm,1},picosec2(2));
        end
    else
        ci=bi+localsec2(1);
        ei=bi+localsec2(2);
        e_val{mm,1}=cat(1,e_val{mm,1},picosec2(2));
    end
    
    if CORR_SHOW==1
        plot(ai,a_val{mm,1}(end)*100,'y*',bi,b_val{mm,1}(end)*100,'g*');
        hold on;
    end
    
    if length(picosec2)>2
        local_dicnotch=bi+localsec2(3)+2;
        dicnotch=ppg_train(local_dicnotch,mm);
    else
        local_dicnotch=ei+2;
        dicnotch=ppg_train(local_dicnotch,mm);
    end
    
    es{mm,1}=cat(1,es{mm,1},ei);
    cs{mm,1}=cat(1,cs{mm,1},ci);
    c_val{mm,1}=cat(1,c_val{mm,1},picosec2(1));
    dicnotch_temp{mm,1}=cat(1,dicnotch_temp{mm,1},dicnotch);
    locaisdicnotch{mm,1}=cat(1,locaisdicnotch{mm,1},local_dicnotch);
    
    
end

%----------------------------------------------------------DN
dicnotch_fin{mm,1}=dicnotch_temp{mm,1}-ppg_train(locaisdiastSIG2{mm,1},mm);

%----------------------------------------------------------RI
% RI_temp{mm,1}=(diast_temp{mm,1}-ppg_train(locaisdiastSIG2{mm,1},mm))./SIG2_amp_fin{mm,1};
RI_temp{mm,1}=(diast_temp{mm,1}-ppg_train(locaisdiastSIG2{mm,1},mm))./SIG2_amp_fin{mm,1};

%----------------------------------------------------------SI
% SI_temp{mm,1}=abs(b_val{mm,1})./a_val{mm,1};
SI_temp{mm,1}=abs(b_val{mm,1})./a_val{mm,1};

%---------------------------------------------- Pulse transit time (PTT)
% PTT{j,1}=locaisppg{j,1}-locaissbp{j,1};

%---------------------------------------------- DELTA T
% DT{mm,1}=abs(t2(locaisspSIG2{mm,1})-t2(locaisdiast{mm,1}));

% Normalized volume
AC_DC=SIG2_amp_fin{mm,1}./DC_ppga_fin{mm,1};

%---------------------------------------------------------- ECG features

%---------------------------------------------- RR series (RR)
%Pan tompkins (R peak detection)
[b,c,ECGrpeak_pos,d,e] = pan_tompkins_qrs(ex_clean_ECG,250,0);
%
RR_fin=diff(ECGrpeak_pos)./250;
locsrr=ECGrpeak_pos(2:end);
HBItemp=diff(locaisspSIG2{mm,1})./125;
[new_HBI,newlocaishbi]=resample(HBItemp,t2(locaisspSIG2{mm,1}(2:end)),1);
locshbi=locaisspSIG2{mm,1}(2:end);
% ---------------------------------------------- Pulse Arrival time (PAT)
%between r peak in ECG and 20% amplitude of PPG systolic
ecg_PAT=t3(ECGrpeak_pos);
[ppg80_values,ppg80_time] = ppgp(SIG2,t2,locaisspSIG2{mm,1},locaisdiastSIG2{mm,1},0.2); % increases
% precision concerning the time features

ppg_PAT=ppg80_time; %t2(locaisdiastSIG2) ou ppg80_time
[patt,locsr,locsppg]=PATcalculator(ecg_PAT,ppg_PAT); % returns patt (in deltasec) and locs (in sec)

% --------------------------------------------- PAT 2.0
[pat2,newlocsr]=resample(patt,locsr,1);

fs=1;
wc1=0.00053;
wc2=0.004;
ordem=5;
newpatt=process(ordem,fs,wc1,wc2,pat2,3);

int5=[t1(1:5*125*60:end) t1(end)];
% exppos=find(t1>11434);
% int5=[t1(1:5*125*60:exppos(1))];

keeppos= [];
for k=1:length(int5)
    
    [val,pos]=min(abs(t1(locaisspSIG1{mm,1})-int5(k)));
    keeppos=[keeppos pos];
    
end
keeppos2= [];
for k=1:length(int5)
    
    [val,pos]=min(abs(locsr-int5(k)));
    keeppos2=[keeppos2 pos];
    
end

interppb=interp1(t1(locaisspSIG1{mm,1}(keeppos)),sbp_fin{mm,:}(keeppos),t1(locaisspSIG1{mm,1}(1)):t1(locaisspSIG1{mm,1}(end)));
interpatb=interp1(locsr(keeppos2),patt(keeppos2),round(locsr(1)):round(locsr(end)));

TEMPPAT=cell(tam,1);
PAT_time=cell(tam,1);
TEMPPAT2=cell(tam,1);
PAT_time2=cell(tam,1);

gama=0.016; %ou 0.018 EXPERIMENTAR GAMAS
sbp2=interppb-2*newpatt(2:end-1)/(gama*interpatb(2:end-1)); %pat 3

TEMPPAT{mm,1}=cat(1,TEMPPAT{mm,1},patt); %patt
PAT_time{mm,1}=cat(1,PAT_time{mm,1},locsr); %locsr

TEMPPAT2{mm,1}=cat(1,TEMPPAT2{mm,1},sbp2); %patt
PAT_time2{mm,1}=cat(1,PAT_time2{mm,1},newlocsr(2:end-1)); %locsr
% PAT_time2{mm,1}=cat(1,PAT_time2{mm,1},newlocsr(1:10224)); %locsr

% %% PPG waveform age test (might help to distinguish waveforms and prove why patient 3 has worse results)
%
% %younger people have bigger DT (in msec)
% %ex: 29y:1346ms, 60y:147ms
%
% me_ba=mean(SI_temp{mm,1})
%
%b/a
%1=0.8307
%3=0.3529
%5=0.7479

%DT
%1=214ms
%3=285ms
%5=183ms

%% PLOT OF ALL PEAKS
second_ppg=(diff(diff(ppg_train)));

wc2=17;
fs=125;
ordem=5;
second_ppg=process(ordem,fs,1,wc2,second_ppg,1);

first_ppg=(diff(ppg_train));

wc2=17;
fs=125;
ordem=5;
first_ppg=process(ordem,fs,1,wc2,first_ppg,1);

num=1;
figure;
ax(1)=subplot (4,1,1);
plot(t2(locaisDC_min{num,1}),DC_min{num,1},'ko',t2(locaisDC_ppga_fin{num,1}),DC_max{mm,1},'go',t2,ex_clean_ppg_DC,'r-');
ylabel('amplitude')
xlabel('Time(sec)')
legend('Diastolic peaks', 'Systolic peaks', 'PPG DC')

ax(2)=subplot (4,1,2);
plot(t2(locaisdiast{mm,1}),diast_temp{mm,1},'bo',t2(3:end),second_ppg*100,'r-',t2(locaisdicnotch{mm,1}),dicnotch_temp{mm,1},'co',t2(as{mm,1}),a_val{mm,1}*100,'g*',t2(bs{mm,1}),b_val{mm,1}*100,'b*',ppg80_time,ppg80_values,'ro',t2(locaisspSIG2{num,1}),SIG2(locaisspSIG2{num,1},num),'ko',t2,SIG2,'g-');
ylabel('amplitude')
xlabel('Time(sec)')
legend('"diastolic" (elgendi)','2nd*100','DN','a','b','20% PPG AC amplitude', 'Systolic peaks', 'PPG AC')

ax(3)=subplot (4,1,3);
plot(t1(locaisdiastSIG1{num,1}),SIG1(locaisdiastSIG1{num,1},num),'m+',t1(locaisspSIG1{num,1}),SIG1(locaisspSIG1{num,1},num),'g*',t1,SIG1,'b-');
ylabel('mmHg')
xlabel('Time(sec)')
legend('diastolic peaks','systolic peaks','BP')

ax(4)=subplot (4,1,4);
plot(t3(ECGrpeak_pos),ex_clean_ECG(ECGrpeak_pos),'m+',t3,ex_clean_ECG,'b-');%,PAT_time{mm,1},TEMPPAT{mm,1},'ro');
ylabel('mV')
xlabel('Time(sec)')
legend('systolic peaks','ECG');%,'PAT')
linkaxes(ax,'x');

%% CONVOLUTION AND DERIVATIVE FILTERING (BP vs PPG)
% ---------------------------------------------------------------------------------------PARAMETERS
convsize=10;
shift=5;

%the following are not that important because it is just changing the
%correlation window, it wont have direct influence on the calculation
%itself, just on the info gathering

correl='Spearman';  %Pearson

CORR_SHOW=0;
%-------------------------------------------------------------------------------------------------------------

amp_AC=cell(tam,1);
amp_AC_DC=cell(tam,1);
amp_DC=cell(tam,1);
MAP=cell(tam,1);
PP=cell(tam,1);
PAT=cell(tam,1);
PAT2=cell(tam,1);
locsPP=cell(tam,1);

MAP_conv=cell(tam,1);
PP_conv=cell(tam,1);
SBP_conv=cell(tam,1);
DBP_conv=cell(tam,1);
RI_conv=cell(tam,1);
DN_conv=cell(tam,1);
SI_conv=cell(tam,1);

locsplotabp=cell(tam,1);
locsplotecg=cell(tam,1);
locsplotppg=cell(tam,1);
HR=cell(tam,1);
HR_conv=cell(tam,1);
RR=cell(tam,1);
HBI_conv=cell(tam,1);

gob=[];
lgob=[];
%10996 11514;
restrict=sort(restrict,1);
bad='false';

for k=ex_begin:shift:ex_end-convsize %experimentar ate tamanho do locaisspSIG2
    
    if ( any( (k>(restrict(:,1)-convsize)) & (k<(restrict(:,2))) ) )
        bad='true';
        temp_loc=find(t2>k);
        gob=[gob 0];
        lgob=[lgob temp_loc(1)];
    else
        gob=[gob 1];
        %PPG
        ppg_time=find(t2>k & t2<k+convsize);
        lgob=[lgob ppg_time(1)];
        %         ppg_time=t2>k & t2<k+convshift;
        %         locs=find(locaisdiastSIG2{j,1}>ppg_time(1) & locaisdiastSIG2{j,1}<ppg_time(end)); % picos cujas amplitudes
        %         correspondem a percentagens da amplitude maxima
        locs=find(locshbi>ppg_time(1) & locshbi<ppg_time(end));
        amp_AC{mm,1}=cat(1,amp_AC{mm,1},mean(SIG2_amp_fin{mm,1}(locs(1):locs(end)),'omitnan'));
        amp_AC_DC{mm,1}=cat(1,amp_AC_DC{mm,1},mean(AC_DC(locs(1):locs(end)),'omitnan'));
        amp_DC{mm,1}=cat(1,amp_DC{mm,1},mean(DC_ppga_fin{mm,1}(locs(1):locs(end)),'omitnan'));
        RI_conv{mm,1}=cat(1,RI_conv{mm,1},mean(RI_temp{mm,1}(locs(1):locs(end)),'omitnan') ); % test PP values
        DN_conv{mm,1}=cat(1,DN_conv{mm,1},mean(dicnotch_fin{mm,1}(locs(1):locs(end)),'omitnan') ); % test PP values
        SI_conv{mm,1}=cat(1,SI_conv{mm,1},mean(SI_temp{mm,1}(locs(1):locs(end)),'omitnan') ); % test PP values
        HBI_conv{mm,1}=cat(1,HBI_conv{mm,1},mean(HBItemp(locs(1):locs(end)),'omitnan') ); % test PP values
        
        locsplotppg{mm,1}=cat(1,locsplotppg{mm,1},ppg_time(1));
        
        if CORR_SHOW==1
            figure;
            ax(1)=subplot (1,3,1);
            plot(locaisspSIG2{mm,1}(locs(1)):locaisspSIG2{mm,1}(locs(end)),SIG2(locaisspSIG2{mm,1}(locs(1)):locaisspSIG2{mm,1}(locs(end))),'b-',locaisspSIG2{mm,1}(locs),SIG2(locaisspSIG2{mm,1}(locs)),'ro',locaisspSIG2{mm,1}(round((locs(1)+locs(end))/2)),amp_AC{mm,1}(end),'r+',locaisspSIG2{mm,1}(round((locs(1)+locs(end))/2)),amp_DC{mm,1}(end),'k+');
            legend('PPG','systolic PPG','ACPPG','DCPPG');
        end
        
        %BP
        abp_time=find(t1>k & t1<k+convsize);
        locs=find(locaisspSIG1{mm,1}>abp_time(1) & locaisspSIG1{mm,1}<abp_time(end));
        locsdiast=find(locaisdiastSIG1{mm,1}>abp_time(1) & locaisdiastSIG1{mm,1}<abp_time(end));
        
        %                 MAP{j,1}=cat(1,MAP{j,1},mean(DER_MAP{j,1}(locs(1):locs(end)-1),'omitnan') );% test MAP variations
        MAP_conv{mm,1}=cat(1,MAP_conv{mm,1},mean(MAP_temp{mm,1}(locs(1):locs(end)),'omitnan') ); % test MAP values
        %                 PP{j,1}=cat(1,PP{j,1},mean(DER_PP{j,1}(locs(1):locs(end)-1),'omitnan')); % test PP variations
        PP_conv{mm,1}=cat(1,PP_conv{mm,1},mean(PP_fin{mm,1}(locs(1):locs(end)),'omitnan') ); % test PP values
        %                 SBP{j,1}=cat(1,SBP{j,1},mean(DER_SBP{j,1}(locs(1):locs(end)-1),'omitnan') );% test MAP variations
        SBP_conv{mm,1}=cat(1,SBP_conv{mm,1},mean(sbp_fin{mm,1}(locs(1):locs(end)),'omitnan') ); % test MAP values
        %                 DBP{j,1}=cat(1,DBP{j,1},mean(DER_DBP{j,1}(locsdiast(1):locsdiast(end)-1),'omitnan')); % test PP variations
        DBP_conv{mm,1}=cat(1,DBP_conv{mm,1},mean(diastolic_fin{mm,1}(locs(1):locs(end)),'omitnan') ); % test PP values
        
        locsplotabp{mm,1}=cat(1,locsplotabp{mm,1},abp_time(1));
        
        if CORR_SHOW==1
            ax(2)=subplot (1,3,2);
            plot(locaisspSIG1{mm,1}(locs(1)):locaisspSIG1{mm,1}(locs(end)),SIG1(locaisspSIG1{mm,1}(locs(1)):locaisspSIG1{mm,1}(locs(end))),'b-',locaisspSIG1{mm,1}(locs),SIG1(locaisspSIG1{mm,1}(locs)),'ro'); %,locaisspSIG1{j,1}(round((locs(1)+locs(end))/2)),MAP_conv{j,1}(end),'r+',locaisspSIG1{j,1}(round((locs(1)+locs(end))/2)),PP_conv{j,1}(end),'k+');
            legend('BP','systolic BP');%,'MAP','PP');
        end
        
        %ECG
        ecg_time=find(t3>k & t3<k+convsize);
        locs=find(PAT_time{mm,1}>k & PAT_time{mm,1}<k+convsize);
        PAT{mm,1}=cat(1,PAT{mm,1},mean(TEMPPAT{mm,1}(1,locs(1):locs(end)),'omitnan') );
        
        locs=find(PAT_time2{mm,1}>k & PAT_time2{mm,1}<k+convsize);
        PAT2{mm,1}=cat(1,PAT2{mm,1},mean(TEMPPAT2{mm,1}(1,locs(1):locs(end)),'omitnan') );
        
        %         locs=find(locsrr>ecg_time(1) & locsrr<ecg_time(end));
        %         RR{mm,1}=cat(1,RR{mm,1},mean(RR_fin(locs(1):locs(end)),'omitnan') );
        %         HR{mm,1}=cat(1,HR{mm,1},60*numel(locs)/convsize); % variacao tambem??????????????????
        locsplotecg{mm,1}=cat(1,locsplotecg{mm,1},ecg_time(1));
        
        if CORR_SHOW==1
            ax(3)=subplot (1,3,3);
            plot(ECGrpeak_pos(locs(1)):ECGrpeak_pos(locs(end)),ex_clean_ECG(ECGrpeak_pos(locs(1)):ECGrpeak_pos(locs(end))),'b-',ECGrpeak_pos(locs),ex_clean_ECG(ECGrpeak_pos(locs)),'ro'); %,exemp_ECG(round((locs(1)+locs(end))/2)),PAT{j,1}(end),'r+');
            legend('ECG','R peak');%,'PAT');
            linkaxes(ax,'x');
            pause
        end
    end
end

% POS DERIVATIVE
MAP{mm,1}=(MAP_conv{mm,1});%./diff(locaisspSIG1{j,1});
PP{mm,1}=(PP_conv{mm,1});%./diff(locaisspSIG1{j,1});
SBP{mm,1}=(SBP_conv{mm,1});%./diff(locaisspSIG1{j,1});
DBP{mm,1}=(DBP_conv{mm,1});%./diff(locaisspSIG1{j,1});
RI{mm,1}=RI_conv{mm,1};
DN{mm,1}=DN_conv{mm,1};
SI=SI_conv{mm,1};
HBI=HBI_conv{mm,1};
% %% normalization zscore
% 
% RI{mm,1}=zscore(RI_conv{mm,1});
% DN{mm,1}=zscore(DN_conv{mm,1});
% SI=zscore(SI_conv{mm,1});
% HBI=zscore(HBI_conv{mm,1});
% PAT{mm,1}= zscore(PAT{mm,1});
% PAT2{mm,1}=zscore(PAT2{mm,1});
% amp_AC{1,1}=zscore(amp_AC{1,1});
% amp_AC_DC{1,1}=zscore(amp_AC_DC{1,1});
% amp_DC{1,1}=zscore(amp_DC{1,1});
% 
% corrpp=[corr(amp_AC_DC{1,1},PP{mm,1},'type','Spearman') corr(amp_AC{1,1},PP{mm,1},'type','Spearman') corr(amp_DC{1,1},PP{mm,1},'type','Spearman') corr(PAT{1,1},PP{mm,1},'type','Spearman') corr(PAT2{1,1},PP{mm,1},'type','Spearman') corr(HBI,PP{mm,1},'type','Spearman') corr(RI{mm,1},PP{mm,1},'type','Spearman') corr(DN{1,1},PP{mm,1},'type','Spearman') corr(SI,PP{mm,1},'type','Spearman')]
% corrsbp=[corr(amp_AC_DC{1,1},SBP{mm,1},'type','Spearman') corr(amp_AC{1,1},SBP{mm,1},'type','Spearman') corr(amp_DC{1,1},SBP{mm,1},'type','Spearman') corr(PAT{1,1},SBP{mm,1},'type','Spearman') corr(PAT2{1,1},SBP{mm,1},'type','Spearman') corr(HBI,SBP{mm,1},'type','Spearman') corr(RI{mm,1},SBP{mm,1},'type','Spearman') corr(DN{1,1},SBP{mm,1},'type','Spearman') corr(SI,SBP{mm,1},'type','Spearman')]
% corrdbp=[corr(amp_AC_DC{1,1},DBP{mm,1},'type','Spearman') corr(amp_AC{1,1},DBP{mm,1},'type','Spearman') corr(amp_DC{1,1},DBP{mm,1},'type','Spearman') corr(PAT{1,1},DBP{mm,1},'type','Spearman') corr(PAT2{1,1},DBP{mm,1},'type','Spearman') corr(HBI,DBP{mm,1},'type','Spearman') corr(RI{mm,1},DBP{mm,1},'type','Spearman') corr(DN{1,1},DBP{mm,1},'type','Spearman') corr(SI,DBP{mm,1},'type','Spearman')]
% corrmap=[corr(amp_AC_DC{1,1},MAP{mm,1},'type','Spearman') corr(amp_AC{1,1},MAP{mm,1},'type','Spearman') corr(amp_DC{1,1},MAP{mm,1},'type','Spearman') corr(PAT{1,1},MAP{mm,1},'type','Spearman') corr(PAT2{1,1},MAP{mm,1},'type','Spearman') corr(HBI,MAP{mm,1},'type','Spearman') corr(RI{mm,1},MAP{mm,1},'type','Spearman') corr(DN{1,1},MAP{mm,1},'type','Spearman') corr(SI,MAP{mm,1},'type','Spearman')]
% 
% %%
% bxpp=[];
% bxmap=[];
% bxdbp=[];
% bxsbp=[];
% iqrpp=[];
% iqrsbp=[];
% iqrdbp=[];
% iqrmap=[];
% 
% mdnpp=[];
% mdnsbp=[];
% mdndbp=[];
% mdnmap=[];
% 
% % correlaçoes
% 
% corrpp=[corr(amp_AC_DC{1,1},PP{mm,1},'type','Spearman') corr(amp_AC{1,1},PP{mm,1},'type','Spearman') corr(amp_DC{1,1},PP{mm,1},'type','Spearman') corr(PAT{1,1},PP{mm,1},'type','Spearman') corr(PAT2{1,1},PP{mm,1},'type','Spearman') corr(HBI,PP{mm,1},'type','Spearman') corr(RI{mm,1},PP{mm,1},'type','Spearman') corr(DN{1,1},PP{mm,1},'type','Spearman') corr(SI,PP{mm,1},'type','Spearman')]
% corrsbp=[corr(amp_AC_DC{1,1},SBP{mm,1},'type','Spearman') corr(amp_AC{1,1},SBP{mm,1},'type','Spearman') corr(amp_DC{1,1},SBP{mm,1},'type','Spearman') corr(PAT{1,1},SBP{mm,1},'type','Spearman') corr(PAT2{1,1},SBP{mm,1},'type','Spearman') corr(HBI,SBP{mm,1},'type','Spearman') corr(RI{mm,1},SBP{mm,1},'type','Spearman') corr(DN{1,1},SBP{mm,1},'type','Spearman') corr(SI,SBP{mm,1},'type','Spearman')]
% corrdbp=[corr(amp_AC_DC{1,1},DBP{mm,1},'type','Spearman') corr(amp_AC{1,1},DBP{mm,1},'type','Spearman') corr(amp_DC{1,1},DBP{mm,1},'type','Spearman') corr(PAT{1,1},DBP{mm,1},'type','Spearman') corr(PAT2{1,1},DBP{mm,1},'type','Spearman') corr(HBI,DBP{mm,1},'type','Spearman') corr(RI{mm,1},DBP{mm,1},'type','Spearman') corr(DN{1,1},DBP{mm,1},'type','Spearman') corr(SI,DBP{mm,1},'type','Spearman')]
% corrmap=[corr(amp_AC_DC{1,1},MAP{mm,1},'type','Spearman') corr(amp_AC{1,1},MAP{mm,1},'type','Spearman') corr(amp_DC{1,1},MAP{mm,1},'type','Spearman') corr(PAT{1,1},MAP{mm,1},'type','Spearman') corr(PAT2{1,1},MAP{mm,1},'type','Spearman') corr(HBI,MAP{mm,1},'type','Spearman') corr(RI{mm,1},MAP{mm,1},'type','Spearman') corr(DN{1,1},MAP{mm,1},'type','Spearman') corr(SI,MAP{mm,1},'type','Spearman')]
% 
%% EVALUATION OF CORRELATION!!

% PROBLEMS WITH PEARSON CORR: SENSITIVE TO OUTLIERS!!

plotar=0; % plotar evoluçao de SCC e features
bplot=0; % plotar boxplot=1

DATA=3; %1=AC_DC(+), 2=AC_ppga(+), 3=DC_ppga(-) %4=PAT(-) %5=PAT2(+) 6=RR(-) %7=RI(+) %8=DN(+) %9=SI(-)
threshold=0.7;
int1=0.7;
int2=1;
good_PP=[];
good_MAP=[];
correlation_coef=[];
Spearman_coef=[];
pos_beg=[];
pos_end=[];
lin_reg_coef=[];
b_reg_coef=[];

l=1;

%chi squared (1: rejeita dist. normal, 0: nao rejeita)
chi_PP1{l,1}=cell(tam,1);
chi_MAP1{l,1}=cell(tam,1);
chi_SBP1{l,1}=cell(tam,1);
chi_DBP1{l,1}=cell(tam,1);

%Coef var (normal units): helps to realise if good correlation is related with
%high feature variations or not
cv_PP1{l,1}=cell(tam,1);
cv_MAP1{l,1}=cell(tam,1);
cv_SBP1{l,1}=cell(tam,1);
cv_DBP1{l,1}=cell(tam,1);

%TEST
wind_size=23; %4,9,14/7,17,27/11,23,35 MEXER NESTE!!!!!!! PARA 10,20,30,60,120,300 SEC
CORR_shift=wind_size;

CORR_SHOW=0;
%---------------------------------------------------------------------
if strcmp(correl,'Spearman')
    
    if DATA==2
        
        dados='AC';
        
        test=amp_AC{mm,1};
        
        for i=1:CORR_shift:length(test)-wind_size % VARY SIZE OF WINDOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            test1=test(i:i+wind_size);
            PP1{l,1}=PP{l,1}(i:i+wind_size);
            SBP1{l,1}=SBP{l,1}(i:i+wind_size);
            DBP1{l,1}=DBP{l,1}(i:i+wind_size);
            MAP1{l,1}=MAP{l,1}(i:i+wind_size);
            
            chi_PP1{l,1}=cat(1,chi_PP1{l,1},chi2gof(PP1{l,1}));
            chi_SBP1{l,1}=cat(1,chi_SBP1{l,1},chi2gof(SBP1{l,1}));
            chi_DBP1{l,1}=cat(1,chi_DBP1{l,1},chi2gof(DBP1{l,1}));
            chi_MAP1{l,1}=cat(1,chi_MAP1{l,1},chi2gof(MAP1{l,1}));
            
            cv_PP1{l,1}=cat(1,cv_PP1{l,1},std(PP1{l,1})./mean(PP1{l,1}));
            cv_MAP1{l,1}=cat(1,cv_MAP1{l,1},std(MAP1{l,1})./mean(MAP1{l,1}));
            cv_SBP1{l,1}=cat(1,cv_SBP1{l,1},std(SBP1{l,1})./mean(SBP1{l,1}));
            cv_DBP1{l,1}=cat(1,cv_DBP1{l,1},std(DBP1{l,1})./mean(DBP1{l,1}));
            
            pos_beg=[pos_beg t2(locsplotppg{mm,1}(i))./60];
            pos_end=[pos_end t2(locsplotppg{mm,1}(i+wind_size))./60];
            
            %--------------------------------------PLOT OF THE REGION IN ANALYSIS
            if CORR_SHOW==1
                
                figure
                ax(1)=subplot (3,1,1);
                plot(t3(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),ex_clean_ECG(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),'b-',t3(ECGrpeak_pos),ex_clean_ECG(ECGrpeak_pos),'g*');
                legend('ecg','systolic peaks')
                ylabel('mV')
                xlabel('Time(sec)')
                ax(2)=subplot (3,1,2);
                plot(t2(locsplotppg{mm,1}(i:i+wind_size)),test1,'ko',t2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),SIG2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),'b-',ppg80_time,ppg80_values,'r*',t2(locaisspSIG2{num,1}),SIG2(locaisspSIG2{num,1},num),'g*');
                legend('AC PPG (AMPLITUDE)','ppg','20% diastolic peaks', 'Systolic peaks')
                ylabel('amplitude')
                xlabel('Time(sec)')
                ax(3)=subplot (3,1,3);
                plot(t1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),SIG1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),'b-',t1(locsplotabp{mm,1}(i:i+wind_size)),MAP1{l,1},'k+',t1(locsplotabp{mm,1}(i:i+wind_size)),PP1{l,1},'go',t1(locaisdiastSIG1{num,1}),SIG1(locaisdiastSIG1{num,1},num),'r*',t1(locaisspSIG1{num,1}),SIG1(locaisspSIG1{num,1},num),'g*');
                linkaxes(ax,'x');
                legend('bp' ,'MAP','PP','diastolic peaks','systolic peaks');
                ylabel('mmHg')
                xlabel('Time(sec)')
            end
            
            %ERROR IN POLYFIT?!?!!?!?!?!?
            [A]=polyfit(test1,PP1{l,1},1); %DBP
            b=A(end);
            a=A(end-1);
            
            [A1]=polyfit(test1,SBP1{l,1},1); %DBP
            b1=A1(end);
            a1=A1(end-1);
            
            [A2]=polyfit(test1,DBP1{l,1},1); %SBP
            b2=A2(end);
            a2=A2(end-1);
            
            [A4]=polyfit(test1,MAP1{l,1},1); %SBP
            b4=A4(end);
            a4=A4(end-1);
            
            R=corr(test1,PP1{l,1},'type','Spearman');
            
            R1=corr(test1,SBP1{l,1},'type','Spearman');
            
            R2=corr(test1,DBP1{l,1},'type','Spearman');
            
            R4=corr(test1,MAP1{l,1},'type','Spearman');
            
            %------------------------------------- PLOTS WITH CORRELATION COEFS
            %AND LIN_REG COEFS
            
            lin_reg_coef=cat(2,lin_reg_coef,[A(1,1);A1(1,1);A2(1,1);A4(1,1)]);
            b_reg_coef=cat(2,b_reg_coef,[A(end);A1(end);A2(end);A4(end)]);
            
            Spearman_coef=cat(2,Spearman_coef,[R;R1;R2;R4]);
            
            rele_wind_ac=cell(1,4);
            
            feat=1;
            nn_pp=find(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)<0);
            feat=2;
            nn_sbp=find(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)<0);
            feat=3;
            nn_dbp=find(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)<0);
            feat=4;
            nn_map=find(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)<0);
            
            rele_wind_ac{1,1}=nn_pp;
            rele_wind_ac{1,2}=nn_sbp;
            rele_wind_ac{1,3}=nn_dbp;
            rele_wind_ac{1,4}=nn_map;
            
            if CORR_SHOW==1
                
                figure;
                
                subplot 422
                
                plot((test1),PP1{l,1},'go');
                xlabel('AC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('PP (mmHg)')
                str1 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a, R);
                title(str1)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,PP1{l,1},txt1)
                
                subplot 424
                
                plot((test1),SBP1{l,1},'go');
                xlabel('AC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('SBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a1, R1);
                title(str)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,SBP1{l,1},txt1)
                
                subplot 426
                
                plot((test1),DBP1{l,1},'go');
                xlabel('AC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('DBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a2, R2);
                title(str)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,DBP1{l,1},txt1)
                
                subplot 428
                plot((test1),MAP1{l,1},'go');
                xlabel('AC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('MAP (mmHg)')
                str2 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a4, R4);
                title(str2)
                
                ax(1)=subplot(4,2,[1,3])
                plot(t3(locsplotecg{mm,1}(i:i+wind_size)),test1)
                ylabel('AC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                xlabel('Time(sec)')
                
                ax(2)=subplot(4,2,[5,7])
                plot(t1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),SIG1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)))
                legend('BP');
                ylabel('mmHg')
                xlabel('Time(sec)')
                linkaxes(ax,'x');
                
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,MAP1{l,1},txt1)
                
                pause
                
            end
            
            
        end
        ACe=Spearman_coef;
        
    elseif DATA==8
        dados='DN';
        
        test=DN{mm,1};%SI
        
        for i=1:CORR_shift:length(test)-wind_size % VARY SIZE OF WINDOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            test1=test(i:i+wind_size);
            
            MAP1{l,1}=MAP{l,1}(i:i+wind_size);
            PP1{l,1}=PP{l,1}(i:i+wind_size);
            SBP1{l,1}=SBP{l,1}(i:i+wind_size);
            DBP1{l,1}=DBP{l,1}(i:i+wind_size);
            
            chi_PP1{l,1}=cat(1,chi_PP1{l,1},chi2gof(PP1{l,1}));
            chi_SBP1{l,1}=cat(1,chi_SBP1{l,1},chi2gof(SBP1{l,1}));
            chi_DBP1{l,1}=cat(1,chi_DBP1{l,1},chi2gof(DBP1{l,1}));
            chi_MAP1{l,1}=cat(1,chi_MAP1{l,1},chi2gof(MAP1{l,1}));
            
            cv_PP1{l,1}=cat(1,cv_PP1{l,1},std(PP1{l,1})./mean(PP1{l,1}));
            cv_MAP1{l,1}=cat(1,cv_MAP1{l,1},std(MAP1{l,1})./mean(MAP1{l,1}));
            cv_SBP1{l,1}=cat(1,cv_SBP1{l,1},std(SBP1{l,1})./mean(SBP1{l,1}));
            cv_DBP1{l,1}=cat(1,cv_DBP1{l,1},std(DBP1{l,1})./mean(DBP1{l,1}));
            
            pos_beg=[pos_beg t2(locsplotppg{mm,1}(i))./60];
            pos_end=[pos_end t2(locsplotppg{mm,1}(i+wind_size))./60];
            
            %ERROR IN POLYFIT?!?!!?!?!?!?
            [A]=polyfit(test1,PP1{l,1},1); %DBP
            b=A(end);
            a=A(end-1);
            
            [A1]=polyfit(test1,SBP1{l,1},1); %DBP
            b1=A1(end);
            a1=A1(end-1);
            
            [A2]=polyfit(test1,DBP1{l,1},1); %SBP
            b2=A2(end);
            a2=A2(end-1);
            
            [A4]=polyfit(test1,MAP1{l,1},1); %SBP
            b4=A4(end);
            a4=A4(end-1);
            
            R=corr(test1,PP1{l,1},'type','Spearman');
            
            R1=corr(test1,SBP1{l,1},'type','Spearman');
            
            R2=corr(test1,DBP1{l,1},'type','Spearman');
            
            R4=corr(test1,MAP1{l,1},'type','Spearman');
            
            %------------------------------------- PLOTS WITH CORRELATION COEFS
            %AND LIN_REG COEFS
            
            lin_reg_coef=cat(2,lin_reg_coef,[A(1,1);A1(1,1);A2(1,1);A4(1,1)]);
            b_reg_coef=cat(2,b_reg_coef,[A(end);A1(end);A2(end);A4(end)]);
            
            Spearman_coef=cat(2,Spearman_coef,[R;R1;R2;R4]);
            
            rele_wind_dn=cell(1,4);
            
            feat=1;
            nn_pp=find(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)<0);
            feat=2;
            nn_sbp=find(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)<0);
            feat=3;
            nn_dbp=find(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)<0);
            feat=4;
            nn_map=find(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)<0);
            
            rele_wind_dn{1,1}=nn_pp;
            rele_wind_dn{1,2}=nn_sbp;
            rele_wind_dn{1,3}=nn_dbp;
            rele_wind_dn{1,4}=nn_map;
            
            
            
            %--------------------------------------PLOT OF THE REGION IN ANALYSIS
            if CORR_SHOW==1
                
                figure
                ax(1)=subplot (3,1,1);
                plot(t3(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),ex_clean_ECG(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),'b-',t3(ECGrpeak_pos),ex_clean_ECG(ECGrpeak_pos),'g*');
                legend('ecg','systolic peaks')
                ylabel('mV')
                xlabel('Time(sec)')
                ax(2)=subplot (3,1,2);
                plot(t2(locsplotppg{mm,1}(i:i+wind_size)),test1,'ko',t2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),SIG2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),'b-',ppg80_time,ppg80_values,'r*',t2(locaisspSIG2{mm,1}),SIG2(locaisspSIG2{mm,1},mm),'g*');
                legend('AC/DC PPG (AMPLITUDE)','ppg','20% diastolic peaks', 'Systolic peaks')
                ylabel('amplitude')
                xlabel('Time(sec)')
                ax(3)=subplot (3,1,3);
                plot(t1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),SIG1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),'b-',t1(locsplotabp{mm,1}(i:i+wind_size)),MAP1{l,1},'k+',t1(locsplotabp{mm,1}(i:i+wind_size)),PP1{l,1},'go',t1(locaisdiastSIG1{mm,1}),SIG1(locaisdiastSIG1{mm,1},mm),'r*',t1(locaisspSIG1{mm,1}),SIG1(locaisspSIG1{mm,1},mm),'g*');
                linkaxes(ax,'x');
                legend('bp' ,'MAP','PP','diastolic peaks','systolic peaks');
                ylabel('mmHg')
                xlabel('Time(sec)')
                
                figure;
                subplot 422
                plot((test1),PP1{l,1},'go');
                xlabel('AC/DC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('PP (mmHg)')
                str1 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a, R);
                title(str1)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,PP1{l,1},txt1)
                
                subplot 424
                
                plot((test1),SBP1{l,1},'go');
                xlabel('AC/DC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('SBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a1, R1);
                title(str)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,SBP1{l,1},txt1)
                
                subplot 426
                
                plot((test1),DBP1{l,1},'go');
                xlabel('AC/DC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('DBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a2, R2);
                title(str)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,DBP1{l,1},txt1)
                
                subplot 428
                plot((test1),MAP1{l,1},'go');
                xlabel('AC/DC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('MAP (mmHg)')
                str2 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a4, R4);
                title(str2)
                
                ax(1)=subplot(4,2,[1,3])
                plot(t3(locsplotecg{mm,1}(i:i+wind_size)),test1)
                ylabel('AC/DC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                xlabel('Time(sec)')
                
                ax(2)=subplot(4,2,[5,7])
                plot(t1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),SIG1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)))
                legend('BP');
                ylabel('mmHg')
                xlabel('Time(sec)')
                linkaxes(ax,'x');
                
                mult=1;
                
                pause
            end
            
        end
        
        DNe=Spearman_coef;
    elseif DATA==1
        dados='AC/DC';
        
        test=amp_AC_DC{mm,1};
        
        for i=1:CORR_shift:length(test)-wind_size % VARY SIZE OF WINDOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            test1=test(i:i+wind_size);
            MAP1{l,1}=MAP{l,1}(i:i+wind_size);
            PP1{l,1}=PP{l,1}(i:i+wind_size);
            SBP1{l,1}=SBP{l,1}(i:i+wind_size);
            DBP1{l,1}=DBP{l,1}(i:i+wind_size);
            
            chi_PP1{l,1}=cat(1,chi_PP1{l,1},chi2gof(PP1{l,1}));
            chi_SBP1{l,1}=cat(1,chi_SBP1{l,1},chi2gof(SBP1{l,1}));
            chi_DBP1{l,1}=cat(1,chi_DBP1{l,1},chi2gof(DBP1{l,1}));
            chi_MAP1{l,1}=cat(1,chi_MAP1{l,1},chi2gof(MAP1{l,1}));
            
            cv_PP1{l,1}=cat(1,cv_PP1{l,1},std(PP1{l,1})./mean(PP1{l,1}));
            cv_MAP1{l,1}=cat(1,cv_MAP1{l,1},std(MAP1{l,1})./mean(MAP1{l,1}));
            cv_SBP1{l,1}=cat(1,cv_SBP1{l,1},std(SBP1{l,1})./mean(SBP1{l,1}));
            cv_DBP1{l,1}=cat(1,cv_DBP1{l,1},std(DBP1{l,1})./mean(DBP1{l,1}));
            
            pos_beg=[pos_beg t2(locsplotppg{mm,1}(i))./60];
            pos_end=[pos_end t2(locsplotppg{mm,1}(i+wind_size))./60];
            
            
            %ERROR IN POLYFIT?!?!!?!?!?!?
            [A]=polyfit(test1,PP1{l,1},1); %DBP
            b=A(end);
            a=A(end-1);
            
            [A1]=polyfit(test1,SBP1{l,1},1); %DBP
            b1=A1(end);
            a1=A1(end-1);
            
            [A2]=polyfit(test1,DBP1{l,1},1); %SBP
            b2=A2(end);
            a2=A2(end-1);
            
            [A4]=polyfit(test1,MAP1{l,1},1); %SBP
            b4=A4(end);
            a4=A4(end-1);
            
            R=corr(test1,PP1{l,1},'type','Spearman');
            
            R1=corr(test1,SBP1{l,1},'type','Spearman');
            
            R2=corr(test1,DBP1{l,1},'type','Spearman');
            
            R4=corr(test1,MAP1{l,1},'type','Spearman');
            
            %------------------------------------- PLOTS WITH CORRELATION COEFS
            %AND LIN_REG COEFS
            
            lin_reg_coef=cat(2,lin_reg_coef,[A(1,1);A1(1,1);A2(1,1);A4(1,1)]);
            b_reg_coef=cat(2,b_reg_coef,[A(end);A1(end);A2(end);A4(end)]);
            
            Spearman_coef=cat(2,Spearman_coef,[R;R1;R2;R4]);
            
            rele_wind_acdc=cell(1,4);
            
            feat=1;
            pp_pp=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            feat=2;
            pp_sbp=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            feat=3;
            pp_dbp=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            feat=4;
            pp_map=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            
            rele_wind_acdc{1,1}=pp_pp;
            rele_wind_acdc{1,2}=pp_sbp;
            rele_wind_acdc{1,3}=pp_dbp;
            rele_wind_acdc{1,4}=pp_map;
            
            
            %--------------------------------------PLOT OF THE REGION IN ANALYSIS
            if CORR_SHOW==1
                
                figure
                ax(1)=subplot (3,1,1);
                plot(t3(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),ex_clean_ECG(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),'b-',t3(ECGrpeak_pos),ex_clean_ECG(ECGrpeak_pos),'g*');
                legend('ecg','systolic peaks')
                ylabel('mV')
                xlabel('Time(sec)')
                ax(2)=subplot (3,1,2);
                plot(t2(locsplotppg{mm,1}(i:i+wind_size)),test1,'ko',t2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),SIG2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),'b-',ppg80_time,ppg80_values,'r*',t2(locaisspSIG2{mm,1}),SIG2(locaisspSIG2{mm,1},mm),'g*');
                legend('AC/DC PPG (AMPLITUDE)','ppg','20% diastolic peaks', 'Systolic peaks')
                ylabel('amplitude')
                xlabel('Time(sec)')
                ax(3)=subplot (3,1,3);
                plot(t1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),SIG1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),'b-',t1(locsplotabp{mm,1}(i:i+wind_size)),MAP1{l,1},'k+',t1(locsplotabp{mm,1}(i:i+wind_size)),PP1{l,1},'go',t1(locaisdiastSIG1{mm,1}),SIG1(locaisdiastSIG1{mm,1},mm),'r*',t1(locaisspSIG1{mm,1}),SIG1(locaisspSIG1{mm,1},mm),'g*');
                linkaxes(ax,'x');
                legend('bp' ,'MAP','PP','diastolic peaks','systolic peaks');
                ylabel('mmHg')
                xlabel('Time(sec)')
                
                figure;
                subplot 422
                plot((test1),PP1{l,1},'go');
                xlabel('AC/DC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('PP (mmHg)')
                str1 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a, R);
                title(str1)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,PP1{l,1},txt1)
                
                subplot 424
                
                plot((test1),SBP1{l,1},'go');
                xlabel('AC/DC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('SBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a1, R1);
                title(str)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,SBP1{l,1},txt1)
                
                subplot 426
                
                plot((test1),DBP1{l,1},'go');
                xlabel('AC/DC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('DBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a2, R2);
                title(str)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,DBP1{l,1},txt1)
                
                subplot 428
                plot((test1),MAP1{l,1},'go');
                xlabel('AC/DC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('MAP (mmHg)')
                str2 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a4, R4);
                title(str2)
                
                ax(1)=subplot(4,2,[1,3])
                plot(t3(locsplotecg{mm,1}(i:i+wind_size)),test1)
                ylabel('AC/DC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                xlabel('Time(sec)')
                
                ax(2)=subplot(4,2,[5,7])
                plot(t1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),SIG1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)))
                legend('BP');
                ylabel('mmHg')
                xlabel('Time(sec)')
                linkaxes(ax,'x');
                
                mult=1;
                
                %                 BlandAltman(mult*PP1{l,1}, mult*test1,{'PP (mmHg)','AC/DC PPG (AMPLITUDE)'},'baStatsMode','Non-parametric','corrInfo',{'RMSE','rho','eq'},'baInfo','IQR');
                % %                 +57
                %                 BlandAltman(mult*SBP1{l,1}, mult*test1,{'SBP (mmHg)','AC/DC PPG (AMPLITUDE)'},'baStatsMode','Non-parametric','corrInfo',{'RMSE','rho','eq'},'baInfo','IQR');
                % %                 +115
                %                 BlandAltman(mult*DBP1{l,1}, mult*test1,{'DBP (mmHg)','AC/DC PPG (AMPLITUDE)'},'baStatsMode','Non-parametric','corrInfo',{'RMSE','rho','eq'},'baInfo','IQR');
                % %                 +57
                %                 BlandAltman(mult*MAP1{l,1}, mult*test1,{'MAP (mmHg)','AC/DC PPG (AMPLITUDE)'},'baStatsMode','Non-parametric','corrInfo',{'RMSE','rho','eq'},'baInfo','IQR');
                % %                 +76
                pause
            end
            
            
        end
        
        ACDCe=Spearman_coef;
    elseif DATA==3
        dados='DC';
        
        test=amp_DC{mm,1};
        
        for i=1:CORR_shift:length(test)-wind_size % VARY SIZE OF WINDOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            test1=test(i:i+wind_size);
            PP1{l,1}=PP{l,1}(i:i+wind_size);
            SBP1{l,1}=SBP{l,1}(i:i+wind_size);
            DBP1{l,1}=DBP{l,1}(i:i+wind_size);
            MAP1{l,1}=MAP{l,1}(i:i+wind_size);
            chi_PP1{l,1}=cat(1,chi_PP1{l,1},chi2gof(PP1{l,1}));
            chi_SBP1{l,1}=cat(1,chi_SBP1{l,1},chi2gof(SBP1{l,1}));
            chi_DBP1{l,1}=cat(1,chi_DBP1{l,1},chi2gof(DBP1{l,1}));
            chi_MAP1{l,1}=cat(1,chi_MAP1{l,1},chi2gof(MAP1{l,1}));
            
            cv_PP1{l,1}=cat(1,cv_PP1{l,1},std(PP1{l,1})./mean(PP1{l,1}));
            cv_MAP1{l,1}=cat(1,cv_MAP1{l,1},std(MAP1{l,1})./mean(MAP1{l,1}));
            cv_SBP1{l,1}=cat(1,cv_SBP1{l,1},std(SBP1{l,1})./mean(SBP1{l,1}));
            cv_DBP1{l,1}=cat(1,cv_DBP1{l,1},std(DBP1{l,1})./mean(DBP1{l,1}));
            
            pos_beg=[pos_beg t2(locsplotppg{mm,1}(i))./60];
            pos_end=[pos_end t2(locsplotppg{mm,1}(i+wind_size))./60];
            
            %--------------------------------------PLOT OF THE REGION IN ANALYSIS
            if CORR_SHOW==1
                
                figure
                ax(1)=subplot (3,1,1);
                plot(t3(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),ex_clean_ECG(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),'b-',t3(ECGrpeak_pos),ex_clean_ECG(ECGrpeak_pos),'g*');
                legend('ecg','systolic peaks')
                ylabel('mV')
                xlabel('Time(sec)')
                ax(2)=subplot (3,1,2);
                plot(t2(locsplotppg{mm,1}(i:i+wind_size)),test1,'ko',t2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),SIG2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),'b-',ppg80_time,ppg80_values,'r*',t2(locaisspSIG2{num,1}),SIG2(locaisspSIG2{num,1},num),'g*');
                legend('DC PPG (AMPLITUDE)','ppg','20% diastolic peaks', 'Systolic peaks')
                ylabel('amplitude')
                xlabel('Time(sec)')
                ax(3)=subplot (3,1,3);
                plot(t1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),SIG1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),'b-',t1(locsplotabp{mm,1}(i:i+wind_size)),MAP1{l,1},'k+',t1(locsplotabp{mm,1}(i:i+wind_size)),PP1{l,1},'go',t1(locaisdiastSIG1{num,1}),SIG1(locaisdiastSIG1{num,1},num),'r*',t1(locaisspSIG1{num,1}),SIG1(locaisspSIG1{num,1},num),'g*');
                linkaxes(ax,'x');
                legend('bp' ,'MAP','PP','diastolic peaks','systolic peaks');
                ylabel('mmHg')
                xlabel('Time(sec)')
            end
            %ERROR IN POLYFIT?!?!!?!?!?!?
            [A]=polyfit(test1,PP1{l,1},1); %DBP
            b=A(end);
            a=A(end-1);
            
            [A1]=polyfit(test1,SBP1{l,1},1); %DBP
            b1=A1(end);
            a1=A1(end-1);
            
            [A2]=polyfit(test1,DBP1{l,1},1); %SBP
            b2=A2(end);
            a2=A2(end-1);
            
            
            [A4]=polyfit(test1,MAP1{l,1},1); %SBP
            b4=A4(end);
            a4=A4(end-1);
            
            R=corr(test1,PP1{l,1},'type','Spearman');
            
            %         [R1]=corrcoef(test1,SBP1{l,1}); %DBP
            R1=corr(test1,SBP1{l,1},'type','Spearman');
            
            %         fprintf ('SDBP vs RI');
            %         [R2]=corrcoef(test1,DBP1{l,1}); %SBP
            R2=corr(test1,DBP1{l,1},'type','Spearman');
            
            
            %         [R4]=corrcoef(test1,MAP1{l,1}); %SBP
            R4=corr(test1,MAP1{l,1},'type','Spearman');
            
            %------------------------------------- PLOTS WITH CORRELATION COEFS
            %AND LIN_REG COEFS
            
            lin_reg_coef=cat(2,lin_reg_coef,[A(1,1);A1(1,1);A2(1,1);A4(1,1)]);
            b_reg_coef=cat(2,b_reg_coef,[A(end);A1(end);A2(end);A4(end)]);
            
            Spearman_coef=cat(2,Spearman_coef,[R;R1;R2;R4]);
            
            rele_wind_dc=cell(1,4);
            
            feat=1;
            pp_pp=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            feat=2;
            pp_sbp=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            feat=3;
            pp_dbp=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            feat=4;
            pp_map=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            
            rele_wind_dc{1,1}=pp_pp;
            rele_wind_dc{1,2}=pp_sbp;
            rele_wind_dc{1,3}=pp_dbp;
            rele_wind_dc{1,4}=pp_map;
            
            if CORR_SHOW==1
                
                figure;
                subplot 422
                plot((test1),PP1{l,1},'go');
                xlabel('DC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('PP (mmHg)')
                str1 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a, R);
                title(str1)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,PP1{l,1},txt1)
                
                subplot 424
                
                plot((test1),SBP1{l,1},'go');
                xlabel('DC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('SBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a1, R1);
                title(str)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,SBP1{l,1},txt1)
                
                subplot 426
                
                plot((test1),DBP1{l,1},'go');
                xlabel('DC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('DBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a2, R2);
                title(str)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,DBP1{l,1},txt1)
                
                subplot 428
                plot((test1),MAP1{l,1},'go');
                xlabel('DC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('MAP (mmHg)')
                str2 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a4, R4);
                title(str2)
                
                ax(1)=subplot(4,2,[1,3])
                plot(t3(locsplotecg{mm,1}(i:i+wind_size)),test1)
                ylabel('DC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                xlabel('Time(sec)')
                
                ax(2)=subplot(4,2,[5,7])
                plot(t1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),SIG1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)))
                legend('BP');
                ylabel('mmHg')
                xlabel('Time(sec)')
                linkaxes(ax,'x');
                
                pause
                
            end
            
            
        end
        DCe=Spearman_coef;
        
    elseif DATA==4
        dados='PAT ';
        test=PAT{mm,1};
        
        for i=1:CORR_shift:length(test)-wind_size % VARY SIZE OF WINDOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            test1=test(i:i+wind_size);
            PP1{l,1}=PP{l,1}(i:i+wind_size);
            SBP1{l,1}=SBP{l,1}(i:i+wind_size);
            DBP1{l,1}=DBP{l,1}(i:i+wind_size);
            MAP1{l,1}=MAP{l,1}(i:i+wind_size);
            
            chi_PP1{l,1}=cat(1,chi_PP1{l,1},chi2gof(PP1{l,1}));
            chi_SBP1{l,1}=cat(1,chi_SBP1{l,1},chi2gof(SBP1{l,1}));
            chi_DBP1{l,1}=cat(1,chi_DBP1{l,1},chi2gof(DBP1{l,1}));
            chi_MAP1{l,1}=cat(1,chi_MAP1{l,1},chi2gof(MAP1{l,1}));
            
            cv_PP1{l,1}=cat(1,cv_PP1{l,1},std(PP1{l,1})./mean(PP1{l,1}));
            cv_MAP1{l,1}=cat(1,cv_MAP1{l,1},std(MAP1{l,1})./mean(MAP1{l,1}));
            cv_SBP1{l,1}=cat(1,cv_SBP1{l,1},std(SBP1{l,1})./mean(SBP1{l,1}));
            cv_DBP1{l,1}=cat(1,cv_DBP1{l,1},std(DBP1{l,1})./mean(DBP1{l,1}));
            
            pos_beg=[pos_beg t3(locsplotecg{mm,1}(i))./60];
            pos_end=[pos_end t3(locsplotecg{mm,1}(i+wind_size))./60];
            %--------------------------------------PLOT OF THE REGION IN ANALYSIS
            if CORR_SHOW==1
                
                figure
                ax(1)=subplot (3,1,1);
                plot(t3(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),ex_clean_ECG(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),'b-',t3(ECGrpeak_pos),ex_clean_ECG(ECGrpeak_pos),'g*',t3(locsplotecg{mm,1}(i:i+wind_size)),test1,'ro');
                legend('ecg','systolic peaks','d(PAT)/dt')
                ylabel('mV')
                xlabel('Time(sec)')
                ax(2)=subplot (3,1,2);
                plot(t2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),SIG2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),'b-',ppg80_time,ppg80_values,'r*',t2(locaisspSIG2{mm,1}),SIG2(locaisspSIG2{mm,1},mm),'g*');
                legend('ppg','20% diastolic peaks', 'Systolic peaks')
                ylabel('amplitude')
                xlabel('Time(sec)')
                ax(3)=subplot (3,1,3);
                plot(t1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),SIG1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),'b-',t1(locsplotabp{mm,1}(i:i+wind_size)),MAP1{l,1},'k+',t1(locsplotabp{mm,1}(i:i+wind_size)),PP1{l,1},'go',t1(locaisdiastSIG1{mm,1}),SIG1(locaisdiastSIG1{mm,1},mm),'r*',t1(locaisspSIG1{mm,1}),SIG1(locaisspSIG1{mm,1},mm),'g*');
                linkaxes(ax,'x');
                legend('bp' ,'d(MAP)/dt','d(PP)/dt','diastolic peaks','systolic peaks');
                ylabel('mmHg')
                xlabel('Time(sec)')
            end
            %ERROR IN POLYFIT?!?!!?!?!?!?
            [A]=polyfit(test1,PP1{l,1},1); %DBP
            b=A(end);
            a=A(end-1);
            
            [A1]=polyfit(test1,SBP1{l,1},1); %DBP
            b1=A1(end);
            a1=A1(end-1);
            
            [A2]=polyfit(test1,DBP1{l,1},1); %SBP
            b2=A2(end);
            a2=A2(end-1);
            
            [A4]=polyfit(test1,MAP1{l,1},1); %SBP
            b4=A4(end);
            a4=A4(end-1);
            
            
            R=corr(test1,PP1{l,1},'type','Spearman');
            
            R1=corr(test1,SBP1{l,1},'type','Spearman');
            
            R2=corr(test1,DBP1{l,1},'type','Spearman');
            
            R4=corr(test1,MAP1{l,1},'type','Spearman');
            
            %------------------------------------- PLOTS WITH CORRELATION COEFS
            %AND LIN_REG COEFS
            
            lin_reg_coef=cat(2,lin_reg_coef,[A(1,1);A1(1,1);A2(1,1);A4(1,1)]);
            b_reg_coef=cat(2,b_reg_coef,[A(end);A1(end);A2(end);A4(end)]);
            
            Spearman_coef=cat(2,Spearman_coef,[R;R1;R2;R4]);
            
            rele_wind_pat=cell(1,4);
            
            
            feat=1;
            pp_pp=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            feat=2;
            pp_sbp=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            feat=3;
            pp_dbp=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            feat=4;
            pp_map=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            
            rele_wind_pat{1,1}=pp_pp;
            rele_wind_pat{1,2}=pp_sbp;
            rele_wind_pat{1,3}=pp_dbp;
            rele_wind_pat{1,4}=pp_map;
            
            
            if CORR_SHOW==1
                
                figure;
                subplot 411
                plot((test1),PP1{l,1},'go');
                xlabel('PAT (SEC)')  %xlabel('AC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('PP (mmHg)')
                str1 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a, R);
                title(str1)
                
                subplot 412
                
                plot((test1),SBP1{l,1},'go');
                xlabel('PAT (SEC)')
                ylabel('SBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a1, R1);
                title(str)
                
                subplot 413
                
                plot((test1),DBP1{l,1},'go');
                xlabel('PAT (SEC)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('DBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a2, R2);
                title(str)
                
                subplot 414
                plot((test1),MAP1{l,1},'go');
                xlabel('PAT (SEC)')% xlabel('AC PPG (AMPLITUDE)') % xlabel('DC PPG')
                ylabel('MAP (mmHg)')
                str2 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a4, R4);
                title(str2)
                
                pause
                
                
            end
            
            
        end
        PATe=Spearman_coef;
        
    elseif DATA==5
        dados='PAT ';
        test=PAT2{mm,1};
        
        for i=1:CORR_shift:length(test)-wind_size % VARY SIZE OF WINDOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            test1=test(i:i+wind_size);
            PP1{l,1}=PP{l,1}(i:i+wind_size);
            SBP1{l,1}=SBP{l,1}(i:i+wind_size);
            DBP1{l,1}=DBP{l,1}(i:i+wind_size);
            MAP1{l,1}=MAP{l,1}(i:i+wind_size);
            
            chi_PP1{l,1}=cat(1,chi_PP1{l,1},chi2gof(PP1{l,1}));
            chi_SBP1{l,1}=cat(1,chi_SBP1{l,1},chi2gof(SBP1{l,1}));
            chi_DBP1{l,1}=cat(1,chi_DBP1{l,1},chi2gof(DBP1{l,1}));
            chi_MAP1{l,1}=cat(1,chi_MAP1{l,1},chi2gof(MAP1{l,1}));
            
            cv_PP1{l,1}=cat(1,cv_PP1{l,1},std(PP1{l,1})./mean(PP1{l,1}));
            cv_MAP1{l,1}=cat(1,cv_MAP1{l,1},std(MAP1{l,1})./mean(MAP1{l,1}));
            cv_SBP1{l,1}=cat(1,cv_SBP1{l,1},std(SBP1{l,1})./mean(SBP1{l,1}));
            cv_DBP1{l,1}=cat(1,cv_DBP1{l,1},std(DBP1{l,1})./mean(DBP1{l,1}));
            
            pos_beg=[pos_beg t3(locsplotecg{mm,1}(i))./60];
            pos_end=[pos_end t3(locsplotecg{mm,1}(i+wind_size))./60];
            %--------------------------------------PLOT OF THE REGION IN ANALYSIS
            if CORR_SHOW==1
                
                figure
                ax(1)=subplot (3,1,1);
                plot(t3(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),ex_clean_ECG(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),'b-',t3(ECGrpeak_pos),ex_clean_ECG(ECGrpeak_pos),'g*',t3(locsplotecg{mm,1}(i:i+wind_size)),test1,'ro');
                legend('ecg','systolic peaks','d(PAT)/dt')
                ylabel('mV')
                xlabel('Time(sec)')
                ax(2)=subplot (3,1,2);
                plot(t2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),SIG2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),'b-',ppg80_time,ppg80_values,'r*',t2(locaisspSIG2{mm,1}),SIG2(locaisspSIG2{mm,1},mm),'g*');
                legend('ppg','20% diastolic peaks', 'Systolic peaks')
                ylabel('amplitude')
                xlabel('Time(sec)')
                ax(3)=subplot (3,1,3);
                plot(t1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),SIG1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),'b-',t1(locsplotabp{mm,1}(i:i+wind_size)),MAP1{l,1},'k+',t1(locsplotabp{mm,1}(i:i+wind_size)),PP1{l,1},'go',t1(locaisdiastSIG1{mm,1}),SIG1(locaisdiastSIG1{mm,1},mm),'r*',t1(locaisspSIG1{mm,1}),SIG1(locaisspSIG1{mm,1},mm),'g*');
                linkaxes(ax,'x');
                legend('bp' ,'d(MAP)/dt','d(PP)/dt','diastolic peaks','systolic peaks');
                ylabel('mmHg')
                xlabel('Time(sec)')
            end
            %ERROR IN POLYFIT?!?!!?!?!?!?
            [A]=polyfit(test1,PP1{l,1},1); %DBP
            b=A(end);
            a=A(end-1);
            
            [A1]=polyfit(test1,SBP1{l,1},1); %DBP
            b1=A1(end);
            a1=A1(end-1);
            
            [A2]=polyfit(test1,DBP1{l,1},1); %SBP
            b2=A2(end);
            a2=A2(end-1);
            
            [A4]=polyfit(test1,MAP1{l,1},1); %SBP
            b4=A4(end);
            a4=A4(end-1);
            
            
            R=corr(test1,PP1{l,1},'type','Spearman');
            
            R1=corr(test1,SBP1{l,1},'type','Spearman');
            
            R2=corr(test1,DBP1{l,1},'type','Spearman');
            
            R4=corr(test1,MAP1{l,1},'type','Spearman');
            
            %------------------------------------- PLOTS WITH CORRELATION COEFS
            %AND LIN_REG COEFS
            
            lin_reg_coef=cat(2,lin_reg_coef,[A(1,1);A1(1,1);A2(1,1);A4(1,1)]);
            b_reg_coef=cat(2,b_reg_coef,[A(end);A1(end);A2(end);A4(end)]);
            
            Spearman_coef=cat(2,Spearman_coef,[R;R1;R2;R4]);
            
            rele_wind_pat=cell(1,4);
            
            
            feat=1;
            pp_pp=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            feat=2;
            pp_sbp=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            feat=3;
            pp_dbp=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            feat=4;
            pp_map=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            
            rele_wind_pat{1,1}=pp_pp;
            rele_wind_pat{1,2}=pp_sbp;
            rele_wind_pat{1,3}=pp_dbp;
            rele_wind_pat{1,4}=pp_map;
            
            
            if CORR_SHOW==1
                
                figure;
                subplot 411
                plot((test1),PP1{l,1},'go');
                xlabel('PAT (SEC)')  %xlabel('AC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('PP (mmHg)')
                str1 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a, R);
                title(str1)
                
                subplot 412
                
                plot((test1),SBP1{l,1},'go');
                xlabel('PAT (SEC)')
                ylabel('SBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a1, R1);
                title(str)
                
                subplot 413
                
                plot((test1),DBP1{l,1},'go');
                xlabel('PAT (SEC)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('DBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a2, R2);
                title(str)
                
                subplot 414
                plot((test1),MAP1{l,1},'go');
                xlabel('PAT (SEC)')% xlabel('AC PPG (AMPLITUDE)') % xlabel('DC PPG')
                ylabel('MAP (mmHg)')
                str2 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a4, R4);
                title(str2)
                
                pause
                
                
            end
            
            
        end
        PAT2e=Spearman_coef;
        
    elseif DATA==7
        test=RI{mm,1};
        dados='RI';
        
        for i=1:CORR_shift:length(test)-wind_size % VARY SIZE OF WINDOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            test1=test(i:i+wind_size);
            PP1{l,1}=PP{l,1}(i:i+wind_size);
            SBP1{l,1}=SBP{l,1}(i:i+wind_size);
            DBP1{l,1}=DBP{l,1}(i:i+wind_size);
            MAP1{l,1}=MAP{l,1}(i:i+wind_size);
            
            chi_PP1{l,1}=cat(1,chi_PP1{l,1},chi2gof(PP1{l,1}));
            chi_SBP1{l,1}=cat(1,chi_SBP1{l,1},chi2gof(SBP1{l,1}));
            chi_DBP1{l,1}=cat(1,chi_DBP1{l,1},chi2gof(DBP1{l,1}));
            chi_MAP1{l,1}=cat(1,chi_MAP1{l,1},chi2gof(MAP1{l,1}));
            
            cv_PP1{l,1}=cat(1,cv_PP1{l,1},std(PP1{l,1})./mean(PP1{l,1}));
            cv_MAP1{l,1}=cat(1,cv_MAP1{l,1},std(MAP1{l,1})./mean(MAP1{l,1}));
            cv_SBP1{l,1}=cat(1,cv_SBP1{l,1},std(SBP1{l,1})./mean(SBP1{l,1}));
            cv_DBP1{l,1}=cat(1,cv_DBP1{l,1},std(DBP1{l,1})./mean(DBP1{l,1}));
            
            pos_beg=[pos_beg t3(locsplotecg{mm,1}(i))./60];
            pos_end=[pos_end t3(locsplotecg{mm,1}(i+wind_size))./60];
            %--------------------------------------PLOT OF THE REGION IN ANALYSIS
            if CORR_SHOW==1
                
                figure
                ax(1)=subplot (3,1,1);
                plot(t3(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),ex_clean_ECG(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),'b-',t3(ECGrpeak_pos),ex_clean_ECG(ECGrpeak_pos),'g*');
                legend('ecg','systolic peaks')
                ylabel('mV')
                xlabel('Time(sec)')
                ax(2)=subplot (3,1,2);
                plot(t2(locsplotppg{mm,1}(i:i+wind_size)),test1,'ko',t2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),SIG2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),'b-',ppg80_time,ppg80_values,'r*',t2(locaisspSIG2{num,1}),SIG2(locaisspSIG2{num,1},num),'g*');
                legend('AC PPG (AMPLITUDE)','ppg','20% diastolic peaks', 'Systolic peaks')
                ylabel('amplitude')
                xlabel('Time(sec)')
                ax(3)=subplot (3,1,3);
                plot(t1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),SIG1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),'b-',t1(locsplotabp{mm,1}(i:i+wind_size)),MAP1{l,1},'k+',t1(locsplotabp{mm,1}(i:i+wind_size)),PP1{l,1},'go',t1(locaisdiastSIG1{num,1}),SIG1(locaisdiastSIG1{num,1},num),'r*',t1(locaisspSIG1{num,1}),SIG1(locaisspSIG1{num,1},num),'g*');
                linkaxes(ax,'x');
                legend('bp' ,'MAP','PP','diastolic peaks','systolic peaks');
                ylabel('mmHg')
                xlabel('Time(sec)')
            end
            %ERROR IN POLYFIT?!?!!?!?!?!?
            [A]=polyfit(test1,PP1{l,1},1); %DBP
            b=A(end);
            a=A(end-1);
            
            [A1]=polyfit(test1,SBP1{l,1},1); %DBP
            b1=A1(end);
            a1=A1(end-1);
            
            [A2]=polyfit(test1,DBP1{l,1},1); %SBP
            b2=A2(end);
            a2=A2(end-1);
            
            [A4]=polyfit(test1,MAP1{l,1},1); %SBP
            b4=A4(end);
            a4=A4(end-1);
            
            %PP
            R=corr(test1,PP1{l,1},'type','Spearman');
            
            %DBP
            R1=corr(test1,SBP1{l,1},'type','Spearman');
            
            %SBP
            R2=corr(test1,DBP1{l,1},'type','Spearman');
            
            %SBP
            R4=corr(test1,MAP1{l,1},'type','Spearman');
            
            
            %------------------------------------- PLOTS WITH CORRELATION COEFS
            %AND LIN_REG COEFS
            
            lin_reg_coef=cat(2,lin_reg_coef,[A(1,1);A1(1,1);A2(1,1);A4(1,1)]);
            b_reg_coef=cat(2,b_reg_coef,[A(end);A1(end);A2(end);A4(end)]);
            
            Spearman_coef=cat(2,Spearman_coef,[R;R1;R2;R4]);
            
            rele_wind_ri=cell(1,4);
            
            feat=1;
            nn_pp=find(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)<0);
            feat=2;
            nn_sbp=find(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)<0);
            feat=3;
            nn_dbp=find(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)<0);
            feat=4;
            nn_map=find(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)<0);
            
            rele_wind_ri{1,1}=nn_pp;
            rele_wind_ri{1,2}=nn_sbp;
            rele_wind_ri{1,3}=nn_dbp;
            rele_wind_ri{1,4}=nn_map;
            
            if CORR_SHOW==1
                
                figure;
                subplot 422
                plot((test1),PP1{l,1},'go');
                xlabel('RI') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('PP (mmHg)')
                str1 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a, R);
                title(str1)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,PP1{l,1},txt1)
                
                subplot 424
                
                plot((test1),SBP1{l,1},'go');
                xlabel('RI') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('SBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a1, R1);
                title(str)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,SBP1{l,1},txt1)
                
                subplot 426
                
                plot((test1),DBP1{l,1},'go');
                xlabel('RI') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('DBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a2, R2);
                title(str)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,DBP1{l,1},txt1)
                
                subplot 428
                plot((test1),MAP1{l,1},'go');
                xlabel('RI') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('MAP (mmHg)')
                str2 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a4, R4);
                title(str2)
                
                ax(1)=subplot(4,2,[1,3])
                plot(t3(locsplotecg{mm,1}(i:i+wind_size)),test1)
                ylabel('AC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                xlabel('Time(sec)')
                
                ax(2)=subplot(4,2,[5,7])
                plot(t1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),SIG1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)))
                legend('BP');
                ylabel('mmHg')
                xlabel('Time(sec)')
                linkaxes(ax,'x');
                
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,MAP1{l,1},txt1)
                
                pause
                
            end
            
            
        end
        RIe=Spearman_coef;
    elseif DATA==9
        test=SI;
        dados='SI';
        
        for i=1:CORR_shift:length(test)-wind_size % VARY SIZE OF WINDOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            test1=test(i:i+wind_size);
            PP1{l,1}=PP{l,1}(i:i+wind_size);
            SBP1{l,1}=SBP{l,1}(i:i+wind_size);
            DBP1{l,1}=DBP{l,1}(i:i+wind_size);
            MAP1{l,1}=MAP{l,1}(i:i+wind_size);
            
            chi_PP1{l,1}=cat(1,chi_PP1{l,1},chi2gof(PP1{l,1}));
            chi_SBP1{l,1}=cat(1,chi_SBP1{l,1},chi2gof(SBP1{l,1}));
            chi_DBP1{l,1}=cat(1,chi_DBP1{l,1},chi2gof(DBP1{l,1}));
            chi_MAP1{l,1}=cat(1,chi_MAP1{l,1},chi2gof(MAP1{l,1}));
            
            cv_PP1{l,1}=cat(1,cv_PP1{l,1},std(PP1{l,1})./mean(PP1{l,1}));
            cv_MAP1{l,1}=cat(1,cv_MAP1{l,1},std(MAP1{l,1})./mean(MAP1{l,1}));
            cv_SBP1{l,1}=cat(1,cv_SBP1{l,1},std(SBP1{l,1})./mean(SBP1{l,1}));
            cv_DBP1{l,1}=cat(1,cv_DBP1{l,1},std(DBP1{l,1})./mean(DBP1{l,1}));
            
            pos_beg=[pos_beg t3(locsplotecg{mm,1}(i))./60];
            pos_end=[pos_end t3(locsplotecg{mm,1}(i+wind_size))./60];
            %--------------------------------------PLOT OF THE REGION IN ANALYSIS
            if CORR_SHOW==1
                
                figure
                ax(1)=subplot (3,1,1);
                plot(t3(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),ex_clean_ECG(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),'b-',t3(ECGrpeak_pos),ex_clean_ECG(ECGrpeak_pos),'g*');
                legend('ecg','systolic peaks')
                ylabel('mV')
                xlabel('Time(sec)')
                ax(2)=subplot (3,1,2);
                plot(t2(locsplotppg{mm,1}(i:i+wind_size)),test1,'ko',t2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),SIG2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),'b-',ppg80_time,ppg80_values,'r*',t2(locaisspSIG2{num,1}),SIG2(locaisspSIG2{num,1},num),'g*');
                legend('AC PPG (AMPLITUDE)','ppg','20% diastolic peaks', 'Systolic peaks')
                ylabel('amplitude')
                xlabel('Time(sec)')
                ax(3)=subplot (3,1,3);
                plot(t1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),SIG1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),'b-',t1(locsplotabp{mm,1}(i:i+wind_size)),MAP1{l,1},'k+',t1(locsplotabp{mm,1}(i:i+wind_size)),PP1{l,1},'go',t1(locaisdiastSIG1{num,1}),SIG1(locaisdiastSIG1{num,1},num),'r*',t1(locaisspSIG1{num,1}),SIG1(locaisspSIG1{num,1},num),'g*');
                linkaxes(ax,'x');
                legend('bp' ,'MAP','PP','diastolic peaks','systolic peaks');
                ylabel('mmHg')
                xlabel('Time(sec)')
            end
            %ERROR IN POLYFIT?!?!!?!?!?!?
            [A]=polyfit(test1,PP1{l,1},1); %DBP
            b=A(end);
            a=A(end-1);
            
            [A1]=polyfit(test1,SBP1{l,1},1); %DBP
            b1=A1(end);
            a1=A1(end-1);
            
            [A2]=polyfit(test1,DBP1{l,1},1); %SBP
            b2=A2(end);
            a2=A2(end-1);
            
            [A4]=polyfit(test1,MAP1{l,1},1); %SBP
            b4=A4(end);
            a4=A4(end-1);
            
            %PP
            R=corr(test1,PP1{l,1},'type','Spearman');
            
            %DBP
            R1=corr(test1,SBP1{l,1},'type','Spearman');
            
            %SBP
            R2=corr(test1,DBP1{l,1},'type','Spearman');
            
            %SBP
            R4=corr(test1,MAP1{l,1},'type','Spearman');
            
            
            %------------------------------------- PLOTS WITH CORRELATION COEFS
            %AND LIN_REG COEFS
            
            lin_reg_coef=cat(2,lin_reg_coef,[A(1,1);A1(1,1);A2(1,1);A4(1,1)]);
            b_reg_coef=cat(2,b_reg_coef,[A(end);A1(end);A2(end);A4(end)]);
            
            Spearman_coef=cat(2,Spearman_coef,[R;R1;R2;R4]);
            
            rele_wind_si=cell(1,4);
            
            feat=1;
            nn_pp=find(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)<0);
            feat=2;
            nn_sbp=find(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)<0);
            feat=3;
            nn_dbp=find(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)<0);
            feat=4;
            nn_map=find(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)<0);
            
            rele_wind_si{1,1}=nn_pp;
            rele_wind_si{1,2}=nn_sbp;
            rele_wind_si{1,3}=nn_dbp;
            rele_wind_si{1,4}=nn_map;
            
            if CORR_SHOW==1
                
                figure;
                subplot 422
                plot((test1),PP1{l,1},'go');
                xlabel('RI') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('PP (mmHg)')
                str1 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a, R);
                title(str1)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,PP1{l,1},txt1)
                
                subplot 424
                
                plot((test1),SBP1{l,1},'go');
                xlabel('RI') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('SBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a1, R1);
                title(str)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,SBP1{l,1},txt1)
                
                subplot 426
                
                plot((test1),DBP1{l,1},'go');
                xlabel('RI') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('DBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a2, R2);
                title(str)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,DBP1{l,1},txt1)
                
                subplot 428
                plot((test1),MAP1{l,1},'go');
                xlabel('RI') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('MAP (mmHg)')
                str2 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a4, R4);
                title(str2)
                
                ax(1)=subplot(4,2,[1,3])
                plot(t3(locsplotecg{mm,1}(i:i+wind_size)),test1)
                ylabel('AC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                xlabel('Time(sec)')
                
                ax(2)=subplot(4,2,[5,7])
                plot(t1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),SIG1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)))
                legend('BP');
                ylabel('mmHg')
                xlabel('Time(sec)')
                linkaxes(ax,'x');
                
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,MAP1{l,1},txt1)
                
                pause
                
            end
            
            
        end
        SIe=Spearman_coef;
        
    elseif DATA==6
        dados='HBI';
        
        test=HBI;
        %         figure
        %         subplot 221
        %         histogram(test-PP{l,1})
        %         legend('HR-PP');
        %
        %         subplot 222
        %         histogram(test-SBP{l,1})
        %         legend('HR-SBP');
        %
        %         subplot 223
        %         histogram(test-DBP{l,1})
        %         legend('HR-DBP');
        %
        %         subplot 224
        %         histogram(test-MAP{l,1})
        %         legend('HR-MAP');
        
        for i=1:CORR_shift:length(test)-wind_size % VARY SIZE OF WINDOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            test1=test(i:i+wind_size);
            
            PP1{l,1}=PP{l,1}(i:i+wind_size);
            SBP1{l,1}=SBP{l,1}(i:i+wind_size);
            DBP1{l,1}=DBP{l,1}(i:i+wind_size);
            MAP1{l,1}=MAP{l,1}(i:i+wind_size);
            
            chi_PP1{l,1}=cat(1,chi_PP1{l,1},chi2gof(PP1{l,1}));
            chi_SBP1{l,1}=cat(1,chi_SBP1{l,1},chi2gof(SBP1{l,1}));
            chi_DBP1{l,1}=cat(1,chi_DBP1{l,1},chi2gof(DBP1{l,1}));
            chi_MAP1{l,1}=cat(1,chi_MAP1{l,1},chi2gof(MAP1{l,1}));
            
            cv_PP1{l,1}=cat(1,cv_PP1{l,1},std(PP1{l,1})./mean(PP1{l,1}));
            cv_MAP1{l,1}=cat(1,cv_MAP1{l,1},std(MAP1{l,1})./mean(MAP1{l,1}));
            cv_SBP1{l,1}=cat(1,cv_SBP1{l,1},std(SBP1{l,1})./mean(SBP1{l,1}));
            cv_DBP1{l,1}=cat(1,cv_DBP1{l,1},std(DBP1{l,1})./mean(DBP1{l,1}));
            
            pos_beg=[pos_beg t3(locsplotecg{mm,1}(i))./60];
            pos_end=[pos_end t3(locsplotecg{mm,1}(i+wind_size))./60];
            
            %ERROR IN POLYFIT?!?!!?!?!?!?
            [A]=polyfit(test1,PP1{l,1},1); %DBP
            b=A(end);
            a=A(end-1);
            
            [A1]=polyfit(test1,SBP1{l,1},1); %DBP
            b1=A1(end);
            a1=A1(end-1);
            
            [A2]=polyfit(test1,DBP1{l,1},1); %SBP
            b2=A2(end);
            a2=A2(end-1);
            
            [A4]=polyfit(test1,MAP1{l,1},1); %SBP
            b4=A4(end);
            
            R=corr(test1,PP1{l,1},'type','Spearman');
            
            R1=corr(test1,SBP1{l,1},'type','Spearman');
            
            R2=corr(test1,DBP1{l,1},'type','Spearman');
            
            R4=corr(test1,MAP1{l,1},'type','Spearman');
            
            %------------------------------------- PLOTS WITH CORRELATION COEFS
            %AND LIN_REG COEFS
            
            lin_reg_coef=cat(2,lin_reg_coef,[A(1,1);A1(1,1);A2(1,1);A4(1,1)]);
            b_reg_coef=cat(2,b_reg_coef,[A(end);A1(end);A2(end);A4(end)]);
            
            Spearman_coef=cat(2,Spearman_coef,[R;R1;R2;R4]);
            
            rele_wind_rr=cell(1,4);
            
            feat=1;
            pp_pp=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            feat=2;
            pp_sbp=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            feat=3;
            pp_dbp=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            feat=4;
            pp_map=find(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0);
            
            rele_wind_rr{1,1}=pp_pp;
            rele_wind_rr{1,2}=pp_sbp;
            rele_wind_rr{1,3}=pp_dbp;
            rele_wind_rr{1,4}=pp_map;
            
            %--------------------------------------PLOT OF THE REGION IN ANALYSIS
            if CORR_SHOW==1
                
                figure
                ax(1)=subplot (3,1,1);
                plot(t3(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),ex_clean_ECG(locsplotecg{mm,1}(i):locsplotecg{mm,1}(i+wind_size)),'b-',t3(ECGrpeak_pos),ex_clean_ECG(ECGrpeak_pos),'g*',t3(locsplotecg{mm,1}(i:i+wind_size)),test1,'ro');
                legend('ecg','systolic peaks','HR')
                ylabel('mV')
                xlabel('Time(sec)')
                
                ax(2)=subplot (3,1,2);
                plot(t2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),SIG2(locsplotppg{mm,1}(i):locsplotppg{mm,1}(i+wind_size)),'b-',ppg80_time,ppg80_values,'r*',t2(locaisspSIG2{mm,1}),SIG2(locaisspSIG2{mm,1},mm),'g*');
                legend('ppg','20% diastolic peaks', 'Systolic peaks')
                ylabel('amplitude')
                xlabel('Time(sec)')
                
                ax(3)=subplot (3,1,3);
                plot(t1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),SIG1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),'b-',t1(locsplotabp{mm,1}(i:i+wind_size)),MAP1{l,1},'k+',t1(locsplotabp{mm,1}(i:i+wind_size)),PP1{l,1},'go',t1(locaisdiastSIG1{mm,1}),SIG1(locaisdiastSIG1{mm,1},mm),'r*',t1(locaisspSIG1{mm,1}),SIG1(locaisspSIG1{mm,1},mm),'g*');
                linkaxes(ax,'x');
                legend('bp' ,'MAP','PP','diastolic peaks','systolic peaks');
                ylabel('mmHg')
                xlabel('Time(sec)')
                
                str1 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a, R);
                title(str1)
                
                figure;
                subplot 422
                plot((test1),PP1{l,1},'go');
                xlabel('RR (sec)')  %xlabel('AC PPG (AMPLITUDE)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('PP (mmHg)')
                str1 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a, R);
                title(str1)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,PP1{l,1},txt1)
                
                subplot 424
                
                plot((test1),SBP1{l,1},'go');
                xlabel('RR (sec)')
                ylabel('SBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a1, R1);
                title(str)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,SBP1{l,1},txt1)
                
                subplot 426
                
                plot((test1),DBP1{l,1},'go');
                xlabel('RR (sec)') %xlabel('AC/DC PPG') xlabel('DC PPG')
                ylabel('DBP (mmHg)')
                str = sprintf('The slope is %f and the correlation coef. is %f',a2, R2);
                title(str)
                %             txt1 = ['Ponto ',num2str(t1(1,locsplotabp{l,1}(i:i+wind_size)))];
                %             text(test1,DBP1{l,1},txt1)
                
                subplot 428
                plot((test1),MAP1{l,1},'go');
                xlabel('RR (sec)')% xlabel('AC PPG (AMPLITUDE)') % xlabel('DC PPG')
                ylabel('MAP (mmHg)')
                str2 = sprintf('The slope is %.5f and the correlation coef. is %.2f',a4, R4);
                title(str2)
                
                ax(1)=subplot(4,2,[1,3])
                plot(t3(locsplotecg{mm,1}(i:i+wind_size)),test1)
                ylabel('RR (sec)')
                xlabel('Time(sec)')
                
                ax(2)=subplot(4,2,[5,7])
                plot(t1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)),SIG1(locsplotabp{mm,1}(i):locsplotabp{mm,1}(i+wind_size)))
                legend('BP');
                ylabel('mmHg')
                xlabel('Time(sec)')
                linkaxes(ax,'x');
                
                pause
                
            end
        end
        RRe=Spearman_coef;
        
    end
    
    %     lin_reg=array2table([lin_reg_coef],'RowNames',{'A';'A1';'A2'; 'A4'})
    inicio=pos_beg;
    final=pos_end;
    %SPEARMAN COEF
    
    wait=convsize+shift*(wind_size-1);
    update=shift*CORR_shift;
    if plotar==1
        
        figure
        ax(1)=subplot (9,2,[1 6]);
        stairs(inicio,Spearman_coef(1,:),'r-','LineWidth',2);hold on;
        stairs(inicio,Spearman_coef(2,:),'g-','LineWidth',2);hold on;
        stairs(inicio,Spearman_coef(3,:),'b-','LineWidth',2);hold on;
        stairs(inicio,Spearman_coef(4,:),'k-','LineWidth',2);hold on;
        set(gca,'fontsize',15)
        ylim([int1 int2])
        ylabel('SCC')
        title([ 'Window size ',num2str((wait-rem(wait,60))/60), 'min',num2str(rem(wait,60)),'sec.'],'fontsize',15)
        
        for m=1:numel(Spearman_coef(1,:))
            str=int2str(m);
            text(final(m),Spearman_coef(1,m),str);
            plot([final(m) final(m)],[min(Spearman_coef(2,:)) max(Spearman_coef(2,:))],'m--'); hold on;
        end
        
        legend('PP','MAP','DBP','SBP');
        
        %      BP FEATURES
        ax(2)=subplot (9,2,[11 12]);
        plot(t1(locsplotabp{l,1})./60,PP{l,1},'r*-');
        set(gca,'fontsize',15)
        
        ylabel('mmHg')
        
        %         for m=1:numel(Spearman_coef(1,:))
        %             str=int2str(m);
        %             str2={str,['CV= ' num2str(cv_PP1{l,1}{m+1,1},'%.1g')],['\chi2= ' num2str(chi_PP1{l,1}{m+1,1})]};
        %             text(final(m),mean(PP{l,1}),str2);
        %             hold on;
        %         end
        legend('PP');
        
        ax(3)=subplot (9,2,[13 14]);
        
        plot(t1(locsplotabp{l,1})./60,MAP{l,1},'g*-');
        set(gca,'fontsize',15)
        
        ylabel('mmHg')
        %
        %           for m=1:numel(Spearman_coef(1,:))
        %             str=int2str(m);
        %             str2={str,['CV= ' num2str(cv_MAP1{l,1}{m+1,1},'%.1g')],['\chi2= ' num2str(chi_MAP1{l,1}{m+1,1})]};
        %             text(final(m),mean(MAP{l,1}),str2);
        %             hold on;
        %           end
        legend('MAP');
        
        
        ax(4)=subplot (9,2,[15 16]);
        
        plot(t1(locsplotabp{l,1})./60,DBP{l,1},'b*-');
        set(gca,'fontsize',15)
        
        ylabel('mmHg')
        
        %         for m=1:numel(Spearman_coef(1,:))
        %             str=int2str(m);
        %             str2={str,['CV= ' num2str(cv_DBP1{l,1}{m+1,1},'%.1g')],['\chi2= ' num2str(chi_DBP1{l,1}{m+1,1})]};
        %             text(final(m),mean(DBP{l,1}),str2);
        %             hold on;
        %         end
        legend('DBP');
        
        ax(5)=subplot (9,2,[17 18]);
        
        plot(t1(locsplotabp{l,1})./60,SBP{l,1},'k*-');
        set(gca,'fontsize',15)
        xlabel('Time(min)')
        ylabel('mmHg')
        %         for m=1:numel(Spearman_coef(1,:))
        %             str=int2str(m);
        %             str2={str,['CV= ' num2str(cv_SBP1{l,1}{m+1,1},'%.1g')],['\chi2= ' num2str(chi_SBP1{l,1}{m+1,1})]};
        %             text(final(m),mean(SBP{l,1}),str2);
        %             hold on;
        %         end
        legend('SBP');
        
        ax(6)=subplot (9,2,[7 8]);
        stairs(inicio,lin_reg_coef(1,:),'r-','LineWidth',2);hold on;
        stairs(inicio,lin_reg_coef(2,:),'g-','LineWidth',2);hold on;
        stairs(inicio,lin_reg_coef(3,:),'b-','LineWidth',2);hold on;
        stairs(inicio,lin_reg_coef(4,:),'k-','LineWidth',2);hold on;
        set(gca,'fontsize',15)
        %plot(final,lin_reg_coef(1,:),'ro-',final,lin_reg_coef(4,:),'g*-',final,lin_reg_coef(3,:),'b+-',final,lin_reg_coef(2,:),'ks-');
        ylabel('LRC (mmHg/a.u.)')
        %
        %         for m=1:numel(Spearman_coef(1,:))
        %             str=int2str(m);
        %             plot([inicio(m) inicio(m)],[min(lin_reg_coef(2,:)) max(lin_reg_coef(2,:))],'m--','LineWidth',2); hold on;
        %             text(inicio(m),Spearman_coef(1,m),str);
        %         end
        
        
        legend('PP','MAP','DBP','SBP');%,'begin');
        
        ax(7)=subplot (9,2,[9 10]);
        
        plot(t3(locsplotecg{l,1})./60,test,'m*-');hold on;
        set(gca,'fontsize',15)
        
        ylabel([dados])  % mNp
        
        
        %         for m=1:numel(Spearman_coef(1,:))
        %             str=int2str(m);
        %             plot([inicio(m) inicio(m)],[min(test) max(test)],'m--'); hold on;
        %             text(inicio(m),mean(test),str);
        %         end
        %         legend(dados,'begin')
        
        %     %    SIGNALS
        figure
        ax(8)=subplot (3,1,1);
        plot(t2./60,SIG2,'r-'); %,t2(as{mm,1})./60,ppg_train(as{mm,1},mm),'yo',t2(bs{mm,1})./60,ppg_train(bs{mm,1},mm),'bo');%, t3./3600,ex_clean_ECG,'r-')
        legend('PPG AC') %,'a','b')%,'ECG');
        ylabel('PPG AC (mNp)')
        xlabel('Time(min)')
        
        ax(9)=subplot (3,1,2);
        plot(t1./60,SIG1);
        legend('BP');
        ylabel('BP (mmHg)')
        xlabel('Time(min)')
        
        ax(10)=subplot (3,1,3);
        plot(t3./60,ex_clean_ECG);
        legend('ECG');
        ylabel('ECG (mV)')
        xlabel('Time(min)')
        
        linkaxes(ax,'x');
    end
    %     convsize,shift,CORR_shift,wind_size
    
    Spearman_coef_tab=array2table([mean(Spearman_coef,2),mean(Spearman_coef,2)],'RowNames',{'PP';'SBP';'DBP';'MAP'},'VariableNames',{'Mean1','varcoef'})
    reg_coef_tab=array2table([mean(lin_reg_coef,2),std(lin_reg_coef,0,2)./mean(lin_reg_coef,2)*100],'RowNames',{'PP';'SBP';'DBP';'MAP'},'VariableNames',{'Mean1','varcoef'});
    %         figure;boxplot([Spearman_coef(1,:)',Spearman_coef(2,:)',Spearman_coef(3,:)',Spearman_coef(4,:)'],'Labels',{'PP','SBP','DBP','MAP'},'datalim',[-1 -0.6]); title('Spearman coefficient');
    %         title([dados ' DATA collected during ',num2str((wait-rem(wait,60))/60), 'min',num2str(rem(wait,60)),'sec (each window). Updates each ', num2str(update), ' sec'], 'FontSize', 17)
end

%---------------------------------------------------------------------- plot final results
%bar graph
% figure;
% bar(1:length(Spearman_coef(1,:)),Spearman_coef(1,:))

% box plot
if bplot==1
    %     figure;
    %     boxplot([Spearman_coef(1,:)',Spearman_coef(2,:)',Spearman_coef(3,:)',Spearman_coef(4,:)'],'Labels',{'PP','SBP','DBP','MAP'}); title('Spearman coefficient');
    bxpp=[bxpp; Spearman_coef(1,:)];
    bxmap=[bxmap; Spearman_coef(4,:)];
    bxsbp=[bxsbp; Spearman_coef(2,:)];
    bxdbp=[bxdbp; Spearman_coef(3,:)];
    
    iqrpp=[iqrpp; iqr(Spearman_coef(1,:))];
    iqrsbp=[iqrsbp; iqr(Spearman_coef(2,:))];
    iqrdbp=[iqrdbp; iqr(Spearman_coef(3,:))];
    iqrmap=[iqrmap; iqr(Spearman_coef(4,:))];
    
    mdnpp=[mdnpp; median(Spearman_coef(1,:))];
    mdnsbp=[mdnsbp; median(Spearman_coef(2,:))];
    mdndbp=[mdndbp; median(Spearman_coef(3,:))];
    mdnmap=[mdnmap; median(Spearman_coef(4,:))];
    
end

% -------------------------------------------------check correlation with lin.regression
results=[];

for feat=1:4
    
    %increase(pos.lincoef) and pos.corr
    pp=sum(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)>0)./length(lin_reg_coef);
    
    %increase(pos.lincoef) and neg.corr
    np=sum(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)>0)./length(lin_reg_coef);
    
    %decrease(neg.lincoef) and pos.corr
    pn=sum(Spearman_coef(feat,:)>threshold & lin_reg_coef(feat,:)<0)./length(lin_reg_coef);
    
    %decrease(neg.lincoef) and neg.corr
    nn=sum(Spearman_coef(feat,:)<-threshold & lin_reg_coef(feat,:)<0)./length(lin_reg_coef);
    
    results=[results; pp*100 np*100; pn*100 nn*100];
    
end
results'

%% load data
rmse_group=[];
scc_group=[];
npp=[];
nsbp=[];
ndbp=[];
nmap=[];
ppba=[];
sbpba=[];
dbpba=[];
mapba=[];
calib=0;
ppaami=[];
sbpaami=[];
dbpaami=[];
mapaami=[];

signal='s_5_1_n.mat';
s_5_1_n=load(signal);
signal='s_1_1_n.mat';
s_1_1_n=load(signal);
signal='s_3_1_n.mat';
s_3_1_n=load(signal);

signal='s_5_0_1_n.mat';
s_5_0_1_n=load(signal);
signal='s_5_0_2_n.mat';
s_5_0_2_n=load(signal);
signal='s_5_0_3_n.mat';
s_5_0_3_n=load(signal);
signal='s_5_0_4_n.mat';
s_5_0_4_n=load(signal);

signal='s_1_0_1_n.mat';
s_1_0_1_n=load(signal);
signal='s_1_0_2_n.mat';
s_1_0_2_n=load(signal);
signal='s_1_0_3_n.mat';
s_1_0_3_n=load(signal);
signal='s_1_0_4_n.mat';
s_1_0_4_n=load(signal);

signal='s_3_0_1_n.mat';
s_3_0_1_n=load(signal);
signal='s_3_0_2_n.mat';
s_3_0_2_n=load(signal);
signal='s_3_0_3_n.mat';
s_3_0_3_n=load(signal);
signal='s_3_0_4_n.mat';
s_3_0_4_n=load(signal);

signal='s_5_1.mat';
s_5_1=load(signal);
signal='s_1_1.mat';
s_1_1=load(signal);
signal='s_3_1.mat';
s_3_1=load(signal);

signal='s_5_0_1.mat';
s_5_0_1=load(signal);
signal='s_5_0_2.mat';
s_5_0_2=load(signal);
signal='s_5_0_3.mat';
s_5_0_3=load(signal);
signal='s_5_0_4.mat';
s_5_0_4=load(signal);

signal='s_1_0_1.mat';
s_1_0_1=load(signal);
signal='s_1_0_2.mat';
s_1_0_2=load(signal);
signal='s_1_0_3.mat';
s_1_0_3=load(signal);
signal='s_1_0_4.mat';
s_1_0_4=load(signal);

signal='s_3_0_1.mat';
s_3_0_1=load(signal);
signal='s_3_0_2.mat';
s_3_0_2=load(signal);
signal='s_3_0_3.mat';
s_3_0_3=load(signal);
signal='s_3_0_4.mat';
s_3_0_4=load(signal);

%% Multiple linear regression (real vs model)
mm=1; l=1;
opt=4;
% with calibration
% 1: 2 min parameters (left out patient)
% 2: 2 min parameters (left out patient) (no PAT1 or PAT2)
% 3: 2 min parameters (left out patient) (no AC/DC or DC)
% 7: 2 min parameters (left out patient) (no PAT2)
wind_size=23;
initi=5*wind_size; %(PAT1:1*wind_size; PAT3:2*wind_size  PAT5:5*wind_size )
initf=6*wind_size; %1*23
%testing
ff=0.99;
BA=1;
load('param_111_35_n.mat')
% 4: final hospital implementation
% 5: final hospital implementation (no PAT1 or PAT2)
% 6: final hospital implementation (no ac/dc, dc)
% 8: final hospital implementation (no PAT2)
info=[s_1_1_n.together(:,1:9)];
pp=[s_1_1_n.together(:,10)];
sbp=[s_1_1_n.together(:,11)];
dbp=[s_1_1_n.together(:,12)];
map=[s_1_1_n.together(:,13)];

if opt==1
    
    X = [info];
    
    obj1 = recursiveLS(10,PPf,'ForgettingFactor',ff);
    obj2 = recursiveLS(10,SBPf,'ForgettingFactor',ff);
    obj3 = recursiveLS(10,DBPf,'ForgettingFactor',ff);
    obj4 = recursiveLS(10,MAPf,'ForgettingFactor',ff);
    
    
    for i = initi:initf %runs first 2 minutes to adapt parameters
        H = [1 X(i,:)];
        [PPcal, ~] = step(obj1,pp(i),H);
        [SBPcal, ~] = step(obj2,sbp(i),H);
        [DBPcal, ~] = step(obj3,dbp(i),H);
        [MAPcal, ~] = step(obj4,map(i),H);
    end
    calib=1;
    
elseif opt==2
    
    k=initf;
    X = [info(:,3:end)];
    
    obj1 = recursiveLS(8,PPf,'ForgettingFactor',ff);
    obj2 = recursiveLS(8,SBPf,'ForgettingFactor',ff);
    obj3 = recursiveLS(8,DBPf,'ForgettingFactor',ff);
    obj4 = recursiveLS(8,MAPf,'ForgettingFactor',ff);
    
    
    for i = initi:initf %runs first 2 minutes to adapt parameters
        H = [1 X(i,:)];
        [PPcal, ~] = step(obj1,pp(i),H);
        [SBPcal, ~] = step(obj2,sbp(i),H);
        [DBPcal, ~] = step(obj3,dbp(i),H);
        [MAPcal, ~] = step(obj4,map(i),H);
    end
        calib=1;
elseif opt==3
    
    k=initf;
    X = [info(:,1:2) info(:,4) info(:,6:9)];
    
    obj1 = recursiveLS(8,PPf,'ForgettingFactor',ff);
    obj2 = recursiveLS(8,SBPf,'ForgettingFactor',ff);
    obj3 = recursiveLS(8,DBPf,'ForgettingFactor',ff);
    obj4 = recursiveLS(8,MAPf,'ForgettingFactor',ff);
    
    
    for i = initi:initf %runs first 2 minutes to adapt parameters
        H = [1 X(i,:)];
        [PPcal, ~] = step(obj1,pp(i),H);
        [SBPcal, ~] = step(obj2,sbp(i),H);
        [DBPcal, ~] = step(obj3,dbp(i),H);
        [MAPcal, ~] = step(obj4,map(i),H);
    end
        calib=1;
elseif opt==7
    
    k=initf;
    X = [info(:,1) info(:,3:9)];
    
    obj1 = recursiveLS(9,PPf,'ForgettingFactor',ff);
    obj2 = recursiveLS(9,SBPf,'ForgettingFactor',ff);
    obj3 = recursiveLS(9,DBPf,'ForgettingFactor',ff);
    obj4 = recursiveLS(9,MAPf,'ForgettingFactor',ff);
    
    
    for i = initi:initf %runs first 2 minutes to adapt parameters
        H = [1 X(i,:)];
        [PPcal, ~] = step(obj1,pp(i),H);
        [SBPcal, ~] = step(obj2,sbp(i),H);
        [DBPcal, ~] = step(obj3,dbp(i),H);
        [MAPcal, ~] = step(obj4,map(i),H);
    end
        calib=1;
elseif opt==4
    
    X = [ones(size(info,1),1) info];
    
    k=1;
    
   if calib==0
        PPPARAM = PPf;
        
        SBPPARAM = SBPf;
        
        DBPPARAM = DBPf;
        
        MAPPARAM = MAPf;
    elseif calib==1
        PPPARAM = PPcal;
        
        SBPPARAM = SBPcal;
        
        DBPPARAM = DBPcal;
        
        MAPPARAM = MAPcal;
   end
    
    %----------------------------------------------------error calculation
    errmap=[];
    errpp=[];
    errdbp=[];
    errsbp=[];
    
    y=pp;
    PPFIT=PPPARAM(1) + PPPARAM(2)*X(k:end,2) + PPPARAM(3)*X(k:end,3) + PPPARAM(4)*X(k:end,4) + PPPARAM(5)*X(k:end,5) + PPPARAM(6)*X(k:end,6) + PPPARAM(7)*X(k:end,7) + PPPARAM(8)*X(k:end,8) + PPPARAM(9)*X(k:end,9) + PPPARAM(10)*X(k:end,10);
    ppmse = immse(PPFIT,y);
    
    y=sbp;
    SBPFIT=SBPPARAM(1) + SBPPARAM(2)*X(k:end,2) + SBPPARAM(3)*X(k:end,3) + SBPPARAM(4)*X(k:end,4) + SBPPARAM(5)*X(k:end,5) + SBPPARAM(6)*X(k:end,6) + SBPPARAM(7)*X(k:end,7) + SBPPARAM(8)*X(k:end,8) + SBPPARAM(9)*X(k:end,9) + SBPPARAM(10)*X(k:end,10);
    sbpmse = immse(SBPFIT,y);
    
    y=dbp;
    DBPFIT=DBPPARAM(1) + DBPPARAM(2)*X(k:end,2) + DBPPARAM(3)*X(k:end,3) + DBPPARAM(4)*X(k:end,4) + DBPPARAM(5)*X(k:end,5) + DBPPARAM(6)*X(k:end,6) + DBPPARAM(7)*X(k:end,7) + DBPPARAM(8)*X(k:end,8) + DBPPARAM(9)*X(k:end,9) + DBPPARAM(10)*X(k:end,10);
    dbpmse = immse(DBPFIT,y);
    
    y=map;
    MAPFIT=MAPPARAM(1) + MAPPARAM(2)*X(k:end,2) + MAPPARAM(3)*X(k:end,3) + MAPPARAM(4)*X(k:end,4) + MAPPARAM(5)*X(k:end,5) + MAPPARAM(6)*X(k:end,6) + MAPPARAM(7)*X(k:end,7) + MAPPARAM(8)*X(k:end,8) + MAPPARAM(9)*X(k:end,9) + MAPPARAM(10)*X(k:end,10);
    mapmse = immse(MAPFIT,y);
    
    errors=[sqrt(ppmse) sqrt(sbpmse) sqrt(dbpmse) sqrt(mapmse)]
    rmse_group=[rmse_group errors'];
    
    corrs=[corr(PPFIT,pp,'type','Spearman') corr(SBPFIT,sbp,'type','Spearman') corr(DBPFIT,dbp,'type','Spearman') corr(MAPFIT,map,'type','Spearman') ]
    scc_group=[scc_group corrs'];
            calib=0;

%     %real vs fit plot
%     f=figure;
%     ax(1)=subplot (4,1,1);
%     plot(t1(locsplotabp{mm,1}(k:end))./60,PPFIT,'r-',t1(locsplotabp{mm,1})./60,PP{l,1},'g-')
%     ylabel('mmHg')
%     xlabel('Time(min)')
%     legend('PP fit','PP real')
%     
%     ax(2)=subplot (4,1,2);
%     plot(t1(locsplotabp{mm,1}(k:end))./60,SBPFIT,'r-',t1(locsplotabp{mm,1})./60,SBP{l,1},'g-')
%     ylabel('mmHg')
%     xlabel('Time(min)')
%     legend('SBP fit','SBP real')
%     
%     ax(3)=subplot (4,1,3);
%     plot(t1(locsplotabp{mm,1}(k:end))./60,DBPFIT,'r-',t1(locsplotabp{mm,1})./60,DBP{l,1},'g-')
%     ylabel('mmHg')
%     xlabel('Time(min)')
%     legend('DBP fit','DBP real')
%     
%     ax(4)=subplot (4,1,4);
%     plot(t1(locsplotabp{mm,1}(k:end))./60,MAPFIT,'r-',t1(locsplotabp{mm,1})./60,MAP{l,1},'g-')
%     ylabel('mmHg')
%     xlabel('Time(min)')
%     legend('MAP fit','MAP real')
%     
%     set(f,'name','Causal system','numbertitle','off')
%     
%     linkaxes(ax,'x');
%   
    
elseif opt==8
    
    X = [ones(size(info,1),1) info(:,1) info(:,3:9)];
    
    k=1;
    
   if calib==0
        
        PPPARAM = PPf;
        
        SBPPARAM = SBPf;
        
        DBPPARAM = DBPf;
        
        MAPPARAM = MAPf;
    elseif calib==1
    
        PPPARAM = PPcal;
        
        SBPPARAM = SBPcal;
        
        DBPPARAM = DBPcal;
        
        MAPPARAM = MAPcal;
    end
    %----------------------------------------------------error calculation
    errmap=[];
    errpp=[];
    errdbp=[];
    errsbp=[];
    
    y=pp;
    PPFIT=PPPARAM(1) + PPPARAM(2)*X(k:end,2) + PPPARAM(3)*X(k:end,3) + PPPARAM(4)*X(k:end,4) + PPPARAM(5)*X(k:end,5) + PPPARAM(6)*X(k:end,6) + PPPARAM(7)*X(k:end,7) + PPPARAM(8)*X(k:end,8) + PPPARAM(9)*X(k:end,9);
    ppmse = immse(PPFIT,y);
    
    y=sbp;
    SBPFIT=SBPPARAM(1) + SBPPARAM(2)*X(k:end,2) + SBPPARAM(3)*X(k:end,3) + SBPPARAM(4)*X(k:end,4) + SBPPARAM(5)*X(k:end,5) + SBPPARAM(6)*X(k:end,6) + SBPPARAM(7)*X(k:end,7) + SBPPARAM(8)*X(k:end,8) + SBPPARAM(9)*X(k:end,9);
    sbpmse = immse(SBPFIT,y);
    
    y=dbp;
    DBPFIT=DBPPARAM(1) + DBPPARAM(2)*X(k:end,2) + DBPPARAM(3)*X(k:end,3) + DBPPARAM(4)*X(k:end,4) + DBPPARAM(5)*X(k:end,5) + DBPPARAM(6)*X(k:end,6) + DBPPARAM(7)*X(k:end,7) + DBPPARAM(8)*X(k:end,8) + DBPPARAM(9)*X(k:end,9);
    dbpmse = immse(DBPFIT,y);
    
    y=map;
    MAPFIT=MAPPARAM(1) + MAPPARAM(2)*X(k:end,2) + MAPPARAM(3)*X(k:end,3) + MAPPARAM(4)*X(k:end,4) + MAPPARAM(5)*X(k:end,5) + MAPPARAM(6)*X(k:end,6) + MAPPARAM(7)*X(k:end,7) + MAPPARAM(8)*X(k:end,8) + MAPPARAM(9)*X(k:end,9);
    mapmse = immse(MAPFIT,y);
    
    errors=[sqrt(ppmse) sqrt(sbpmse) sqrt(dbpmse) sqrt(mapmse)]
    rmse_group=[rmse_group errors'];
    
    corrs=[corr(PPFIT,pp,'type','Spearman') corr(SBPFIT,sbp,'type','Spearman') corr(DBPFIT,dbp,'type','Spearman') corr(MAPFIT,map,'type','Spearman') ]
    scc_group=[scc_group corrs'];
            calib=0;

%     %real vs fit plot
%     f=figure;
%     ax(1)=subplot (4,1,1);
%     plot(t1(locsplotabp{mm,1}(k:end))./60,PPFIT,'r-',t1(locsplotabp{mm,1})./60,PP{l,1},'g-')
%     ylabel('mmHg')
%     xlabel('Time(min)')
%     legend('PP fit','PP real')
%     
%     ax(2)=subplot (4,1,2);
%     plot(t1(locsplotabp{mm,1}(k:end))./60,SBPFIT,'r-',t1(locsplotabp{mm,1})./60,SBP{l,1},'g-')
%     ylabel('mmHg')
%     xlabel('Time(min)')
%     legend('SBP fit','SBP real')
%     
%     ax(3)=subplot (4,1,3);
%     plot(t1(locsplotabp{mm,1}(k:end))./60,DBPFIT,'r-',t1(locsplotabp{mm,1})./60,DBP{l,1},'g-')
%     ylabel('mmHg')
%     xlabel('Time(min)')
%     legend('DBP fit','DBP real')
%     
%     ax(4)=subplot (4,1,4);
%     plot(t1(locsplotabp{mm,1}(k:end))./60,MAPFIT,'r-',t1(locsplotabp{mm,1})./60,MAP{l,1},'g-')
%     ylabel('mmHg')
%     xlabel('Time(min)')
%     legend('MAP fit','MAP real')
%     
%     set(f,'name','Causal system','numbertitle','off')
%     
%     linkaxes(ax,'x');
elseif opt==5

    X = [ones(size(info,1),1) info(:,3:9)];
    
    k=1;
   
    if calib==0
        PPPARAM = PPf;
        
        SBPPARAM = SBPf;
        
        DBPPARAM = DBPf;
        
        MAPPARAM = MAPf;
    elseif calib==1
        PPPARAM = PPcal;
        
        SBPPARAM = SBPcal;
        
        DBPPARAM = DBPcal;
        
        MAPPARAM = MAPcal;
   end
    
    %----------------------------------------------------error calculation
    errmap=[];
    errpp=[];
    errdbp=[];
    errsbp=[];
    
    y=pp;
    PPFIT=PPPARAM(1) + PPPARAM(2)*X(k:end,2) + PPPARAM(3)*X(k:end,3) + PPPARAM(4)*X(k:end,4) + PPPARAM(5)*X(k:end,5) + PPPARAM(6)*X(k:end,6) + PPPARAM(7)*X(k:end,7)+ PPPARAM(8)*X(k:end,8);
    ppmse = immse(PPFIT,y);
    
    y=sbp;
    SBPFIT=SBPPARAM(1) + SBPPARAM(2)*X(k:end,2) + SBPPARAM(3)*X(k:end,3) + SBPPARAM(4)*X(k:end,4) + SBPPARAM(5)*X(k:end,5) + SBPPARAM(6)*X(k:end,6) + SBPPARAM(7)*X(k:end,7)+ SBPPARAM(8)*X(k:end,8);
    sbpmse = immse(SBPFIT,y);
    
    y=dbp;
    DBPFIT=DBPPARAM(1) + DBPPARAM(2)*X(k:end,2) + DBPPARAM(3)*X(k:end,3) + DBPPARAM(4)*X(k:end,4) + DBPPARAM(5)*X(k:end,5) + DBPPARAM(6)*X(k:end,6) + DBPPARAM(7)*X(k:end,7)+ DBPPARAM(8)*X(k:end,8);
    dbpmse = immse(DBPFIT,y);
    
    y=map;
    MAPFIT=MAPPARAM(1) + MAPPARAM(2)*X(k:end,2) + MAPPARAM(3)*X(k:end,3) + MAPPARAM(4)*X(k:end,4) + MAPPARAM(5)*X(k:end,5) + MAPPARAM(6)*X(k:end,6) + MAPPARAM(7)*X(k:end,7)+ MAPPARAM(8)*X(k:end,8);
    mapmse = immse(MAPFIT,y);
    
    errors=[sqrt(ppmse) sqrt(sbpmse) sqrt(dbpmse) sqrt(mapmse)]
    rmse_group=[rmse_group errors'];

    corrs=[corr(PPFIT,pp,'type','Spearman') corr(SBPFIT,sbp,'type','Spearman') corr(DBPFIT,dbp,'type','Spearman') corr(MAPFIT,map,'type','Spearman') ]
    scc_group=[scc_group corrs'];
        calib=0;

    %real vs fit plot
%     f=figure;
%     ax(1)=subplot (4,1,1);
%     plot(t1(locsplotabp{mm,1}(k:end))./60,PPFIT,'r-',t1(locsplotabp{mm,1})./60,PP{l,1},'g-')
%     ylabel('mmHg')
%     xlabel('Time(min)')
%     legend('PP fit','PP real')
%     
%     ax(2)=subplot (4,1,2);
%     plot(t1(locsplotabp{mm,1}(k:end))./60,SBPFIT,'r-',t1(locsplotabp{mm,1})./60,SBP{l,1},'g-')
%     ylabel('mmHg')
%     xlabel('Time(min)')
%     legend('SBP fit','SBP real')
%     
%     ax(3)=subplot (4,1,3);
%     plot(t1(locsplotabp{mm,1}(k:end))./60,DBPFIT,'r-',t1(locsplotabp{mm,1})./60,DBP{l,1},'g-')
%     ylabel('mmHg')
%     xlabel('Time(min)')
%     legend('DBP fit','DBP real')
%     
%     ax(4)=subplot (4,1,4);
%     plot(t1(locsplotabp{mm,1}(k:end))./60,MAPFIT,'r-',t1(locsplotabp{mm,1})./60,MAP{l,1},'g-')
%     ylabel('mmHg')
%     xlabel('Time(min)')
%     legend('MAP fit','MAP real')
%     set(f,'name','Causal system','numbertitle','off')
%     linkaxes(ax,'x');
%     
    %repeatibility
    %     rep=2.77*[std(PP{mm,1}) std(SBP{mm,1}) std(DBP{mm,1}) std(MAP{mm,1})]
    %     repfit=2.77*[std(DBPFIT) std(SBPFIT) std(PPFIT) std(MAPFIT)]
elseif opt==6
    
    X = [ones(size(info,1),1) info(:,1:2) info(:,4) info(:,6:9)];
    
    k=1;
    if calib==0
        PPPARAM = PPf;
        
        SBPPARAM = SBPf;
        
        DBPPARAM = DBPf;
        
        MAPPARAM = MAPf;
    elseif calib==1
        PPPARAM = PPcal;
        
        SBPPARAM = SBPcal;
        
        DBPPARAM = DBPcal;
        
        MAPPARAM = MAPcal;
    end
    %----------------------------------------------------error calculation
    errmap=[];
    errpp=[];
    errdbp=[];
    errsbp=[];
    
    y=pp(k:end);
    PPFIT=PPPARAM(1) + PPPARAM(2)*X(k:end,2) + PPPARAM(3)*X(k:end,3) + PPPARAM(4)*X(k:end,4) + PPPARAM(5)*X(k:end,5) + PPPARAM(6)*X(k:end,6) + PPPARAM(7)*X(k:end,7) + PPPARAM(8)*X(k:end,8);
    ppmse = immse(PPFIT,y);
    
    y=sbp(k:end);
    SBPFIT=SBPPARAM(1) + SBPPARAM(2)*X(k:end,2) + SBPPARAM(3)*X(k:end,3) + SBPPARAM(4)*X(k:end,4) + SBPPARAM(5)*X(k:end,5) + SBPPARAM(6)*X(k:end,6) + SBPPARAM(7)*X(k:end,7) + SBPPARAM(8)*X(k:end,8);
    sbpmse = immse(SBPFIT,y);
    
    y=dbp(k:end);
    DBPFIT=DBPPARAM(1) + DBPPARAM(2)*X(k:end,2) + DBPPARAM(3)*X(k:end,3) + DBPPARAM(4)*X(k:end,4) + DBPPARAM(5)*X(k:end,5) + DBPPARAM(6)*X(k:end,6) + DBPPARAM(7)*X(k:end,7) + DBPPARAM(8)*X(k:end,8);
    dbpmse = immse(DBPFIT,y);
    
    y=map(k:end);
    MAPFIT=MAPPARAM(1) + MAPPARAM(2)*X(k:end,2) + MAPPARAM(3)*X(k:end,3) + MAPPARAM(4)*X(k:end,4) + MAPPARAM(5)*X(k:end,5) + MAPPARAM(6)*X(k:end,6) + MAPPARAM(7)*X(k:end,7) + MAPPARAM(8)*X(k:end,8);
    mapmse = immse(MAPFIT,y);
    
    errors=[sqrt(ppmse) sqrt(sbpmse) sqrt(dbpmse) sqrt(mapmse)]
    rmse_group=[rmse_group errors'];
    corrs=[corr(PPFIT,pp,'type','Spearman') corr(SBPFIT,sbp,'type','Spearman') corr(DBPFIT,dbp,'type','Spearman') corr(MAPFIT,map,'type','Spearman') ]
    scc_group=[scc_group corrs'];
           calib=0;
 
    %real vs fit plot
%     f=figure;
%     ax(1)=subplot (4,1,1);
%     plot(t1(locsplotabp{mm,1}(k:end))./60,PPFIT,'r-',t1(locsplotabp{mm,1})./60,pp,'g-')
%     ylabel('mmHg')
%     xlabel('Time(min)')
%     legend('PP fit','PP real')
%     
%     ax(2)=subplot (4,1,2);
%     plot(t1(locsplotabp{mm,1}(k:end))./60,SBPFIT,'r-',t1(locsplotabp{mm,1})./60,sbp,'g-')
%     ylabel('mmHg')
%     xlabel('Time(min)')
%     legend('SBP fit','SBP real')
%     
%     ax(3)=subplot (4,1,3);
%     plot(t1(locsplotabp{mm,1}(k:end))./60,DBPFIT,'r-',t1(locsplotabp{mm,1})./60,dbp,'g-')
%     ylabel('mmHg')
%     xlabel('Time(min)')
%     legend('DBP fit','DBP real')
%     
%     ax(4)=subplot (4,1,4);
%     plot(t1(locsplotabp{mm,1}(k:end))./60,MAPFIT,'r-',t1(locsplotabp{mm,1})./60,map,'g-')
%     ylabel('mmHg')
%     xlabel('Time(min)')
%     legend('MAP fit','MAP real')
%     set(f,'name','Causal system','numbertitle','off')
%     linkaxes(ax,'x');
    
    %repeatibility
    %     rep=2.77*[std(PP{mm,1}) std(SBP{mm,1}) std(DBP{mm,1}) std(MAP{mm,1})]
    %     repfit=2.77*[std(DBPFIT) std(SBPFIT) std(PPFIT) std(MAPFIT)]
end
       
       npp=[npp; PPFIT-mean(PPFIT-pp)];
       nsbp=[nsbp; SBPFIT-mean(SBPFIT-sbp)];
       ndbp=[ndbp; DBPFIT-mean(DBPFIT-dbp)];
       nmap=[nmap; MAPFIT-mean(MAPFIT-map)];

       ppba=[ppba; pp];
       sbpba=[sbpba; sbp];
       dbpba=[dbpba; dbp];
       mapba=[mapba; map];
       
       ppaami=[ppaami; PPFIT-pp];
       sbpaami=[sbpaami; SBPFIT-sbp];
       dbpaami=[dbpaami; DBPFIT-dbp];
       mapaami=[mapaami; MAPFIT-map];
       
%%       
if BA==1
%        BlandAltman(ppba, npp,{'PP','pp'},'corrInfo',{'eq','RMSE','rho'},'baInfo',{'IQR'},'data1Mode','Truth');
%     
%        BlandAltman(sbpba, nsbp, {'SBP','sbp'},'corrInfo',{'eq','RMSE','rho'},'baInfo',{'IQR'},'data1Mode','Truth');
%    
%        BlandAltman(dbpba, ndbp, {'DBP','dbp'},'corrInfo',{'eq','RMSE','rho'},'baInfo',{'IQR'},'data1Mode','Truth');
%     
%        BlandAltman(mapba, nmap,{'MAP','map'},'corrInfo',{'eq','RMSE','rho'},'baInfo',{'IQR'},'data1Mode','Truth');

       mediapp=mean(npp-ppba)
       sdpp=std(npp-ppba)
       
       mediasbp=mean(nsbp-sbpba)
       sdsbp=std(nsbp-sbpba)
       mediadbp=mean(ndbp-dbpba)
       sddbp=std(ndbp-dbpba)
       mediamap=mean(nmap-mapba)
       sdmap=std(nmap-mapba)
%        
%        mediapp=mean(ppaami)
%        sdpp=std(ppaami)
%        
%        mediasbp=mean(sbpaami)
%        sdsbp=std(sbpaami)
%        mediadbp=mean(dbpaami)
%        sddbp=std(dbpaami)
%        mediamap=mean(mapaami)
%        sdmap=std(mapaami)
%        
       pp_5=100*sum(abs(npp-ppba)<5)/length(npp-ppba);  
       pp_10=100*sum(abs(npp-ppba)<10)/length(npp-ppba);  
       pp_15=100*sum(abs(npp-ppba)<15)/length(npp-ppba);  
       
       sbp_5=100*sum(abs(nsbp-sbpba)<5)/length(nsbp-sbpba);  
       sbp_10=100*sum(abs(nsbp-sbpba)<10)/length(nsbp-sbpba);  
       sbp_15=100*sum(abs(nsbp-sbpba)<15)/length(nsbp-sbpba);  
       
       dbp_5=100*sum(abs(ndbp-dbpba)<5)/length(ndbp-dbpba);  
       dbp_10=100*sum(abs(ndbp-dbpba)<10)/length(ndbp-dbpba);  
       dbp_15=100*sum(abs(ndbp-dbpba)<15)/length(ndbp-dbpba);  
       
       map_5=100*sum(abs(nmap-mapba)<5)/length(nmap-mapba);  
       map_10=100*sum(abs(nmap-mapba)<10)/length(nmap-mapba);  
       map_15=100*sum(abs(nmap-mapba)<15)/length(nmap-mapba);  
       
%        tab=[pp_5 pp_10 pp_15; sbp_5 sbp_10 sbp_15; dbp_5 dbp_10 dbp_15; map_5 map_10 map_15]
end
%% save SIGNALS

entire=1;
feat=3;

%saving
% 1: save parts of the signal with SCC>70% (entire=0) and entire signal (entire=1)
% (feat=1 (PPvsAC/DC) feat=2 (SBPvsAC/DC) feat=3 (DBPvsAC/DC) feat=4 (MAPvsAC/DC)): PATIENT 1
% (feat=1 (PPvsB/A) feat=2 (SBPvsB/A) feat=3 (DBPvsAC) feat=4 (MAPvsB/A)): PATIENT 3
% (feat=1 (PPvsDC) feat=2 (SBPvsDC) feat=3 (DBPvsAC/DC) feat=4 (MAPvsDC)): PATIENT 5

signal='s_5_1_n.mat';

if entire==0
    threshold=0.7;
    wws=find(abs(Spearman_coef(feat,:))>threshold);
    
    pat=[];
    pat2=[];
    amp_ac_dc=[];
    amp_ac=[];
    amp_dc=[];
    ri=[];
    hbi=[];
    si=[];
    dn=[];
    
    pp=[];
    sbp=[];
    dbp=[];
    map=[];
    
    for k=1:length(Spearman_coef(feat,:))-1 %experimentar ate tamanho do locaisspSIG2
        
        if ( any(wws==k) )
            
            pat=[pat; PAT{mm,1}(k*23:(k+1)*23)];
            pat2=[pat2; PAT2{mm,1}(k*23:(k+1)*23)];
            
            amp_ac_dc=[amp_ac_dc; amp_AC_DC{mm,1}(k*23:(k+1)*23)];
            amp_dc=[amp_dc; amp_DC{mm,1}(k*23:(k+1)*23)];
            amp_ac=[amp_ac; amp_AC{mm,1}(k*23:(k+1)*23)];
            ri=[ri; RI{mm,1}(k*23:(k+1)*23)];
            hbi=[hbi; HBI(k*23:(k+1)*23)];
            si=[si; SI(k*23:(k+1)*23)];
            dn=[dn; DN{mm,1}(k*23:(k+1)*23)];
            
            pp=[pp; PP{mm,1}(k*23:(k+1)*23)];
            sbp=[sbp; SBP{mm,1}(k*23:(k+1)*23)];
            dbp=[dbp; DBP{mm,1}(k*23:(k+1)*23)];
            map=[map; MAP{mm,1}(k*23:(k+1)*23)];
        else
        end
    end
    
    pat1=pat;
    pat21=pat2;
    amp_ac_dc1=amp_ac_dc;
    amp_ac1=amp_ac;
    amp_dc1=amp_dc;
    ri1=ri;
    hbi1=hbi;
    si1=si;
    dn1=dn;
    
    pp1=pp;
    sbp1=sbp;
    dbp1=dbp;
    map1=map;
    
    together=[pat1,pat21,amp_ac_dc1,amp_ac1,amp_dc1,ri1,hbi1,si1,dn1,pp1,sbp1,dbp1,map1];
    save(signal,'together')
    
elseif entire==1
    
    pat1=PAT{mm,1};
    pat21=PAT2{mm,1};
    amp_ac_dc1=amp_AC_DC{mm,1};
    amp_ac1= amp_AC{mm,1};
    amp_dc1= amp_DC{mm,1};
    ri1= RI{mm,1} ;
    hbi1=HBI;
    si1=SI;
    dn1=DN{mm,1};
    
    pp1=PP{mm,1};
    sbp1=SBP{mm,1};
    dbp1=DBP{mm,1};
    map1=MAP{mm,1};
    
    together=[pat1,pat21,amp_ac_dc1,amp_ac1,amp_dc1,ri1,hbi1,si1,dn1,pp1,sbp1,dbp1,map1];
    
    save(signal,'together')
end


%% save parameters
param='param_101_13_n.mat';
%'param_xyz_ij.mat';

% ASSEMBLE INFO
info_pp=[s_1_1_n.together(:,1:2) s_1_1_n.together(:,4) s_1_1_n.together(:,6:9); s_3_1_n.together(:,1:2) s_3_1_n.together(:,4) s_3_1_n.together(:,6:9)];
info_sbp=[s_1_1_n.together(:,1:2) s_1_1_n.together(:,4) s_1_1_n.together(:,6:9); s_3_1_n.together(:,1:2) s_3_1_n.together(:,4) s_3_1_n.together(:,6:9)];
info_dbp=[s_1_1_n.together(:,1:2) s_1_1_n.together(:,4) s_1_1_n.together(:,6:9); s_3_1_n.together(:,1:2) s_3_1_n.together(:,4) s_3_1_n.together(:,6:9)];
info_map=[s_1_1_n.together(:,1:2) s_1_1_n.together(:,4) s_1_1_n.together(:,6:9); s_3_1_n.together(:,1:2) s_3_1_n.together(:,4) s_3_1_n.together(:,6:9)];
pp=[s_1_1_n.together(:,10); s_3_1_n.together(:,10)];
sbp=[s_1_1_n.together(:,11); s_3_1_n.together(:,11)];
dbp=[s_1_1_n.together(:,12); s_3_1_n.together(:,12)];
map=[s_1_1_n.together(:,13); s_3_1_n.together(:,13)];

%CALCULATION
pp_test = [ones(size(info_pp,1),1) info_pp];
sbp_test = [ones(size(info_sbp,1),1) info_sbp];
dbp_test = [ones(size(info_dbp,1),1) info_dbp];
map_test = [ones(size(info_map,1),1) info_map];

[PPf,~,~]=regress(pp,pp_test);
[SBPf,~,~]=lscov(sbp_test,sbp);
[DBPf,~,~]=lscov(dbp_test,dbp);
[MAPf,~,~]=regress(map,map_test)

save(param,'PPf','SBPf','DBPf','MAPf')

%% kruskal wallis
figure
info=[s_5_1.together(:,1:9)];
[p,tbl,stats] = kruskalwallis(info,[],'off');
c = multcompare(stats)
set(gca,'fontsize',25)

% Two group means are significantly different if their intervals are disjoint;
