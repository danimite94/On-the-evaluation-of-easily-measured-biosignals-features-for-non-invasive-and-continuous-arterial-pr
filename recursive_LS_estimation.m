load('param_111_13.mat')
mm=1;
l=1;
obj1 = recursiveLS(10,PPf,'ForgettingFactor',0.98);
obj2 = recursiveLS(10,SBPf,'ForgettingFactor',0.98);
obj3 = recursiveLS(10,DBPf,'ForgettingFactor',0.98);
obj4 = recursiveLS(10,MAPf,'ForgettingFactor',0.98);

% Set 'ForgettingFactor' less than 1 for estimating time-varying coefficients.
input=s_5_1.together(:,1:9);
output1=s_5_1.together(:,10);
output2=s_5_1.together(:,11);
output3=s_5_1.together(:,12);
output4=s_5_1.together(:,13);
out=1;
for i = 1:2*23 %runs first 2 minutes to adapt parameters
    H = [out info(i,:)];
    [PPcal, ~] = step(obj1,pp,H);
    [SBPcal, ~] = step(obj2,sbp,H);
    [DBPcal, ~] = step(obj3,dbp,H);
    [MAPcal, ~] = step(obj4,map,H);
    
end

% PPPARAM = PPf;
% 
% SBPPARAM = SBPf;
% 
% DBPPARAM = DBPf;
% 
% MAPPARAM = MAPf;

k=1;
X = [ones(size(s_5_1.together(:,1:9),1),1) s_5_1.together(:,1:9)];
errmap=[];
errpp=[];
errdbp=[];
errsbp=[];

y=output1;
PPFIT=PPPARAM(1) + PPPARAM(2)*X(k:end,2) + PPPARAM(3)*X(k:end,3) + PPPARAM(4)*X(k:end,4) + PPPARAM(5)*X(k:end,5) + PPPARAM(6)*X(k:end,6) + PPPARAM(7)*X(k:end,7) + PPPARAM(8)*X(k:end,8) + PPPARAM(9)*X(k:end,9);
ppmse = immse(PPFIT,y);

y=output2;
SBPFIT=SBPPARAM(1) + SBPPARAM(2)*X(k:end,2) + SBPPARAM(3)*X(k:end,3) + SBPPARAM(4)*X(k:end,4) + SBPPARAM(5)*X(k:end,5) + SBPPARAM(6)*X(k:end,6) + SBPPARAM(7)*X(k:end,7) + SBPPARAM(8)*X(k:end,8) + SBPPARAM(9)*X(k:end,9);
sbpmse = immse(SBPFIT,y);

y=output3;
DBPFIT=DBPPARAM(1) + DBPPARAM(2)*X(k:end,2) + DBPPARAM(3)*X(k:end,3) + DBPPARAM(4)*X(k:end,4) + DBPPARAM(5)*X(k:end,5) + DBPPARAM(6)*X(k:end,6) + DBPPARAM(7)*X(k:end,7) + DBPPARAM(8)*X(k:end,8) + DBPPARAM(9)*X(k:end,9);
dbpmse = immse(DBPFIT,y);

y=output4;
MAPFIT=MAPPARAM(1) + MAPPARAM(2)*X(k:end,2) + MAPPARAM(3)*X(k:end,3) + MAPPARAM(4)*X(k:end,4) + MAPPARAM(5)*X(k:end,5) + MAPPARAM(6)*X(k:end,6) + MAPPARAM(7)*X(k:end,7) + MAPPARAM(8)*X(k:end,8) + MAPPARAM(9)*X(k:end,9);
mapmse = immse(MAPFIT,y);

errors=[sqrt(ppmse) sqrt(sbpmse) sqrt(dbpmse) sqrt(mapmse)]
rmse_group=[rmse_group errors'];

corrs=[corr(PPFIT,output1,'type','Spearman') corr(SBPFIT,output2,'type','Spearman') corr(DBPFIT,output3,'type','Spearman') corr(MAPFIT,output4,'type','Spearman') ]
scc_group=[scc_group corrs'];

%real vs fit plot
f=figure;
ax(1)=subplot (4,1,1);
plot(t1(locsplotabp{mm,1}(k:end))./60,PPFIT,'r-',t1(locsplotabp{mm,1})./60,PP{l,1},'g-')
ylabel('mmHg')
xlabel('Time(min)')
legend('PP fit','PP real')

ax(2)=subplot (4,1,2);
plot(t1(locsplotabp{mm,1}(k:end))./60,SBPFIT,'r-',t1(locsplotabp{mm,1})./60,SBP{l,1},'g-')
ylabel('mmHg')
xlabel('Time(min)')
legend('SBP fit','SBP real')

ax(3)=subplot (4,1,3);
plot(t1(locsplotabp{mm,1}(k:end))./60,DBPFIT,'r-',t1(locsplotabp{mm,1})./60,DBP{l,1},'g-')
ylabel('mmHg')
xlabel('Time(min)')
legend('DBP fit','DBP real')

ax(4)=subplot (4,1,4);
plot(t1(locsplotabp{mm,1}(k:end))./60,MAPFIT,'r-',t1(locsplotabp{mm,1})./60,MAP{l,1},'g-')
ylabel('mmHg')
xlabel('Time(min)')
legend('MAP fit','MAP real')

set(f,'name','Causal system','numbertitle','off')
