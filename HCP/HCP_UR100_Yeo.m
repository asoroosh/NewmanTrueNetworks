clear


N = 114;
nsub = 99;

T = 1200; 

bias = @(E1,E0) ((E1-E0)./E0.*100);

%Sim me some nets
%A = binornd(1,0.5,[N N M]);
%for i = 1:M
%    A(:,:,i)=triu(A(:,:,i),1)+triu(A(:,:,i),1)';
%end
addpath(genpath('/Users/sorooshafyouni/Home/GitClone/DVARS'))
addpath(genpath('/Users/sorooshafyouni/Home/GitClone/xDF'))
addpath(genpath('/Users/sorooshafyouni/Home/GitClone/NewmanTrueNetworks'))
addpath(genpath('/Users/sorooshafyouni/Home/matlab/Ext/BCT/2017_01_15_BCT'))

load('/Users/sorooshafyouni/Home/BCF/BCFAnal/FC/S/HCP_100Unrel_SubList.mat')
SubList = HCP_10Unrel_SubList;
SubList([3]) = [];

for i = 1:nsub
    disp([num2str(i) ' - On subject: ' SubList{i}])
    load(['/Users/sorooshafyouni/Home/BCF/BCFAnal/FC/100HCPTimeSeries/Yeo/HCP_FPP_' SubList{i} '_OnlyMTS.mat'])
    %load(['/Users/sorooshafyouni/Home/BCF/BCFAnal/FC/ICA200/FPP/HCP_FPP_' SubList{i} '_OnlyMTS.mat'])
    
%     mts = GSRme(mts,T);
    
%     [~,Stat] = xDF(mts',T);
%     mat_tmp = Stat.z.rzf;
%     mat_tmp(mat_tmp<0) = 0; 
%     mat_tmp = z2p_bon(mat_tmp);
    
    %mat_tmp = corr(mts);
    %mat_tmp(mat_tmp<0) = 0; 
    %mat_tmp = threshold_proportional(mat_tmp,0.15);
    %mat_tmp(mat_tmp>0) = 1;
    
    mat_tmp = corr(mts);
    mat_tmp(mat_tmp<0) = 0; 
    mat_tmp = atanh(mat_tmp).*sqrt(T-3);
    mat_tmp = z2p_fdr(mat_tmp);
    mat_tmp(mat_tmp>0) = 1;    
    
    mat_tmp(1:N+1:end) = 0;
    A(:,:,i) = mat_tmp;
    
    Naive_dgr(:,i) = degrees_und(mat_tmp);
    
    clear mat_tmp;
    
    [~,Stat] = xDF(mts',T,'truncate','adaptive','TVOff');
    mat_tmp = Stat.z;
    mat_tmp(mat_tmp<0) = 0;
    mat_tmp = z2p_fdr(mat_tmp);
    mat_tmp(mat_tmp>0) = 1; 
    
    mat_tmp(1:N+1:end) = 0;
    Axdf(:,:,i) = mat_tmp;
    
    xDF_dgr(:,i) = degrees_und(mat_tmp);
    
    
%     [~,mat_tmp] = CRBCF(mts',T);
%     mat_tmp(mat_tmp<0) = 0;
%     mat_tmp = z2p_fdr(mat_tmp);
%     mat_tmp(mat_tmp>0) = 1;   
%     mat_tmp(1:N+1:end) = 0;
%     Acr(:,:,i) = mat_tmp;
end

disp(['Alpha-beta on Corr z-scores'])
[Q_nodal,bet1_nodal,alp1_nodal,rho1_nodal,FDR_nodal,itr_nodal] = EM_nodal_bu(A);
%[Q,cr,bet1,alp1,rho1,itr] = EM_bu(A);

disp(['Alpha-beta on xDF z-scores'])
[Q_xdf_nodal,bet1_xdf_nodal,alp1_xdf_nodal,rho1_xdf_nodal,FDR_xdf_nodal,itr_xdf_nodal] = EM_nodal_bu(Axdf);

% disp(['Alpha-beta on xDF z-scores'])
% [Q_cr_nodal,bet1_cr_nodal,alp1_cr_nodal,rho1_cr_nodal,FDR_cr_nodal,itr_cr_nodal] = EM_nodal_bu(Acr);

%%%
xDFd = mean(xDF_dgr,2);
Nd = mean(Naive_dgr,2);
[r,p] = corr(alp1_xdf_nodal',xDFd)

figure; 
subplot(2,1,1)
hold on;
g = bar([bias(bet1_xdf_nodal,bet1_nodal)]);
title('FPR')
legend('xDF','BH')
ylim([-50 50])

subplot(2,1,2)
hold on;
g = bar([bias(alp1_xdf_nodal,alp1_nodal)]);
title('TPR')
xlabel('Nodes')
legend('xDF','BH')
ylim([-50 50])

%=====================================

figure; 
subplot(3,1,1)
hold on;
g = bar([bet1_nodal;bet1_xdf_nodal]');
g(1).FaceColor = 'r';
g(2).FaceColor = 'b';
title('FPR')
legend('Naive','xDF')
ylim([0 1])

subplot(3,1,2)
hold on;
g = bar([alp1_nodal;alp1_xdf_nodal]');
g(1).FaceColor = 'r';
g(2).FaceColor = 'b';
title('TPR')
xlabel('Nodes')
legend('Naive','xDF')
ylim([0 1])

subplot(3,1,3)
hold on;
g = bar([alp1_nodal-alp1_xdf_nodal; bet1_nodal-bet1_xdf_nodal]');
g(1).FaceColor = 'r';
g(2).FaceColor = 'b';
title('TPR')
xlabel('Nodes')
legend('Difference In TPR','Difference in FPR')




