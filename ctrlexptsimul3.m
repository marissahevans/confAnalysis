% Simulation of the control experiment
%clear ; close all
load('MeCont_metaReach_exp_S1_2021-08-24_15-45_controlresults.mat')
%load('ME_metaReach_exp_S1_2021-07-06_15-34_controlresults.mat')
%load('HP_metaReach_exp_S1_2021-07-07_14-28_controlresults.mat')
invtform = invert(displayInfo.tform);

%X = 0; Y = 0;       % target coordinates
%target = [X Y];     % target location for all trials
%target = contResultsMat.targetLoc/10*4.65;
[target(1) target(2)] = (transformPointsForward(invtform,contResultsMat.targetLoc(1),contResultsMat.targetLoc(2)));
[target(1) target(2)] = (transformPointsForward(tformMM,target(1),target(2)));

%sigma_m = 9;        % true motor noise
%sigma_p = 20;       % true proprioceptive noise
p_vec = 1:200;      % vector of possible proprioceptive noise values
m_vec = 1:200;      % vector of possible motor noise values

% Reliabilites - unknown to the experimenter
%Rm = 1/sigma_m^2;
%Rp = 1/sigma_p^2; 

%numTrials = 200;    % number of simulated trials 

% Reach end points - known to the experimenter
%reaches = target + randn(numTrials,2)*sigma_m;
%reaches = contResultsMat.wacScreenEnd/100*4.65;
reaches = [contResultsMat.wacEndPoint(61:end,1),contResultsMat.wacEndPoint(61:end,2)];
[reaches(:,1) reaches(:,2)] = transformPointsForward(tformMM,reaches(:,1),reaches(:,2));

% Perceived locations - unknown to the experimenter
%sensed = reaches + randn(numTrials,2)*sigma_p;

% Participant-reported endpoint - known to the experimenter
%indicated = (Rm/(Rp+Rm))*target + (Rp/(Rp+Rm))*sensed;
%indicated = contResultsMat.mouseEndPt/100*4.65;
mouse = [contResultsMat.mouseEndPt(61:end,1),contResultsMat.mouseEndPt(61:end,2)];
[indicated(:,1) indicated(:,2)] = (transformPointsForward(invtform,mouse(:,1),mouse(:,2)));
[indicated(:,1) indicated(:,2)] = (transformPointsForward(tformMM,indicated(:,1),indicated(:,2)));

%fprintf(1,'True motor SD: %.1f\n',sigma_m);
%fprintf(1,'True proprioceptive SD: %.1f\n',sigma_p);

% Simultaneously estimate sigma_m and sigma_p by ML

LL=zeros(length(m_vec),length(p_vec));
LLigivene=zeros(length(m_vec),length(p_vec));
LLe=zeros(length(m_vec),1);
for ii = 1:length(m_vec)        %loop over all sigma_m options
    RmTemp = 1/m_vec(ii)^2;
% log likelihood of sigma_m given target and endpoint:
    LLe(ii) = log2disotnormal(reaches,target,m_vec(ii));   
    for jj = 1:length(p_vec)    %loop over all sigma_pﬂ options
        RpTemp = 1/p_vec(jj)^2;
        meanigivene = (RpTemp/(RpTemp+RmTemp))*reaches + ...
            (RmTemp/(RpTemp+RmTemp))*target;
        SDigivene = (RpTemp/(RpTemp+RmTemp))*p_vec(jj);
% log likelihood of sigma_p given sensed location and endpoint:
        LLigivene(ii,jj) = log2disotnormal(indicated,meanigivene,SDigivene);
% log likelihood of the sigma_p/sigma_m pair:
        LL(ii,jj) = LLe(ii) + LLigivene(ii,jj);
    end
end

[i,j] = ind2sub(size(LL),find(LL == max(LL(:))));
%fprintf(1,'Max LL sigma_m = %.1f, sigma_p = %.1f\n',m_vec(i),p_vec(j));

% treat LL like a log posterior (i.e., treat prior as flat over the grid)
% and calculate marginals. First add a constant to all LL values so that
% the maximum is one (to minimize underflows) and normalize afterward.

NormPost = exp(LL - max(LL(:)));
mmarg = sum(NormPost,2);
mmarg = mmarg/sum(mmarg);
sigMmarg = m_vec(find(mmarg == max(mmarg)));

% fprintf(1,'Max marginal sigma_m =%.1f\n', ...
%     m_vec(find(mmarg == max(mmarg))));

pmarg = sum(NormPost,1);
pmarg = pmarg/sum(pmarg);
sigPmarg = p_vec(find(pmarg == max(pmarg)));

% fprintf(1,'Max marginal sigma_p =%.1f\n', ...
%     p_vec(find(pmarg == max(pmarg))));

figure; hold on
scatter(reaches(:,1),reaches(:,2))
scatter(indicated(:,1),indicated(:,2))
scatter(target(1),target(2),'filled')
legend('reaches','indicated','target')
title('Reach end points vs Indicated End Points')
axis([0 476 0 268])
xlabel('mm')
ylabel('mm')
%axis equal
set(gca,'fontsize',18)
text(50,50,['SigmaM = ' num2str(sigMmarg) 'mm'])
text(50,70,['SigmaP = ' num2str(sigPmarg) 'mm'])

%TESTING SECTION%
reachesOver = reaches(indicated(:,2) <= target(:,2),:);
reachesUnder = reaches(indicated(:,2) >= target(:,2),:);

indiOver = indicated(indicated(:,2) <= target(:,2),:);
indiUnder = indicated(indicated(:,2) >= target(:,2),:);

reachesLeft = reaches(indicated(:,1) <= target(:,1),:);
reachesRight = reaches(indicated(:,1) >= target(:,1),:);

indiLeft = indicated(indicated(:,1) <= target(:,1),:);
indiRight = indicated(indicated(:,1) >= target(:,1),:);

figure
subplot(2,2,1); hold on;
scatter(reachesOver(:,1),reachesOver(:,2))
scatter(indiOver(:,1),indiOver(:,2))
scatter(target(1),target(2),'filled')
scatter(median(reaches(:,1)),median(reaches(:,2)),'filled')
legend('reaches','indicated','target','median')
title('Indications are overshooting')
axis([0 476 0 268])
xlabel('mm')
ylabel('mm')
set(gca,'fontsize',18)
text(50,50,['SigmaM = ' num2str(sigMmarg) 'mm'])
text(50,70,['SigmaP = ' num2str(sigPmarg) 'mm'])

subplot(2,2,2); hold on;
scatter(reachesUnder(:,1),reachesUnder(:,2))
scatter(indiUnder(:,1),indiUnder(:,2))
scatter(target(1),target(2),'filled')
scatter(median(reaches(:,1)),median(reaches(:,2)),'filled')
legend('reaches','indicated','target','median')
title('Indications are undershooting')
axis([0 476 0 268])
xlabel('mm')
ylabel('mm')
set(gca,'fontsize',18)
text(50,50,['SigmaM = ' num2str(sigMmarg) 'mm'])
text(50,70,['SigmaP = ' num2str(sigPmarg) 'mm'])

subplot(2,2,3); hold on;
scatter(reachesLeft(:,1),reachesLeft(:,2))
scatter(indiLeft(:,1),indiLeft(:,2))
scatter(target(1),target(2),'filled')
scatter(median(reaches(:,1)),median(reaches(:,2)),'filled')
legend('reaches','indicated','target','median')
title('Indications to the left')
axis([0 476 0 268])
xlabel('mm')
ylabel('mm')
set(gca,'fontsize',18)
text(50,50,['SigmaM = ' num2str(sigMmarg) 'mm'])
text(50,70,['SigmaP = ' num2str(sigPmarg) 'mm'])

subplot(2,2,4); hold on;
scatter(reachesRight(:,1),reachesRight(:,2))
scatter(indiRight(:,1),indiRight(:,2))
scatter(target(1),target(2),'filled')
scatter(median(reaches(:,1)),median(reaches(:,2)),'filled')
legend('reaches','indicated','target','median')
title('Indications to the right')
axis([0 476 0 268])
xlabel('mm')
ylabel('mm')
set(gca,'fontsize',18)
text(50,50,['SigmaM = ' num2str(sigMmarg) 'mm'])
text(50,70,['SigmaP = ' num2str(sigPmarg) 'mm'])



figure
subplot(1,2,1)
histogram(sigMmarg)
subplot(1,2,2)
histogram(sigPmarg)

% Plot the log likelihood surface for sigma_m + sigma_p
figure
surf(m_vec,p_vec,LL)
title('Log Likelihood of \sigma_p and \sigma_m')
xlabel('\sigma_m')
ylabel('\sigma_p')
zlabel('Log Likelihood')

% Log likelihood for sigma_m only
figure
plot(m_vec,LLe,'LineWidth',2)
title('Log Likelihood of \sigma_m')
xlabel('\sigma_m')
ylabel('Log Likelihood')

% log likelihood for sigma_p only 
figure
surf(m_vec,p_vec,LLigivene,'LineWidth',2)
title('Log Likelihood of \sigma_p')
xlabel('\sigma_m')
ylabel('\sigma_p')
zlabel('Log Likelihood')

% marginal m posterior
figure
plot(m_vec,log(mmarg),'LineWidth',2)
title('Marginal Log Posterior of \sigma_m')
xlabel('\sigma_m')
ylabel('Marginal Log Posterior')

% marginal p posterior
figure
plot(p_vec',log(pmarg)','LineWidth',2)
title('Marginal Log Posterior of \sigma_p')
xlabel('\sigma_p')
ylabel('Marginal Log Posterior')

% direct computation of log(p(data|mu,sigma))
% data is Nx2
% mu is 1x2 or Nx2
% sigma is a scalar SD of the isotropic 2-d Gaussian

function ll = log2disotnormal(data,mu,sigma)

centered = data - mu;
N = size(data,1);
ll = -N*log(2*pi) - 2*N*log(sigma) - (sum(centered(:).^2)/(2*sigma^2));
end
