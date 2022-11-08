% Simulation of the control experiment
%clear ; close all
X = 0; Y = 0;       % target coordinates
target = [X Y];     % target location for all trials
%target = contResultsMat.targetLoc/100*4.65;

sigma_m = 15;        % true motor noise
sigma_p = 25;       % true proprioceptive noise
p_vec = 1:100;      % vector of possible proprioceptive noise values
m_vec = 1:100;      % vector of possible motor noise values

% Reliabilites - unknown to the experimenter
Rm = 1/sigma_m^2;
Rp = 1/sigma_p^2; 

numTrials = 300;    % number of simulated trials 

for ss = 1:100
% Reach end points - known to the experimenter
reaches = target + randn(numTrials,2)*sigma_m;
%reaches = contResultsMat.wacScreenEnd/100*4.65;

% Perceived locations - unknown to the experimenter
sensed = reaches + randn(numTrials,2)*sigma_p;

% Participant-reported endpoint - known to the experimenter
indicated = (Rm/(Rp+Rm))*target + (Rp/(Rp+Rm))*sensed;
%indicated = contResultsMat.mouseEndPt/100*4.65;
%fprintf(1,'True motor SD: %.1f\n',sigma_m);
%fprintf(1,'True proprioceptive SD: %.1f\n',sigma_p);

% figure; hold on
% scatter(reaches(:,1),reaches(:,2))
% scatter(indicated(:,1),indicated(:,2))
% scatter(target(1),target(2),'filled')
% legend('reaches','indicated','target')
% title('Reach end points vs Indicated End Points')
% axis([0 47.6 0 26.8])
% xlabel('cm')
% ylabel('cm')
% %axis equal
% set(gca, 'YDir','reverse','fontsize',18)
% text(5,5,'SigmaM = 16.7mm')
% text(5,7,'SigmaP = 29.8mm')

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
sigMmarg(ss) = m_vec(find(mmarg == max(mmarg)));

% fprintf(1,'Max marginal sigma_m =%.1f\n', ...
%     m_vec(find(mmarg == max(mmarg))));

pmarg = sum(NormPost,1);
pmarg = pmarg/sum(pmarg);
sigPmarg(ss) = p_vec(find(pmarg == max(pmarg)));

% fprintf(1,'Max marginal sigma_p =%.1f\n', ...
%     p_vec(find(pmarg == max(pmarg))));
end

figure
set(gcf,'color','w');
subplot(1,2,1)
histogram(sigMmarg,[0:2:40])
xline(sigma_m,'k--','LineWidth',2);
xlim([0 40])
ylim([0 110])
xlabel('sigma m')
ylabel('counts')
xticks([0:5:40])
legend('Simualtions','True sigma m','location','bestoutside')
box off
set(gca, 'TickDir', 'out', 'FontSize', 18)

subplot(1,2,2)
histogram(sigPmarg,[0:2:40])
xline(sigma_p,'k--','LineWidth',2);
xlim([0 40])
ylim([0 110])
xticks([0:5:40])
xlabel('sigma p')
ylabel('counts')
legend('Simualtions','True sigma p','location','bestoutside')
box off
set(gca, 'TickDir', 'out', 'FontSize', 18)



% % Plot the log likelihood surface for sigma_m + sigma_p
% figure
% surf(m_vec,p_vec,LL)
% title('Log Likelihood of \sigma_p and \sigma_m')
% xlabel('\sigma_m')
% ylabel('\sigma_p')
% zlabel('Log Likelihood')
% 
% % Log likelihood for sigma_m only
% figure
% plot(m_vec,LLe,'LineWidth',2)
% title('Log Likelihood of \sigma_m')
% xlabel('\sigma_m')
% ylabel('Log Likelihood')
% 
% % log likelihood for sigma_p only 
% figure
% surf(m_vec,p_vec,LLigivene,'LineWidth',2)
% title('Log Likelihood of \sigma_p')
% xlabel('\sigma_m')
% ylabel('\sigma_p')
% zlabel('Log Likelihood')
% 
% % marginal m posterior
% figure
% plot(m_vec,log(mmarg),'LineWidth',2)
% title('Marginal Log Posterior of \sigma_m')
% xlabel('\sigma_m')
% ylabel('Marginal Log Posterior')
% 
% % marginal p posterior
% figure
% plot(p_vec',log(pmarg)','LineWidth',2)
% title('Marginal Log Posterior of \sigma_p')
% xlabel('\sigma_p')
% ylabel('Marginal Log Posterior')
% 
% % direct computation of log(p(data|mu,sigma))
% % data is Nx2
% % mu is 1x2 or Nx2
% % sigma is a scalar SD of the isotropic 2-d Gaussian

function ll = log2disotnormal(data,mu,sigma)

centered = data - mu;
N = size(data,1);
ll = -N*log(2*pi) - 2*N*log(sigma) - (sum(centered(:).^2)/(2*sigma^2));
end
