%Simulating control experiment
clear ; close all
X = 0; Y = 0;       %target coordinates
target = [X Y];     %target location for all trials

sigma_m = 9;        %motor noise
sigma_p = 20;       %proprioceptive noise
p_vec = 1:200;       %vector of possible proprioceptive noise values
m_vec = 1:200;       %vector of possible motor noise values

numTrials = 200;    %number of times the fit runs / number simulated trials 

%Reach end points - known to experimenter
reaches = target + randn(numTrials,2)*sigma_m;

%perceived locations - unknown to experimenter
sensed = reaches + randn(numTrials,2)*sigma_p;

%reliabilites - unknown to experimenter
Rm = 1/sigma_m^2;
Rp = 1/sigma_p^2;

%participant reported endpoint - known to experimenter
indicated = (Rm/(Rp+Rm))*target + (Rp/(Rp+Rm))*sensed;
reacherrs = reaches - target;
sigmamest = std(reacherrs(:),1);  % ML SD of motor
fprintf(1,'True motor SD: %.1f\n',sigma_m);
fprintf(1,'True proprioceptive SD: %.1f\n',sigma_p);
fprintf(1,'Estimated motor SD from noisy reaches alone: %.1f\n',sigmamest);
Rmest = 1/sigmamest^2;

% Version 1: Assume sigmamest and use the variability of the
% indicated location to solve for sigma_p
%
% In particular: sigma_i^2 = (R_p/(R_p+R_m))^2(sigma_p^2+\sigma_m^2)
% Substituting for sigma_p leads to a quadratic with the
% two solutions below

indicsd = std(indicated(:),1); % SD of indicated positions
a = indicsd^2/sigmamest^4;
b =(2*indicsd^2/sigmamest^2) - 1;
c = indicsd^2 - sigmamest^2;
root = sqrt(b^2-4*a*c);
soln1 = (-a-root)/(2*a);
soln2 = (-a+root)/(2*a);
fprintf(1, ...
    'Proprioception variance estimates based on indicated-location variability: ');
prtnum(soln1);
fprintf(1,', ');
prtnum(soln2);
fprintf(1,'\nCorresponding SD values: ');
prtroot(soln1);
fprintf(1,', ');
prtroot(soln2);
fprintf(1,'\n');

% Version 2: Estimate proprioceptive SD by ML, fixing the motor SD estimate

LLsigp2 = zeros(length(p_vec),1);
for jj = 1:length(p_vec)    %loop over all sigma_p options
    RpTemp = 1/p_vec(jj)^2;
    senseest = ((RpTemp+Rmest)*indicated - Rmest*target)/RpTemp; %calculating the sensed location
    LLsigp2(jj) = log2disotnormal(senseest,reaches,p_vec(jj)); %log likelihood of sigma_p given sensed location and endpoint
end

figure
plot(p_vec,LLsigp2,'LineWidth',2)
title('Log Likelihood of \sigma_p holding \sigma_m fixed')
xlabel('\sigma_p')
ylabel('Log Likelihood')
pest2 = p_vec(find(LLsigp2 == max(LLsigp2)));
fprintf(1,'sigma_p estimate using fixed sigma_m estimate: %.1f\n',pest2);

% Version 3: Simultaneously estimate sigma_m and sigma_p by ML

LLsigFit=zeros(length(m_vec),length(p_vec));
LLsigp=zeros(length(m_vec),length(p_vec));
LLsigm=zeros(length(m_vec),1);
%Estimate sigma_p and sigma_m
for ii = 1:length(m_vec)        %loop over all sigma_m options
    RmTemp = 1/m_vec(ii)^2;
    LLsigm(ii) = log2disotnormal(reaches,target,m_vec(ii)); %log likelihood of sigma_m given target and endpoint   
    for jj = 1:length(p_vec)    %loop over all sigma_pß options
        RpTemp = 1/p_vec(jj)^2;
        senseest = ((RpTemp+RmTemp)*indicated - RmTemp*target)/RpTemp; %calculating the sensed location
        LLsigp(ii,jj) = log2disotnormal(senseest,reaches,p_vec(jj)); %log likelihood of sigma_p given sensed location and endpoint
        LLsigFit(ii,jj) = LLsigm(ii) + LLsigp(ii,jj); %Put all summed values for a sigma_p/sigma_m pair in a matrix
    end
end

[i,j] = ind2sub(size(LLsigFit),find(LLsigFit == max(LLsigFit(:))));
fprintf(1,'Max LL sigma_m = %.1f, sigma_p = %.1f\n',m_vec(i),p_vec(j));

% treat LL like a log posterior (i.e., treat prior as flat over the grid)
% and calculate marginals. First add a constant to all LL values so that
% the maximum is one (to minimize underflows) and normalize afterward.

sigFitNormPost = exp(LLsigFit - max(LLsigFit(:)));
mmarg = sum(sigFitNormPost,2);
mmarg = mmarg/sum(mmarg);
fprintf(1,'Max marginal sigma_m =%.1f\n', ...
    m_vec(find(mmarg == max(mmarg))));
pmarg = sum(sigFitNormPost,1);
pmarg = pmarg/sum(pmarg);
fprintf(1,'Max marginal sigma_p =%.1f\n', ...
    p_vec(find(pmarg == max(pmarg))));

% Plot the log likelihood surface for sigma_m + sigma_p
figure
surf(m_vec,p_vec,LLsigFit)
title('Log Likelihood of \sigma_p and \sigma_m')
xlabel('\sigma_m')
ylabel('\sigma_p')
zlabel('Log Likelihood')

% Log likelihood for sigma_m only
figure
plot(m_vec,LLsigm,'LineWidth',2)
title('Log Likelihood of \sigma_m')
xlabel('\sigma_m')
ylabel('Log Likelihood')

% log likelihood for sigma_p only 
figure
surf(m_vec,p_vec,LLsigp,'LineWidth',2)
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
% mu is 1x2
% sigma is a scalar SD of the isotropic 2-d Gaussian

function ll = log2disotnormal(data,mu,sigma)

centered = data - mu;
N = size(data,1);
ll = -N*log(2*pi) - 2*N*log(sigma) - (sum(centered(:).^2)/(2*sigma^2));
end

function prtnum(x)

fprintf(1,'%.1f',real(x));
if ~isreal(x)
    fprintf(1,' + %.1fi',imag(x));
end
end

function prtroot(x)

if isreal(x) && x >= 0
    fprintf(1,'%.1f',sqrt(x));
else
    fprintf(1,'N.A.');
end
end