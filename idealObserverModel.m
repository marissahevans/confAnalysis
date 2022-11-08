function [modelFitCircleNoise,sigma_m,sigma_p,sigma_s,ePt1,ePt2,ePtFB1,ePtFB2] = idealObserverModel(sigma_mAve,sigma_pAve,sigma_sAve, trialNum,lookUpTab1,endPoints, eptFB)

%This function creates sample data for model testing for the Meta Reaching
%Confidence experiment by Marissa Evans - 09/01/21
%% PARAMETERS

%Tablet specs
xCen = 238;
yCen = 134;
tabSize = [268 476];
matNan = nan(tabSize);

%Target locations
t = [xCen yCen];                             % target location for all trials

if nargin > 5
    %IF USING PARTICIPANT VALUES TO SIMULATE
    sigma_m = sigma_mAve;
    sigma_p = sigma_pAve;
    sigma_s = sigma_sAve;
    ePt = endPoints;
    ePtFB = eptFB;
    
else
    varM = 5^2;
    varP = 20^2;
    varS = 3^2;
    
    %IF SIMULATING RAW DATA FOR MODELS
    sigma_m = lognrnd(log((sigma_mAve^2)/sqrt(varM+sigma_mAve^2)),sqrt(log(varM/(sigma_mAve^2)+1)));            %motor noise
    sigma_p = lognrnd(log((sigma_pAve^2)/sqrt(varP+sigma_pAve^2)),sqrt(log(varP/(sigma_pAve^2)+1)));            %proprioceptive noise
    sigma_s = lognrnd(log((sigma_sAve^2)/sqrt(varS+sigma_sAve^2)),sqrt(log(varS/(sigma_sAve^2)+1)));              %setting noise
    ePt = t + randn(trialNum,2)*sigma_m;       %randomized for all trials
    ePtFB = t + randn(trialNum*2,2)*sigma_m;       %trials where confidence not reported
end
%Reach end points
ePt1 = ePt(:,1);
ePt2 = ePt(:,2);

%Reach end points NOT confidence reports
ePtFB1 = ePtFB(:,1);
ePtFB2 = ePtFB(:,2);

%Sensed locations based on enpoints and sigma_p
sensed = ePt + randn(trialNum,2)*sigma_p;

%% Simualtion
[X, Y] = meshgrid(1:100,1:100);
for ii = 1:length(lookUpTab1)
circSizes(ii) = interp2(X,Y,squeeze(lookUpTab1(ii,:,:)),sigma_p,sigma_m); %max expected gain circle size for test variables
end

maxDist = 1:length(circSizes);

%distances of sensed enpoints from the target
Rm = 1/sigma_m^2;
Rp = 1/sigma_p^2;
mu_pos = (Rm/(Rp+Rm)).*t + (Rp/(Rp+Rm)).*sensed;

%Distance from posterior to target
eucDistPos = sqrt((mu_pos(:,1) - t(1)).^2 + (mu_pos(:,2)-t(2)).^2);
eucDistPos(eucDistPos < 1) = 1;

%selected circle size from the look up table
circSelect = interp1(maxDist,circSizes,eucDistPos);

%Adding setting noise to all trials 
modelFitCircleNoise = circSelect + sigma_s*randn(trialNum,1);
modelFitCircleNoise(modelFitCircleNoise < 1) = 1;

end

