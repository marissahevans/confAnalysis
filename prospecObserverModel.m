function [modelFitCircleNoise,sigma_m,sigma_p,sigma_s,ePt1,ePt2,ePtFB1,ePtFB2] = prospecObserverModel(sigma_mAve,sigma_pAve,sigma_sAve, trialNum,radiusWac,endPoints,eptFB)
%PROSPECITIVE OBSERVER - PROSPREOCEPTIVE NOISE NOT USED
%This function creates sample data for model testing for the Meta Reaching
%Confidence experiment by Marissa Evans - 09/01/21
%% PARAMETERS

%Tablet specs
xCen = 238;
yCen = 134;
r = [10*ones(1,5),linspace(10,0,length(6:150))];

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
    
    %IF SIMULATING RAW DATA
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

%% Simulation

%Find best fitting circle based on knowledge of sigma_m only 

phit = raylcdf(0:length(radiusWac)-1,sigma_m);

gain = phit.*r; %probability of a hit vs points possible at each distance

%Radius point for each distance with the highest gain
maxGainCirc = radiusWac(max(find(gain == max(gain))));

modelFitCircleNoise = maxGainCirc + sigma_s*randn(trialNum,1);
modelFitCircleNoise(modelFitCircleNoise < 1) = 1;
end


