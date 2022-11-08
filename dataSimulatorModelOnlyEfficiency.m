%% SIMULATING DATA (Confidence Reports) FROM MODELS

% This section of the code takes ground truth mean values for sigma_m, sigma_p
% and sigma_s and runs three functions which select a circle size based on
% the observers perspective. The output is a (simulations x trials) matrix
% of the circle sizes selected, a vector each of sigma_m, sigma_p and
% sigma_s selections (currently distributed at 15, 15 and 1 respectively),
% and the endpoint coordinates for each trial based on the target (same for
% all trials) and sigma_m. 

% idealObserverModel.m, retroObserverModel.m, and prospecObserverModel.m
% can be updated depending if you want to simulate data for a model fit
% (based on an underlying ground truth) or from parameters based on the
% participants data. If this script does NOT run due to an error in ePT
% then these scripts are configured to take parameters and need to be
% switched by commenting out the appropriate section at the top. 

load model1LookUpTab %look up table for model 1
load model2LookUpTab %look up table for model 2
load model3LookUpTab %look up table for model 3

radiusWac = linspace(0,69,150);
sigma_m = 35;
sigma_p = 30;
sigma_s = 6;
numTrials = 1000;
numSims = 1;
target1 = [0,0];


for ii = 1:numSims
    [modelFitCircleNoise1(ii,:),sigma_m1(ii),sigma_p1(ii),sigma_s1(ii),ePts1X(ii,:),ePts1Y(ii,:)] = idealObserverModel(sigma_m,sigma_p,sigma_s,numTrials,fit1LookUpMat);
    
    %Control experiment
    Rm = 1/sigma_m1(ii)^2;
    Rp = 1/sigma_p1(ii)^2;
  
    reaches = target1 + randn(numTrials,2)*sigma_m1(ii);
    sensed = reaches + randn(numTrials,2)*sigma_p1(ii);
    indicated = (Rm/(Rp+Rm))*target1 + (Rp/(Rp+Rm))*sensed;
    
    reach1X(ii,:) = reaches(:,1);
    reach1Y(ii,:) = reaches(:,2);
    indic1X(ii,:) = indicated(:,1);
    indic1Y(ii,:) = indicated(:,2);
    
%     [modelFitCircleNoise2(ii,:),sigma_m2(ii),sigma_p2(ii),sigma_s2(ii),ePts2X(ii,:),ePts2Y(ii,:)] = retroObserverModel(sigma_m,sigma_p,sigma_s,numTrials,fit2LookUpMat);
%     
%     %Control experiment
%     Rm = 1/sigma_m2(ii)^2;
%     Rp = 1/sigma_p2(ii)^2;
%   
%     reaches = target1 + randn(numTrials,2)*sigma_m2(ii);
%     sensed = reaches + randn(numTrials,2)*sigma_p2(ii);
%     indicated = (Rm/(Rp+Rm))*target1 + (Rp/(Rp+Rm))*sensed;
%     
%     reach2X(ii,:) = reaches(:,1);
%     reach2Y(ii,:) = reaches(:,2);
%     indic2X(ii,:) = indicated(:,1);
%     indic2Y(ii,:) = indicated(:,2);
    
    [modelFitCircleNoise3(ii,:),sigma_m3(ii),sigma_p3(ii),sigma_s3(ii),ePts3X(ii,:),ePts3Y(ii,:)] = prospecObserverModel(sigma_m,sigma_p,sigma_s,numTrials,radiusWac);

    %Control experiment
    Rm = 1/sigma_m3(ii)^2;
    Rp = 1/sigma_p3(ii)^2;
  
    reaches = target1 + randn(numTrials,2)*sigma_m3(ii);
    sensed = reaches + randn(numTrials,2)*sigma_p3(ii);
    indicated = (Rm/(Rp+Rm))*target1 + (Rp/(Rp+Rm))*sensed;
    
    reach3X(ii,:) = reaches(:,1);
    reach3Y(ii,:) = reaches(:,2);
    indic3X(ii,:) = indicated(:,1);
    indic3Y(ii,:) = indicated(:,2);
end

if sum(isnan(modelFitCircleNoise1)) > 0
    modelFitCircleNoise1(isnan(modelFitCircleNoise1)) = nanmean(modelFitCircleNoise1);
%elseif sum(isnan(modelFitCircleNoise2)) > 0
    %modelFitCircleNoise2(isnan(modelFitCircleNoise2)) = nanmean(modelFitCircleNoise2);
elseif sum(isnan(modelFitCircleNoise3)) > 0
    modelFitCircleNoise3(isnan(modelFitCircleNoise3)) = nanmean(modelFitCircleNoise3);
end

    
%Points earned by simulated data
dist1 = sqrt((ePts1X - t(1)).^2 + (ePts1Y-t(2)).^2);
%dist2 = sqrt((ePts2X - t(1)).^2 + (ePts2Y-t(2)).^2);
dist3 = sqrt((ePts3X - t(1)).^2 + (ePts3Y-t(2)).^2);
pointsPos1 = 10 - (10/ 67.67 * (modelFitCircleNoise1-2.33));
%pointsPos2 = 10 - (10/ 67.67 * (modelFitCircleNoise2-2.33));
pointsPos3 = 10 - (10/ 67.67 * (modelFitCircleNoise3-2.33));

pointsEarn1 = nan(size(modelFitCircleNoise1));
%pointsEarn2 = nan(size(modelFitCircleNoise1));
pointsEarn3 = nan(size(modelFitCircleNoise1));

for ii = 1:numSims
    for jj = 1:numTrials
        
        %model1
        if modelFitCircleNoise1(ii,jj) > 150
            pointsEarn1(ii,jj) = 0;
        elseif modelFitCircleNoise1(ii,jj) <= 2.33 && dist1(ii,jj)-1.87 <= 2.33
            pointsEarn1(ii,jj) = 10;
        elseif modelFitCircleNoise1(ii,jj) >= dist1(ii,jj)-1.87
            pointsEarn1(ii,jj) = pointsPos1(ii,jj);
        else
            pointsEarn1(ii,jj) = 0;
        end
        
%         %model2
%         if modelFitCircleNoise2(ii,jj) > 150
%             pointsEarn2(ii,jj) = 0;
%         elseif modelFitCircleNoise2(ii,jj) <= 2.33 && dist2(ii,jj)-1.87 <= 2.33
%             pointsEarn2(ii,jj) = 10;
%         elseif modelFitCircleNoise2(ii,jj) >= dist2(ii,jj)-1.87
%             pointsEarn2(ii,jj) = pointsPos2(ii,jj);
%         else
%             pointsEarn2(ii,jj) = 0;
%         end
        
        %model3
        if modelFitCircleNoise3(ii,jj) > 150
            pointsEarn3(ii,jj) = 0;
        elseif modelFitCircleNoise3(ii,jj) <= 2.33 && dist3(ii,jj)-1.87 <= 2.33
            pointsEarn3(ii,jj) = 10;
        elseif modelFitCircleNoise3(ii,jj) >= dist3(ii,jj)-1.87
            pointsEarn3(ii,jj) = pointsPos3(ii,jj);
        else
            pointsEarn3(ii,jj) = 0;
        end
    end
end

mean1 = mean(pointsEarn1)
mean3 = mean(pointsEarn3)

effRatio = mean3/mean1

probCap1 = sum(pointsEarn1~=0)/numTrials
probCap2 = sum(pointsEarn3~=0)/numTrials

%save simOutputEfficiency.mat ePts1X ePts1Y ePts2X ePts2Y ePts3X ePts3Y sigma_m1 sigma_m2 sigma_m3 sigma_p1 sigma_p2 sigma_p3 sigma_s1 sigma_s2 sigma_s3 modelFitCircleNoise1 modelFitCircleNoise2 modelFitCircleNoise3...
    %target1 reach1X reach1Y reach2X reach2Y reach3X reach3Y indic1X indic1Y indic2X indic2Y indic3X indic3Y pointsEarn1 pointsEarn2 pointsEarn3
    
save simOutputEfficiency.mat ePts1X ePts1Y ePts3X ePts3Y sigma_m1 sigma_m3 sigma_p1 sigma_p3 sigma_s1 sigma_s3 modelFitCircleNoise1 modelFitCircleNoise3...
    target1 reach1X reach1Y reach3X reach3Y indic1X indic1Y indic3X indic3Y pointsEarn1 pointsEarn3

