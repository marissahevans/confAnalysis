%% fmincon for prospective model #2
parpool(20)

load simOutput
load model2LookUpTab

%Tablet specs
radiusWac = linspace(0,69,150);             %max circle size is 69mm (at 70mm score drops to 0)
                                            %using 150 steps because
                                            %originally it was 150 pixels
                                            %in length

xCen = 238;                                 %center of tablet x coordinate
yCen = 134;                                 %cetner of tablet y coordinate
t = [xCen yCen];                            %target location for all trials
tabSize = [268 476];                        %outside bounds of tablet space
matNan = nan(tabSize);                      %matrix of nans the size of the tablet
sigPmax = 1:100;                           %vector of test sigma_p values
sigMmax = 1:100;                           %vector of test sigma_m values

distTestEndpts1 = sqrt((ePts1X - t(1)).^2 + (ePts1Y-t(2)).^2); %distances from endpoints to target
distTestEndpts2 = sqrt((ePts2X - t(1)).^2 + (ePts2Y-t(2)).^2); %distances from endpoints to target
distTestEndpts3 = sqrt((ePts3X - t(1)).^2 + (ePts3Y-t(2)).^2); %distances from endpoints to target

endPtsFBdist1 = sqrt((ePtsFB1X - t(1)).^2 + (ePtsFB1Y-t(2)).^2);
endPtsFBdist2 = sqrt((ePtsFB2X - t(1)).^2 + (ePtsFB2Y-t(2)).^2);
endPtsFBdist3 = sqrt((ePtsFB3X - t(1)).^2 + (ePtsFB3Y-t(2)).^2);

modelFitCircleNoise1(modelFitCircleNoise1 <1) = 1;
modelFitCircleNoise2(modelFitCircleNoise2 <1) = 1;
modelFitCircleNoise3(modelFitCircleNoise3 <1) = 1;

% Find all distances to target on tablet 
for nn = 1:size(matNan,1)
    for mm = 1:size(matNan,2)
        distFromTarget(nn,mm) = sqrt((mm - t(1))^2 + (nn-t(2))^2);
    end
end

maxDistAll = 0:max(distFromTarget(:)); %range of distance to max distance on the tablet. 
distFromTarget(distFromTarget > max(maxDistAll)) = max(maxDistAll);

%fmincon parameters
lb = [1,1,1];                            %lower bound
ub = [100,100,40];                         %upper bound
init = rand*(ub-lb)+lb;                  %initiation
options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off'); %fmincon options

numTrials = 300;
numSims = 20;


estP2A = zeros(numSims,3);
estP2B = zeros(numSims,3);
estP2C = zeros(numSims,3);

minNLL2A = zeros(1,numSims);
minNLL2B = zeros(1,numSims);
minNLL2C = zeros(1,numSims);

parfor ii = 1:numSims
    
    m2nLogLA = @(p) likelihoodFit2(p(1),p(2),p(3),numTrials,[ePts1X(ii,1:numTrials)',ePts1Y(ii,1:numTrials)'],modelFitCircleNoise1(ii,1:numTrials),t,distTestEndpts1(ii,1:numTrials),fit2LookUpMat,maxDistAll,sigPmax,distFromTarget,endPtsFBdist1(ii,:),target1,[reach1X(ii,:); reach1Y(ii,:)]',[indic1X(ii,:); indic1Y(ii,:)]');
    m2nLogLB = @(p) likelihoodFit2(p(1),p(2),p(3),numTrials,[ePts2X(ii,1:numTrials)',ePts2Y(ii,1:numTrials)'],modelFitCircleNoise2(ii,1:numTrials),t,distTestEndpts2(ii,1:numTrials),fit2LookUpMat,maxDistAll,sigPmax,distFromTarget,endPtsFBdist2(ii,:),target1,[reach2X(ii,:); reach2Y(ii,:)]',[indic2X(ii,:); indic2Y(ii,:)]');
    m2nLogLC = @(p) likelihoodFit2(p(1),p(2),p(3),numTrials,[ePts3X(ii,1:numTrials)',ePts3Y(ii,1:numTrials)'],modelFitCircleNoise3(ii,1:numTrials),t,distTestEndpts3(ii,1:numTrials),fit2LookUpMat,maxDistAll,sigPmax,distFromTarget,endPtsFBdist3(ii,:),target1,[reach3X(ii,:); reach3Y(ii,:)]',[indic3X(ii,:); indic3Y(ii,:)]');
    
    [estP2A(ii,:), minNLL2A(ii)] = fmincon(m2nLogLA,init,[],[],[],[],lb,ub,[],options); %find minimums - model 2 w/ data 1
    [estP2B(ii,:), minNLL2B(ii)] = fmincon(m2nLogLB,init,[],[],[],[],lb,ub,[],options); %find minimums - model 2 w/ data 2
    [estP2C(ii,:), minNLL2C(ii)] = fmincon(m2nLogLC,init,[],[],[],[],lb,ub,[],options); %find minimums - model 2 w/ data 3
    
end

save model2output.mat estP2A estP2B estP2C minNLL2A minNLL2B minNLL2C
