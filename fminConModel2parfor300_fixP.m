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
lb = 1;                            %lower bound
ub = 50;                         %upper bound
init = rand*(ub-lb)+lb;                  %initiation
options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off'); %fmincon options

numTrials = 300;
numSims = 20;


sigS2A = zeros(1,numSims);
sigS2B = zeros(1,numSims);
sigS2C = zeros(1,numSims);

minNLL2A = zeros(1,numSims);
minNLL2B = zeros(1,numSims);
minNLL2C = zeros(1,numSims);

for ii = 1:numSims
    [sigP2A(ii),sigM2A(ii)] = ctrlFit(target1,[reach1X(ii,:); reach1Y(ii,:)]',[indic1X(ii,:); indic1Y(ii,:)]',sigMmax,sigPmax);
    [sigP2B(ii),sigM2B(ii)] = ctrlFit(target1,[reach2X(ii,:); reach2Y(ii,:)]',[indic2X(ii,:); indic2Y(ii,:)]',sigMmax,sigPmax);
    [sigP2C(ii),sigM2C(ii)] = ctrlFit(target1,[reach3X(ii,:); reach3Y(ii,:)]',[indic3X(ii,:); indic3Y(ii,:)]',sigMmax,sigPmax);
end

parfor ii = 1:numSims
    
    m2nLogLA = @(p) likelihoodFit2(sigM2A(ii),sigP2A(ii),p,numTrials,[ePts1X(ii,1:numTrials)',ePts1Y(ii,1:numTrials)'],modelFitCircleNoise1(ii,1:numTrials),t,distTestEndpts1(ii,1:numTrials),fit2LookUpMat,maxDistAll,sigPmax,distFromTarget);
    m2nLogLB = @(p) likelihoodFit2(sigM2B(ii),sigP2B(ii),p,numTrials,[ePts2X(ii,1:numTrials)',ePts2Y(ii,1:numTrials)'],modelFitCircleNoise2(ii,1:numTrials),t,distTestEndpts2(ii,1:numTrials),fit2LookUpMat,maxDistAll,sigPmax,distFromTarget);
    m2nLogLC = @(p) likelihoodFit2(sigM2C(ii),sigP2C(ii),p,numTrials,[ePts3X(ii,1:numTrials)',ePts3Y(ii,1:numTrials)'],modelFitCircleNoise3(ii,1:numTrials),t,distTestEndpts3(ii,1:numTrials),fit2LookUpMat,maxDistAll,sigPmax,distFromTarget);
    
    [sigS2A(ii), minNLL2A(ii)] = fmincon(m2nLogLA,init,[],[],[],[],lb,ub,[],options); %find minimums - model 2 w/ data 1
    [sigS2B(ii), minNLL2B(ii)] = fmincon(m2nLogLB,init,[],[],[],[],lb,ub,[],options); %find minimums - model 2 w/ data 2
    [sigS2C(ii), minNLL2C(ii)] = fmincon(m2nLogLC,init,[],[],[],[],lb,ub,[],options); %find minimums - model 2 w/ data 3
  
end

estP2A = [sigM2A', sigP2A', sigS2A'];
estP2B = [sigM2B', sigP2B', sigS2B'];
estP2C = [sigM2C', sigP2C', sigS2C'];

save model2output_fixP.mat estP2A estP2B estP2C minNLL2A minNLL2B minNLL2C
