%% fmincon for ideal model #1
parpool(20) %number of workers

load simOutput %imports the endpoints and circle sizes from the simualted data
load model1LookUpTab %imports the look up table previously calculated for this model

%Tablet specs
radiusWac = linspace(0,69,150);             %max circle size is 69mm (at 70mm score drops to 0)
                                            %using 150 steps because
                                            %originally it was 150 pixels
                                            %in length

xCen = 238;                                 %center of tablet x coordinate in mm
yCen = 134;                                 %cetner of tablet y coordinate in mm
t = [xCen yCen];                            %target location for all trials
tabSize = [268 476];                        %outside bounds of tablet space in mm
matNan = nan(tabSize);                      %matrix of nans the size of the tablet
sigPmax = 1:100;                             %vector of test sigma_p values
sigMmax = 1:100;                             %vector of test sigma_m values

distFromTarget = matNan;
distTestEndpts1 = sqrt((ePts1X - t(1)).^2 + (ePts1Y-t(2)).^2); %distances from endpoints to target
distTestEndpts2 = sqrt((ePts2X - t(1)).^2 + (ePts2Y-t(2)).^2); %distances from endpoints to target
distTestEndpts3 = sqrt((ePts3X - t(1)).^2 + (ePts3Y-t(2)).^2); %distances from endpoints to target

modelFitCircleNoise1(modelFitCircleNoise1 <1) = 1;
modelFitCircleNoise2(modelFitCircleNoise2 <1) = 1;
modelFitCircleNoise3(modelFitCircleNoise3 <1) = 1;

% Find all distances to target on tablet in mm
for nn = 1:size(matNan,1)
    for mm = 1:size(matNan,2)
        distFromTarget(nn,mm) = sqrt((mm - t(1))^2 + (nn-t(2))^2);
    end
end

maxDistAll = 0:max(distFromTarget(:)); %range of distance to max distance on the tablet. 
distFromTarget(distFromTarget > max(maxDistAll)) = max(maxDistAll);

%fmincon parameters
lb = 1;                          %lower bound
ub = 50;                         %upper bound
init = rand*(ub-lb)+lb;                  %initiation
options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off'); %fmincon options

numTrials = 300;
numSims = 20;

LNT = length(numTrials);

sigS1A = zeros(1,numSims);
sigS1B = zeros(1,numSims);
sigS1C = zeros(1,numSims);

minNLL1A = zeros(1,numSims);
minNLL1B = zeros(1,numSims);
minNLL1C = zeros(1,numSims);

for ii = 1:numSims
    [sigP1A(ii),sigM1A(ii)] = ctrlFit(target1,[reach1X(ii,:); reach1Y(ii,:)]',[indic1X(ii,:); indic1Y(ii,:)]',sigMmax,sigPmax);
    [sigP1B(ii),sigM1B(ii)] = ctrlFit(target1,[reach2X(ii,:); reach2Y(ii,:)]',[indic2X(ii,:); indic2Y(ii,:)]',sigMmax,sigPmax);
    [sigP1C(ii),sigM1C(ii)] = ctrlFit(target1,[reach3X(ii,:); reach3Y(ii,:)]',[indic3X(ii,:); indic3Y(ii,:)]',sigMmax,sigPmax);
end

parfor ii = 1:numSims
    
    m1nLogLA = @(p) likelihoodFit1(sigM1A(ii),sigP1A(ii),p,numTrials,[ePts1X(ii,1:numTrials)',ePts1Y(ii,1:numTrials)'],modelFitCircleNoise1(ii,1:numTrials),t,distTestEndpts1(ii,1:numTrials),fit1LookUpMat,maxDistAll,sigMmax,sigPmax,distFromTarget);
    m1nLogLB = @(p) likelihoodFit1(sigM1B(ii),sigP1B(ii),p,numTrials,[ePts2X(ii,1:numTrials)',ePts2Y(ii,1:numTrials)'],modelFitCircleNoise2(ii,1:numTrials),t,distTestEndpts2(ii,1:numTrials),fit1LookUpMat,maxDistAll,sigMmax,sigPmax,distFromTarget);
    m1nLogLC = @(p) likelihoodFit1(sigM1C(ii),sigP1C(ii),p,numTrials,[ePts3X(ii,1:numTrials)',ePts3Y(ii,1:numTrials)'],modelFitCircleNoise3(ii,1:numTrials),t,distTestEndpts3(ii,1:numTrials),fit1LookUpMat,maxDistAll,sigMmax,sigPmax,distFromTarget);
    
    [sigS1A(ii), minNLL1A(ii)] = fmincon(m1nLogLA,init,[],[],[],[],lb,ub,[],options); %find minimums - model 1 w/ data 1
    [sigS1B(ii), minNLL1B(ii)] = fmincon(m1nLogLB,init,[],[],[],[],lb,ub,[],options); %find minimums - model 1 w/ data 2
    [sigS1C(ii), minNLL1C(ii)] = fmincon(m1nLogLC,init,[],[],[],[],lb,ub,[],options); %find minimums - model 1 w/ data 3
    
end

estP1A = [sigM1A', sigP1A', sigS1A'];
estP1B = [sigM1B', sigP1B', sigS1B'];
estP1C = [sigM1C', sigP1C', sigS1C'];

save model1output_fixP.mat estP1A estP1B estP1C minNLL1A minNLL1B minNLL1C
