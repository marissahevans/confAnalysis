%% Particpinat Data analysis Fit

participants = [{'RE'},{'PL'},{'HP'},{'LL'},{'ME'},{'MK'},{'ZL'},{'SX'},{'FH'},{'CS'},{'YK'},{'SM'},{'OX'},{'ET'},{'MP'},{'FL'}];

load model1LookUpTab %look up table for model 1
load model2LookUpTab %look up table for model 2
load model3LookUpTab %look up table for model 3

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
r = [10*ones(1,10),linspace(10,0,length(11:150))]; %points earned for each circle size
sigPmax = 1:100;                           %vector of test sigma_p values
sigMmax = 1:100;                           %vector of test sigma_m values

distFromTarget = matNan;

% Find all distances to target on tablet 
for nn = 1:size(matNan,1)
    for mm = 1:size(matNan,2)
        distFromTarget(nn,mm) = sqrt((mm - t(1))^2 + (nn-t(2))^2);
    end
end

maxDistAll = 0:max(distFromTarget(:)); %range of distance to max distance on the tablet. 
distFromTarget(distFromTarget > max(maxDistAll)) = max(maxDistAll);

numTrials = 300;

for ii = 1:length(participants)
    
subj = participants{ii};
path = sprintf('/Local/Users/marissa/Documents/Landy Lab/metaReachExperiment/data_metaReach/%s',subj);

load(sprintf('%s/%s_output.mat',path,subj)); %loads end points, circle sizes, and marginals for sigma m&p
if subj == 'FL'
endPoints(136,2) = 133;
end
distTestEndpts = sqrt((endPoints(:,1) - t(1)).^2 + (endPoints(:,2)-t(2)).^2);
endPtsFBdist = sqrt((endPtsFB(:,1) - t(1)).^2 + (endPtsFB(:,2)-t(2)).^2);
confCirc(confCirc <1) = 1;


%fmincon parameters 
lb = [1,1,1];                                               %lower bound
ub = [50,100,20];                          %upper bound using marginal from control
init = [20,40,10]; %, rand*(ub(3)-lb(3))+lb(3)];                                     %initiation
options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off'); %fmincon options

minLogLRetro = @(p) likelihoodFitRetro(p(1),p(2),p(3),numTrials,endPoints,confCirc,t,distTestEndpts,fit2LookUpMat,maxDistAll,sigPmax,distFromTarget,endPtsFBdist,contTar,contReach,contInd);
 
[estP2ret, minNLL2ret] = fmincon(minLogLRetro,init,[],[],[],[],lb,ub,[],options); %find minimums - model 2 

%Calculate the BIC
BIC = 2*minNLL2ret + 3*log(numTrials);

BICmin = BIC - min(BIC);

%Save output
filename = sprintf('%s_dataFittingOutput.mat',subj);
save(fullfile(path,filename), 'estP2ret', 'minNLL2ret','BIC')

subj

end 
