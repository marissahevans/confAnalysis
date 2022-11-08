%% CIRCLE LOOK UP TABLE FOR MODEL #3 - PROSPECTIVE MODEL - ONLY USES SIGMA_M

radiusWac = linspace(0,69,150);

xCen = 238;                                 %center of tablet x coordinate
yCen = 134;                                 %cetner of tablet y coordinate
t = [xCen yCen];                            %target location for all trials
tabSize = [268 476];                        %outside bounds of tablet space
matNan = nan(tabSize);                      %matrix of nans the size of the tablet
r = [10*ones(1,5),linspace(10,0,length(6:150))]; %points earned for each circle size
sigMmax = 1:100;                           %vector of test sigma_m values


% Find all distances to target on tablet
for ii = 1:size(matNan,1)
    for jj = 1:size(matNan,2)
        distFromTarget(ii,jj) = sqrt((jj - t(1))^2 + (ii-t(2))^2);
    end
end

maxDistAll = 0:max(round(distFromTarget(:))); %range of distance to max distance on the tablet.

fit3LookUpMat = nan(1,length(sigMmax)); %initalize vector of circle sizes based on sigma_m options

for mm = 1:length(sigMmax)              %loop over all tested sigma_m values
    phit = raylcdf(1:length(radiusWac),sigMmax(mm)); %CDF of a normal distribution centered on the target with variance of sigma_m
    
    gain = phit.*r;                     %probability of a hit vs points possible at each distance
    
    fit3LookUpMat(mm) = radiusWac(max(find(gain == max(gain)))); %circle size that matches the max expected gain
    
end

save model3LookUpTab.mat fit3LookUpMat