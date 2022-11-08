%% PARAMETERS USED BY ALL MODELS

%Tablet specs
parpool(5)

radiusWac = linspace(0,69,150);

xCen = 238;                                 %center of tablet x coordinate
yCen = 134;                                 %cetner of tablet y coordinate
t = [xCen yCen];                            %target location for all trials
tabSize = [268 476];                        %outside bounds of tablet space
matNan = nan(tabSize);                      %matrix of nans the size of the tablet
r = [10*ones(1,5),linspace(10,0,length(6:150))]; %points earned for each circle size
sigPmax = 1:100;                           %vector of test sigma_p values
sigMmax = 1:100;                           %vector of test sigma_m values


% Find all distances to target on tablet
for ii = 1:size(matNan,1)
    for jj = 1:size(matNan,2)
        distFromTarget(ii,jj) = sqrt((jj - t(1))^2 + (ii-t(2))^2);
    end
end

maxDistAll = 0:max(round(distFromTarget(:))); %range of distance to max distance on the tablet.

%% CIRCLE LOOK UP TABLE FOR MODEL #1 - IDEAL MODEL - USES BOTH SIGMA_M AND SIGMA_P

fit1LookUpMat = nan(length(maxDistAll),length(sigPmax),length(sigMmax)); %initalize variable for lookup table

L1 = length(maxDistAll); %Parpool looping variable for distance
P1 = length(sigPmax); %Parpool looping variable for sigma_p
M1 = length(sigMmax);%Parpool looping variable for sigma_m
R1 = length(radiusWac); %Parpool looping variable for circle size

parfor pp = 1:P1                            %Looping through all test values of sigma_p
    for mm = 1:M1                           %Looping through all test values of sigma_m
        for dist = 1:L1                      %loop over all distances
            
            %Posterior distribution variance
            Rm = 1/sigMmax(mm)^2;
            Rp = 1/sigPmax(pp)^2;
            mu_pos = (Rm/(Rp+Rm)).*t + (Rp/(Rp+Rm)).*[xCen+maxDistAll(dist),yCen];
            var_pos = 1/(Rp+Rm);
            
            %Distance from posterior mean to target
            eucDistPos = sqrt((mu_pos(1) - t(1))^2 + (mu_pos(2)-t(2))^2);

            %Integrate over a circle to caclulate the probability of a
            %hit at a given distance
            ninterval = 1000;
            sd = sqrt(var_pos);
            c1 = 1/(sqrt(2*pi)*sd);
            c2 = 2*sd^2;
            c3 = sd*sqrt(2);
            
            phit = zeros(1,150);
            
            %Loop over all possible circle sizes to find probability of
            %points earned at a given distance
            for ii = 1:R1
                radius = radiusWac(ii);
                r2 = radius^2;
                dy = radius/ninterval;
                midy = (((1:ninterval) - .5)*dy).^2;
                maxx = sqrt(r2 - midy);
                yg = c1*exp(-midy/c2);
                phit(ii) = dy*sum(yg.*(erf((maxx-eucDistPos)/c3)-erf((-maxx-eucDistPos)/c3)));
            end %ii
            
            gain = phit.*r;  %probability of a hit vs points possible at each distance
            
            %Radius point for each distance with the highest gain
            fit1LookUpMat(dist,pp,mm) = radiusWac(max(find(gain == max(gain))));
  
        end %dist
        
    end %mm
    pp
end %pp

save model1LookUpTab.mat fit1LookUpMat