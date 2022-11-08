%% LOOK UP TABLE GENERATOR FOR MODEL 2 - RETROSPECTIVE MODEL

parpool(1)

%Tablet specs
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

%% CIRCLE LOOK UP TABLE FOR MODEL #2 - RETROPECTIVE MODEL - MAINLY USES SIGMA_P

fit2LookUpMat = nan(length(maxDistAll),length(sigPmax)); %initalize variable for lookup table

L1 = length(maxDistAll);
P1 = length(sigPmax);
R1 = length(radiusWac);

parfor dist = 1:L1                  %loop over all x axis points
        for pp = 1:P1          %loop through all possible values of sigma_p
            
            %Integrate over a circle to find probability of a hit at each
            %distance
            ninterval = 1000;
            sd = sigPmax(pp);
            c1 = 1/(sqrt(2*pi)*sd);
            c2 = 2*sd^2;
            c3 = sd*sqrt(2);
            
            phit = zeros(1,150);
            
            for ii = 1:R1
                radius = radiusWac(ii);
                r2 = radius^2;
                dy = radius/ninterval;
                midy = (((1:ninterval) - .5)*dy).^2;
                maxx = sqrt(r2 - midy);
                yg = c1*exp(-midy/c2);
                phit(ii) = dy*sum(yg.*(erf((maxx-maxDistAll(dist))/c3)-erf((-maxx-maxDistAll(dist))/c3)));
            end
            
            gain = phit.*r; %probability of a hit vs points possible at each distance
            
            %Radius point for each distance with the highest gain
            fit2LookUpMat(dist,pp) = radiusWac(max(find(gain == max(gain))));
            
        end %pp
    dist
end %dist

save model2LookUpTab.mat fit2LookUpMat