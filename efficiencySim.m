%% SIMULATION FOR MODEL EFFICENCY

%Tablet specs

radiusWac = linspace(0,69,150);

sigma_m = 20;
sigma_s = 4;
sigma_p = [18 44];
numTrials = 1;
numSims = 1;

xCen = 238;                                 %center of tablet x coordinate
yCen = 134;                                 %cetner of tablet y coordinate
t = [xCen yCen];                            %target location for all trials
tabSize = [268 476];                        %outside bounds of tablet space
r = [10*ones(1,5),linspace(10,0,length(6:150))]; %points earned for each circle size


%% SIMULATION FOR EFFICIENCY IDEAL AND PROSPECTIVE MODELS

R1 = length(radiusWac); %looping variable for circle size

for pp = 1:length(sigma_p)
    Rp = 1/sigma_p(pp)^2;
    
    for mm = 1:length(sigma_m)
        figure
        Rm = 1/sigma_m(mm)^2;
        ePt = t + randn(numTrials,2)*sigma_m(mm);       %randomized for all trials
        
        %Reach end points
        ePt1 = ePt(:,1);
        ePt2 = ePt(:,2);
        
        endDist = (sqrt((ePt1 - t(1)).^2 + (ePt2-t(2)).^2));
        
        
        for tt = 1:numTrials
            
            sensed = ePt(tt,:) + randn(1,2)*sigma_p(pp);
            
            %Posterior distribution variance
            mu_pos = (Rm/(Rp+Rm)).*t + (Rp/(Rp+Rm)).*sensed;
            var_pos = 1/(Rp+Rm);
            
            %Distance from posterior mean to target
            eucDistPos(tt) = sqrt((mu_pos(1) - t(1))^2 + (mu_pos(2)-t(2))^2);
            
            %Integrate over a circle to caclulate the probability of a
            %hit at a given distance
            ninterval = 1000;
            sd = sqrt(var_pos);
            c1 = 1/(sqrt(2*pi)*sd);
            c2 = 2*sd^2;
            c3 = sd*sqrt(2);
            
            phit1 = zeros(1,150);
            
            %Loop over all possible circle sizes to find probability of
            %points earned at a given distance
            for ii = 1:R1
                radius = radiusWac(ii);
                r2 = radius^2;
                dy = radius/ninterval;
                midy = (((1:ninterval) - .5)*dy).^2;
                maxx = sqrt(r2 - midy);
                yg = c1*exp(-midy/c2);
                phit1(ii) = dy*sum(yg.*(erf((maxx-eucDistPos(tt))/c3)-erf((-maxx-eucDistPos(tt))/c3)));
            end %ii
            
            gain1 = phit1.*r;  %probability of a hit vs points possible at each distance
            
            %Radius point for each distance with the highest gain
            temp1 = round(max(find(gain1 == max(gain1))));
            circ1 = radiusWac(temp1)+sigma_s*randn;
            
            if circ1 >= endDist(tt)
                gain1LookUpMat(tt,pp) = gain1(temp1);
            else
                gain1LookUpMat(tt,pp) = 0;
            end
            
            %PROSPECTIVE MODEL
            %loop over all tested sigma_m values
            phit3 = raylcdf(1:length(radiusWac),sigma_m(mm)); %CDF of a normal distribution centered on the target with variance of sigma_m
            
            gain3 = phit3.*r;                     %probability of a hit vs points possible at each distance
            
            temp3 = round(max(find(gain3 == max(gain3))));
            circ3 = radiusWac(temp3)+sigma_s*randn;
            
            if circ3 >= endDist(tt)
                gain3LookUpMat(tt,pp) = gain3(temp3);
            else
                gain3LookUpMat(tt,pp) = 0;
            end %if
            
            subplot(2,2,1)
            plot(radiusWac, phit1)
            hold on
            
            subplot(2,2,2)
            plot(radiusWac,phit3)
            hold on
            
            subplot(2,2,3)
            plot(radiusWac,gain1)
            xline(circ1);
            xline(endDist(tt),'r');
            xline(eucDistPos(tt),'c');
            hold on
            
            subplot(2,2,4)
            plot(radiusWac,gain3)
            xline(circ3);
            xline(endDist(tt),'r');
            xline(eucDistPos(tt),'c');
            hold on
        end %tt
        
    end %mm
end %pp

save efficiencySim.mat gain1LookUpMat gain3LookUpMat

%%

%mean point earnings
low1 = mean(gain1LookUpMat(:,1))
high1 = mean(gain1LookUpMat(:,2))
low3 = mean(gain3LookUpMat(:,1))
high3 = mean(gain3LookUpMat(:,2))

%efficiency ratios
ratioLow = low3/low1
ratioHigh = high3/high1

%Rate of endpoint capture
idealCap = sum(gain1LookUpMat~=0,1)/numTrials
prosCap = sum(gain3LookUpMat~=0,1)/numTrials