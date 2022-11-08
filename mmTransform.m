%% Whats happening in this script -
%This script cleans, organizes and centers the raw data. First it
%transforms all points into mm from Wacom tablet coordinates (A.U.), then
%rotates the endpoints for trials where a confidence judgement was made and
%centers them all on the same target location.

%For the control experiment a quick fit is done of the motor error and
%propioceptive error marginasls as an idiot check before running the full
%model to make sure nothing wacky is going on there.


%A file is saved with the confidence reports, points earned, points
%possible, shifted X,Y endpoint coordinates, original enpoint coordnates,
%original target coordinates for all confidence trials, and control
%experiment endpoints, and control experiment indicated points, plus the
%marginals for reach error and proprioceptive error.


%% Setting up tablet specifics
radiusWac = linspace(0,69,150);

xCen = 238;                                 %center of tablet x coordinate
yCen = 134;                                 %cetner of tablet y coordinate
startX = 238;                               %reach start point X
startY = 43;                                %reach start point Y
t = [xCen yCen];                            %target location for all trials
tabSize = [268 476];                        %outside bounds of tablet space
matNan = nan(tabSize);                      %matrix of nans the size of the tablet

%% Find all distances to target on tablet
for ii = 1:size(matNan,1)
    for jj = 1:size(matNan,2)
        distFromTarget(ii,jj) = sqrt((jj - t(1))^2 + (ii-t(2))^2);
    end
end

maxDistAll = 0:max(round(distFromTarget(:))); %range of distance to max distance on the tablet.
testDist = max(maxDistAll);                 %distance range radius from target tested for s
possibleCircles = 0:testDist;               %based on the maximum radius possible that returns points

%% Load participant data - looping over sessions
subjAll = [{'PL'},{'LL'},{'ZL'},{'FH'},{'FL'},{'RE'},{'HP'},{'ME'},{'MK'},{'SX'},{'CS'},{'YK'},{'SM'},{'OX'},{'ET'},{'MP'}];
for ss = 1:length(subjAll)
    subj = subjAll{ss};
    path = sprintf('/Local/Users/marissa/Documents/Landy Lab/metaReachExperiment/data_metaReach/%s',subj);
    
    if subj == 'FL'
        R1 = 'FL_metaReach_exp_S1_2022-02-04_16-07_controlresults';
        D1 = 'FL_metaReach_exp_S1_2022-02-04_16-07_dispInfo';
        R2 = 'FL_metaReach_exp_S2_2022-02-09_16-00_results';
        D2 = 'FL_metaReach_exp_S2_2022-02-09_16-00_dispInfo';
        R3 = 'FL_metaReach_exp_S3_2022-03-01_13-26_results';
        D3 = 'FL_metaReach_exp_S3_2022-03-01_13-26_dispInfo';
        R4 = 'FL_metaReach_exp_S4_2022-03-04_17-29_results';
        D4 = 'FL_metaReach_exp_S4_2022-03-04_17-29_dispInfo';
        
    elseif subj == 'PL'
        R1 = 'PL_metaReach_exp_S1_2021-10-18_16-20_controlresults';
        D1 = 'PL_metaReach_exp_S2_2021-10-20_13-54_dispInfo';
        R2 = 'PL_resultsMat_S2';
        D2 = 'PL_metaReach_exp_S2_2021-10-20_13-54_dispInfo';
        R3 = 'PL_metaReach_exp_S3_2021-10-20_15-38_results';
        D3 = 'PL_metaReach_exp_S3_2021-10-20_15-38_dispInfo';
        R4 = 'PL_metaReach_exp_S4_2021-10-22_14-04_results';
        D4 = 'PL_metaReach_exp_S4_2021-10-22_14-04_dispInfo';
        
    elseif subj == 'LL'
        R1 = 'LL_metaReach_exp_S1_2021-10-27_17-34_controlresults';
        D1 = 'LL_metaReach_exp_S1_2021-10-27_17-34_dispInfo';
        R2 = 'LL_metaReach_exp_S2_2021-10-28_17-50_results';
        D2 = 'LL_metaReach_exp_S2_2021-10-28_17-50_dispInfo';
        R3 = 'LL_metaReach_exp_S3_2021-10-29_15-34_results';
        D3 = 'LL_metaReach_exp_S3_2021-10-29_15-34_dispInfo';
        R4 = 'LL_metaReach_exp_S4_2021-10-29_16-46_results';
        D4 = 'LL_metaReach_exp_S4_2021-10-29_16-46_dispInfo';
        
    elseif subj == 'ZL'
        R1 = 'ZL_metaReach_exp_S1_2021-11-02_17-54_controlresults';
        D1 = 'ZL_metaReach_exp_S1_2021-11-02_17-54_dispInfo';
        R2 = 'ZL_metaReach_exp_S2_2021-11-05_15-56_results';
        D2 = 'ZL_metaReach_exp_S2_2021-11-05_15-56_dispInfo';
        R3 = 'ZL_metaReach_exp_S3_2021-12-02_20-36_results';
        D3 = 'ZL_metaReach_exp_S3_2021-12-02_20-36_dispInfo';
        R4 = 'ZL_metaReach_exp_S4_2021-12-03_16-15_results';
        D4 = 'ZL_metaReach_exp_S4_2021-12-03_16-15_dispInfo';
        
    elseif subj == 'FH'
        R1 = 'FH_metaReach_exp_S1_2022-01-25_11-07_controlresults';
        D1 = 'FH_metaReach_exp_S1_2022-01-25_11-07_dispInfo';
        R2 = 'FH_metaReach_exp_S2_2022-01-27_11-06_results';
        D2 = 'FH_metaReach_exp_S2_2022-01-27_11-06_dispInfo';
        R3 = 'FH_metaReach_exp_S3_2022-02-01_17-31_results';
        D3 = 'FH_metaReach_exp_S3_2022-02-01_17-31_dispInfo';
        R4 = 'FH_metaReach_exp_S4_2022-02-03_10-49_results';
        D4 = 'FH_metaReach_exp_S4_2022-02-03_10-49_dispInfo';
        
    elseif subj == 'RE'
        R1 = 'RE_metaReach_exp_S1_2021-10-21_11-28_controlresults';
        D1 = 'RE_metaReach_exp_S1_2021-10-21_11-28_dispInfo';
        R2 = 'RE_metaReach_exp_S2_2021-10-21_15-50_results';
        D2 = 'RE_metaReach_exp_S2_2021-10-21_15-50_dispInfo';
        R3 = 'RE_metaReach_exp_S3_2021-10-27_14-09_results';
        D3 = 'RE_metaReach_exp_S3_2021-10-27_14-09_dispInfo';
        R4 = 'RE_metaReach_exp_S4_2021-10-28_15-33_results';
        D4 = 'RE_metaReach_exp_S4_2021-10-28_15-33_dispInfo';
        
    elseif subj == 'HP'
        R1 = 'HP_metaReach_exp_S1_2021-10-18_12-12_controlresults';
        D1 = 'HP_metaReach_exp_S2_2021-10-19_11-10_dispInfo';
        R2 = 'HP_metaReach_exp_S2_2021-10-19_11-10_results';
        D2 = 'HP_metaReach_exp_S2_2021-10-19_11-10_dispInfo';
        R3 = 'HP_metaReach_exp_S3_2021-10-19_12-20_results';
        D3 = 'HP_metaReach_exp_S3_2021-10-19_12-20_dispInfo';
        R4 = 'HP_metaReach_exp_S4_2021-10-19_15-35_results';
        D4 = 'HP_metaReach_exp_S4_2021-10-19_15-35_dispInfo';
        
    elseif subj == 'ME'
        R1 = 'ME_metaReach_exp_S1_2021-10-26_12-52_controlresults';
        D1 = 'ME_metaReach_exp_S1_2021-10-26_12-52_dispInfo';
        R2 = 'ME_metaReach_exp_S2_2021-10-26_17-10_results';
        D2 = 'ME_metaReach_exp_S2_2021-10-26_17-10_dispInfo';
        R3 = 'ME_metaReach_exp_S3_2021-10-27_11-59_results';
        D3 = 'ME_metaReach_exp_S3_2021-10-27_11-59_dispInfo';
        R4 = 'ME_metaReach_exp_S4_2021-11-02_12-41_results';
        D4 = 'ME_metaReach_exp_S4_2021-11-02_12-41_dispInfo';
        
    elseif subj == 'MK'
        R1 = 'MK_metaReach_exp_S1_2021-10-27_16-04_controlresults';
        D1 = 'MK_metaReach_exp_S1_2021-10-27_16-04_dispInfo';
        R2 = 'MK_metaReach_exp_S2_2021-10-28_10-43_results';
        D2 = 'MK_metaReach_exp_S2_2021-10-28_10-43_dispInfo';
        R3 = 'MK_metaReach_exp_S3_2021-11-01_17-11_results';
        D3 = 'MK_metaReach_exp_S3_2021-11-01_17-11_dispInfo';
        R4 = 'MK_metaReach_exp_S4_2021-11-03_11-46_results';
        D4 = 'MK_metaReach_exp_S4_2021-11-03_11-46_dispInfo';
        
    elseif subj == 'SX'
        R1 = 'SX_metaReach_exp_S1_2022-01-18_14-01_controlresults';
        D1 = 'SX_metaReach_exp_S1_2022-01-18_14-01_dispInfo';
        R2 = 'SX_metaReach_exp_S2_2022-01-27_12-52_results';
        D2 = 'SX_metaReach_exp_S2_2022-01-27_12-52_dispInfo';
        R3 = 'SX_metaReach_exp_S3_2022-01-31_13-39_results';
        D3 = 'SX_metaReach_exp_S3_2022-01-31_13-39_dispInfo';
        R4 = 'SX_metaReach_exp_S4_2022-02-02_14-16_results';
        D4 = 'SX_metaReach_exp_S4_2022-02-02_14-16_dispInfo';
        
    elseif subj == 'CS'
        R1 = 'CS_metaReach_exp_S1_2022-01-31_14-30_controlresults';
        D1 = 'CS_metaReach_exp_S1_2022-01-31_14-30_dispInfo';
        R2 = 'CS_metaReach_exp_S2_2022-02-02_11-02_results';
        D2 = 'CS_metaReach_exp_S2_2022-02-02_11-02_dispInfo';
        R3 = 'CS_metaReach_exp_S3_2022-02-02_12-03_results';
        D3 = 'CS_metaReach_exp_S3_2022-02-02_12-03_dispInfo';
        R4 = 'CS_metaReach_exp_S4_2022-02-08_11-54_results';
        D4 = 'CS_metaReach_exp_S4_2022-02-08_11-54_dispInfo';
        
    elseif subj == 'YK'
        R1 = 'YK_metaReach_exp_S1_2022-02-09_13-52_controlresults';
        D1 = 'YK_metaReach_exp_S1_2022-02-09_13-52_dispInfo';
        R2 = 'YK_metaReach_exp_S2_2022-02-15_16-33_results';
        D2 = 'YK_metaReach_exp_S2_2022-02-15_16-33_dispInfo';
        R3 = 'YK_metaReach_exp_S3_2022-02-18_11-09_results';
        D3 = 'YK_metaReach_exp_S3_2022-02-18_11-09_dispInfo';
        R4 = 'YK_metaReach_exp_S4_2022-02-24_12-57_results';
        D4 = 'YK_metaReach_exp_S4_2022-02-24_12-57_dispInfo';
        
    elseif subj == 'SM'
        R1 = 'SM_metaReach_exp_S1_2022-02-17_14-03_controlresults';
        D1 = 'SM_metaReach_exp_S1_2022-02-17_14-03_dispInfo';
        R2 = 'SM_metaReach_exp_S2_2022-02-22_11-58_results';
        D2 = 'SM_metaReach_exp_S2_2022-02-22_11-58_dispInfo';
        R3 = 'SM_metaReach_exp_S3_2022-02-23_13-52_results';
        D3 = 'SM_metaReach_exp_S3_2022-02-23_13-52_dispInfo';
        R4 = 'SM_metaReach_exp_S4_2022-02-28_11-38_results';
        D4 = 'SM_metaReach_exp_S4_2022-02-28_11-38_dispInfo';
        
    elseif subj == 'OX'
        R1 = 'OX_metaReach_exp_S1_2022-02-08_14-59_controlresults';
        D1 = 'OX_metaReach_exp_S1_2022-02-08_14-59_dispInfo';
        R2 = 'OX_metaReach_exp_S2_2022-02-16_16-15_results';
        D2 = 'OX_metaReach_exp_S2_2022-02-16_16-15_dispInfo';
        R3 = 'OX_metaReach_exp_S3_2022-02-28_14-41_results';
        D3 = 'OX_metaReach_exp_S3_2022-02-28_14-41_dispInfo';
        R4 = 'OX_metaReach_exp_S4_2022-03-01_12-00_results';
        D4 = 'OX_metaReach_exp_S4_2022-03-01_12-00_dispInfo';
        
    elseif subj == 'ET'
        R1 = 'ET_metaReach_exp_S1_2022-02-17_16-58_controlresults';
        D1 = 'ET_metaReach_exp_S1_2022-02-17_16-58_dispInfo';
        R2 = 'ET_metaReach_exp_S2_2022-02-23_17-38_results';
        D2 = 'ET_metaReach_exp_S2_2022-02-23_17-38_dispInfo';
        R3 = 'ET_metaReach_exp_S3_2022-02-28_13-48_results';
        D3 = 'ET_metaReach_exp_S3_2022-02-28_13-48_dispInfo';
        D4 = 'ET_metaReach_exp_S4_2022-03-02_14-45_dispInfo';
        R4 = 'ET_metaReach_exp_S4_2022-03-02_14-45_results';
        
    elseif subj == 'MP'
        R1 = 'MP_metaReach_exp_S1_2022-02-24_15-25_controlresults';
        D1 = 'MP_metaReach_exp_S1_2022-02-24_15-25_dispInfo';
        R2 = 'MP_metaReach_exp_S2_2022-02-28_12-35_results';
        D2 = 'MP_metaReach_exp_S2_2022-02-28_12-35_dispInfo';
        R3 = 'MP_metaReach_exp_S3_2022-03-01_14-49_results';
        D3 = 'MP_metaReach_exp_S3_2022-03-01_14-49_dispInfo';
        R4 = 'MP_metaReach_exp_S4_2022-03-03_11-02_results';
        D4 = 'MP_metaReach_exp_S4_2022-03-03_11-02_dispInfo';
    end
    
    confCirc = [];
    endPoints = [];
    endPtsFB = [];
    rawEndpoints = [];
    rawTarg = [];
    
    for ii = 1:4 %number of sessions
        if ii == 1
            %control session
            load(sprintf('%s',path,'/',R1,'.mat'));
            load(sprintf('%s',path,'/',D1,'.mat'));
            
        elseif ii == 2
            %experimental sessions
            load(sprintf('%s',path,'/',R2,'.mat'));
            load(sprintf('%s',path,'/',D2,'.mat'));
        elseif ii == 3
            load(sprintf('%s',path,'/',R3,'.mat'));
            load(sprintf('%s',path,'/',D3,'.mat'));
        elseif ii == 4
            load(sprintf('%s',path,'/',R4,'.mat'));
            load(sprintf('%s',path,'/',D4,'.mat'));
        end
        %% Creating transforms from calibration
        
        %Inverse transform to get from pixel to tablet space
        invtform = invert(displayInfo.tform);           %Inverse transformation from participant's calibration
        
        %Creating a matrix for the calibration points in tablet space
        wacX = displayInfo.calibration{1,6};
        wacY = displayInfo.calibration{1,7};
        wac = [wacX' wacY'];
        
        %Calibration point locations in mm measurments from edge of tablet active
        %area
        mmWacSpace = [23,244; 23,134; 23,27; 238,244; 238,134; 238,27; 446, 244; 446,134; 446,27;];
        
        %affine transform from tablet space to mm
        M = affine_least_square(wac(3,1),wac(3,2),wac(5,1),wac(5,2),wac(6,1),wac(6,2), mmWacSpace(3,1),mmWacSpace(3,2),mmWacSpace(5,1),mmWacSpace(5,2),mmWacSpace(6,1),mmWacSpace(6,2));
        tformMM = affine2d(M');
        
        
        %% CONTROL EXPERIMENT (session 1)
        if ii == 1
            
            %converting target to mm space
            [target1(1), target1(2)] = (transformPointsForward(invtform,contResultsMat.targetLoc(1),contResultsMat.targetLoc(2)));
            [target1(1), target1(2)] = (transformPointsForward(tformMM,target1(1),target1(2)));
            
            p_vec = 1:100;      % vector of possible proprioceptive noise values
            m_vec = 1:100;      % vector of possible motor noise values
            
            % Reach end points converted to mm space - known to the experimenter
            if subj == 'MP'
                reaches = [contResultsMat.wacEndPoint(1:240,1),contResultsMat.wacEndPoint(1:240,2)];
            else
                reaches = [contResultsMat.wacEndPoint(61:end,1),contResultsMat.wacEndPoint(61:end,2)];
            end
            [reaches(:,1), reaches(:,2)] = transformPointsForward(tformMM,reaches(:,1),reaches(:,2));
            
            
            % Participant-reported endpoint converted to mm space - known to the experimenter
            mouse = [contResultsMat.mouseEndPt(61:end,1),contResultsMat.mouseEndPt(61:end,2)];
            [indicated(:,1), indicated(:,2)] = (transformPointsForward(invtform,mouse(:,1),mouse(:,2)));
            [indicated(:,1), indicated(:,2)] = (transformPointsForward(tformMM,indicated(:,1),indicated(:,2)));
            
            % Simultaneously estimate sigma_m and sigma_p by ML
            
            LL=zeros(length(m_vec),length(p_vec));
            LLigivene=zeros(length(m_vec),length(p_vec));
            LLe=zeros(length(m_vec),1);
            
            for vv = 1:length(m_vec)        %loop over all sigma_m options
                RmTemp = 1/m_vec(vv)^2;
                
                % log likelihood of sigma_m given target and endpoint:
                LLe(vv) = log2disotnormal(reaches,target1,m_vec(vv));
                for jj = 1:length(p_vec)    %loop over all sigma_p options
                    RpTemp = 1/p_vec(jj)^2;
                    meanigivene = (RpTemp/(RpTemp+RmTemp))*reaches + ...
                        (RmTemp/(RpTemp+RmTemp))*target1;
                    SDigivene = (RpTemp/(RpTemp+RmTemp))*p_vec(jj);
                    
                    % log likelihood of sigma_p given sensed location and endpoint:
                    LLigivene(vv,jj) = log2disotnormal(indicated,meanigivene,SDigivene);
                    
                    % log likelihood of the sigma_p/sigma_m pair:
                    LL(vv,jj) = LLe(vv) + LLigivene(vv,jj);
                end
            end
            
            % treat LL like a log posterior (i.e., treat prior as flat over the grid)
            % and calculate marginals. First add a constant to all LL values so that
            % the maximum is one (to minimize underflows) and normalize afterward.
            
            NormPost = exp(LL - max(LL(:)));
            
            %Motor Error Marginal
            mmarg = sum(NormPost,2);
            mmarg = mmarg/sum(mmarg);
            sigMmarg = m_vec(find(mmarg == max(mmarg)));
            
            %Proprioception Marginal
            pmarg = sum(NormPost,1);
            pmarg = pmarg/sum(pmarg);
            sigPmarg = p_vec(find(pmarg == max(pmarg)));
            
            %output variables
            contReach = reaches;
            contTar = target1;
            contInd = indicated;
            
            %% CONFIDENCE EXPERIMENT (session 2-4)
        else
            
            %% Converting participant data to mm space
            %create index of confidence judgement trials
            confTrial = resultsMat.confRad ~= 0;        %select only trials where a confidence judgment was made
            conf = resultsMat.confRad(confTrial);      %Confidence reports from data
            pointsEarned = resultsMat.pointsEarned;
            pointsPoss = resultsMat.pointsPossible;
            
            for jj = 1:length(conf)
                temp = conf(jj);
                [temp1, temp2] = transformPointsForward(invtform,temp+displayInfo.xCenter,displayInfo.yCenter);
                dataConf(jj,:) = (transformPointsForward(tformMM,temp1,temp2))-xCen;
            end
            
            endPts = resultsMat.wacEndPoint;
            [endPts(:,1), endPts(:,2)] = transformPointsForward(tformMM,endPts(:,1),endPts(:,2));
            
            [target(:,1), target(:,2)] = (transformPointsForward(invtform,resultsMat.targetLoc(:,1),resultsMat.targetLoc(:,2))); %transform target locations into wacom space
            [target(:,1), target(:,2)] = transformPointsForward(tformMM,target(:,1),target(:,2));
            
            %Rotate points to center coordinates
            for tt = 1:length(target)
                
                %Finding the angle from the starting point, the center of the tablet,
                %and the target location
                cosTheta = dot(target(tt,:)-[startX,startY],[xCen,yCen]-[startX,startY])/((sqrt((startX - target(tt,1))^2 + (startY-target(tt,2))^2))*(sqrt((startX - xCen)^2 + (startY-yCen)^2)));
                theta = real(acosd(cosTheta));
                
                %choose if rotation is clockwise or counter clockwise
                if target(tt,1) <= 238
                    R = [cosd(theta) sind(theta); -sind(theta) cosd(theta)];
                else
                    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
                end
                
                % Rotate target points
                targetShift = R*(target(tt,:)'-[startX; startY]) + [startX; startY];
                
                %Collapse across y space
                yShift = targetShift(2)-yCen;
                
                %rotate data points by target angle then shift by difference in Y
                endPtsShift(tt,:) = R*(endPts(tt,:)'-[startX; startY]) + [startX; startY] - [0; yShift];
                
            end
            
            confCirc = [confCirc; dataConf];
            endPoints = [endPoints; endPtsShift(confTrial,:)];
            endPtsFB = [endPtsFB; endPtsShift(~confTrial,:)];
            rawEndpoints = [rawEndpoints; endPts];
            rawTarg = [rawTarg; target];
            
            
            
        end
    end

% figure
% scatter(endPoints(:,1),endPoints(:,2))
% ylim([0 tabSize(1)])
% xlim([0 tabSize(2)])

filename = sprintf('%s_output.mat',subj);
save(fullfile(path,filename), 'confCirc', 'endPoints','endPtsFB','sigMmarg','sigPmarg','rawEndpoints','rawTarg','contTar','contReach','contInd','pointsEarned','pointsPoss')

clear indicated
end
% direct computation of log(p(data|mu,sigma))
% data is Nx2
% mu is 1x2 or Nx2
% sigma is a scalar SD of the isotropic 2-d Gaussian

function ll = log2disotnormal(data,mu,sigma)

centered = data - mu;
N = size(data,1);
ll = -N*log(2*pi) - 2*N*log(sigma) - (sum(centered(:).^2)/(2*sigma^2));
end