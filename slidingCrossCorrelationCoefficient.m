function [r,lag] = slidingCrossCorrelationCoefficient(X,Y,truncateBy)
% [R,LAG] = SLIDINGCROSSCORRELATIONCOEFFICIENT(X,Y) is a function that
% computes a correlation cofficient for each possible time-shifted
% comparison of sequences X and Y (if X = Y, this is a normalised
% autocorrelation). In this analysis, X remains fixed, and Y is adjusted by
% the amount described in the 'lag' output. The r value is the correlation
% coefficient, which must be between -1 and 1. This value is computed by
% mean-subtracting and then scaling by the inverse sd of each overlapping
% chunk of X and Y for each individual time-shift. The final input,
% 'truncateBy' simply limits the specified number of time-shift
% calculations by not performing those with a small number of overlapping
% elements.
%
% Created by SML Sept 2017

% Defaults:
if nargin < 3
    truncateBy = 0;
end

% Transpose if row vectors:
if size(X,1) == 1
    X = X';
end
if size(Y,1) == 1
    Y = Y';
end

% Check the inputs
assert(truncateBy >=0, 'Please specify a postive integer for the variable truncate by!')
assert(size(X,2)==size(Y,2),'X and Y do not have the same number of signals!')
assert(size(X,1)==size(Y,1),'This code is written with the assumption that X and Y have the same length!')
nSignals = size(X,2);
N = size(X,1);

% Lag axis:
lag = -(N-truncateBy-1):(N-truncateBy-1);
nLag = length(lag);
r = NaN(nLag,nSignals);

% Truncate the signals and compute each correlation:
for ss = 1:nSignals % EACH signal  
    for mm = 1:nLag % EACH lag
        lagBy = lag(mm);
        nOverlap = N - abs(lagBy);
        if lagBy < 0
            X_t = X(1:nOverlap,ss);
            Y_t = Y(end-nOverlap+1:end,ss);
        elseif lagBy > 0
            X_t = X(end-nOverlap+1:end,ss);
            Y_t = Y(1:nOverlap,ss);
        else
            X_t = X(:,ss);
            Y_t = Y(:,ss);
        end
        r(mm,ss) = computeSingleCoeff(X_t,Y_t);
    end
end

end

function [r] = computeSingleCoeff(A,B)
% Takes truncated segments and computes correlation.

A = A - nanmean(A);
B = B - nanmean(B);
r = dot(A/norm(A),B/norm(B));

end