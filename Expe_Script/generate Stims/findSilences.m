function true_timepoints = findSilences(signal,sf,detection,detectionlevel,silsize)

%%% findSilences allows to extract silences from an audio signal.
%%%
%%% It takes SIGNAL in input with its sampling frequency ('SF').
%%% 
%%%
%%% Two methods of detection can be used. Silences can be taken between :
%%%
%%%     - mean+DETECTIONLEVEL*std and mean-DETECTIONLEVEL*std (type 'STD')
%%%
%%%     - 0+DETECTIONLEVEL and 0-DETECTIONLEVEL (type 'FIXED').
%%%
%%% For silences at exactly zero, type 'FIXED' and '0' in DETECTIONLEVEL.
%%%
%%% SILSIZE is an optional parameter that allows to capture only the
%%% silences exceeding a certain time. This time is in seconds.
%%%
%%% findSilences return a matrix containing all the silences as well as
%%% their duration. The duration is return in units which must be divided
%%% by the sampling rate in order to get it in seconds.


if nargin < 3
    detection = 'std';
end
if nargin < 4
    if strcmpi(detection,'std')==1
        detectionlevel = 0.01;
    elseif strcmpi(detection,'fixed')==1
        detectionlevel = 0.1;
    end        
end
if nargin < 5
    silsize = 49;
end

% Application of the detection method
if strcmpi(detection,'std')==1
    silences = (signal<(mean(signal)+detectionlevel*std(signal)))+(signal>(mean(signal)-detectionlevel*std(signal)));
elseif strcmpi(detection,'fixed')==1
    silences = (signal<=0+detectionlevel).*(signal>=0-detectionlevel);
end


% Reshape the input signal in a row
silences = reshape(silences,numel(silences),1);

% Find the time points (beginning and end) of silences
timepoints = find(diff([0;silences;0])~=0);

% Reshape the time points in two rows : [   beginning;
%                                           end       ]
timepoints = reshape(timepoints,2,numel(timepoints)/2);
if timepoints(end)>numel(signal)
    timepoints(end) = timepoints(end)-1;
end

% Find true silences (> 3 time steps)

timeThr = round(silsize/1000*sf);
true_timepoints = timepoints(:,(timepoints(2,:)-timepoints(1,:))>timeThr);
% 
% silmat = [true_timepoints(1,:);true_timepoints(2,:)-true_timepoints(1,:)];
% end