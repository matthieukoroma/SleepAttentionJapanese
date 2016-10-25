% function [SessionStructureLearning, SummaryLearning] = createSessionStructureLearningFinal_Sleep(SubjectID,gr,nEEGTraining,nTraining,nwake,nN2,nN3,nNoMemoryTest)
% Create a structure containing all trials informations. It requires the
% ID of the subject and its group as well as the number of trials of the
% sessions. In these n trials, nwake specify the number of trials where the
% subject will learn, ie wake or N2. The argument nNoMemoryTest specify the number of
% trials played continuously during the falling asleep period.

cd(stimPath)
rand('state',sum(100*clock));
load('lat_sqr_side')
sides = lat_sqr_side(gr,:);
tech='LKFS';
voices = {'f','m'};

% totTrials = totPaires/2;
% nNoMemoryTest = totTrials-(nwake+nsleep+nTraining)*2;

%stimPath=[pwd filesep '..' filesep 'Expe_Material'];
% stimPath=['D:\LSCPData\SleepAttentionJapanese\Expe_Folder\Expe_Material'];
% stimPath='C:\Data\SleepAttentionJapanese\ExpeMaterial';

[Codes,Categories,NamesJ,NameF,Freq,accrosslistpairs,attentionpairs] = load_stim_info(stimPath) ;
[CodesEEGT,TypesEEGT,ExcerptEEGT,VoiceEEGT] = load_EEGtraining_info(stimPath); 

% fid = fopen('C:\\Data\\SleepAttentionJapanese\\Expe_Folder\\Expe_Material\\CorrespondancesTextest.csv');
% Scaned = textscan(fid,'%s%s%s%s%s');

cd(stimPath)
load('trad_pairs.mat');
load('latin_square_EEGtraining.mat')
load('lat_sqr_behavTraining.mat')
load('lat_sqr.mat')
load('lat_sqr_falling.mat')
nTraining = nbehavTraining+nwakeTest;

Codes = cell2mat(trad_pairs(:,1));%Scaned{1}(2:end);%reference number of the stimulus
Categories = cell2mat(trad_pairs(:,4));%Scaned{4}(2:end);%category of the word
NamesJ = trad_pairs(:,3);%Scaned{2}(2:end);%names in japanese
NamesF = trad_pairs(:,2);%Scaned{3}(2:end);%names in french
Freq = cell2mat(trad_pairs(:,5));%Scaned{5}(2:end);%freq of the french name
acrosslistpairs = cell2mat(trad_pairs(:,7));%Scaned{7}(2:end);
attentionpairs = cell2mat(trad_pairs(:,6));%Scaned{6}(2:end);

% fid = fopen([stimPath, '\\CorrespondancesTextest.csv']);
% Scaned = textscan(fid,'%s%s%s%s%s');
% 
% Codes = Scaned{1}(2:end);%reference number of the stimulus
% Categories = Scaned{4}(2:end);%category of the word
% NamesJ = Scaned{2}(2:end);%names in japanese
% NamesF = Scaned{3}(2:end);%names in french
% Freq = Scaned{5}(2:end);%freq of the french name
%%
% Building falling Trials
%attended : trials attended side
%ignored : trials unattended side


permattended = randperm(nNoMemoryTest);
permcond = randperm(nNoMemoryTest);

fallingvoice = [repmat('f',nNoMemoryTest/2,1);repmat('h',nNoMemoryTest/2,1)];

ncond = numel(voices);
ncounterb = floor(nNoMemoryTest/ncond)*ncond;
nrest = nNoMemoryTest-ncounterb;

for nT=1:numel(voices)
    iinv = mod(nT,2)+1;
    subvoice = ncounterb/numel(voices);
    fallingattendedvoice((nT-1)*subvoice+1:nT*subvoice,1) = repmat(voices{nT},subvoice,1);
    fallingignoredvoice((nT-1)*subvoice+1:nT*subvoice,1) = repmat(voices{iinv},subvoice,1);
    fallingattendedside((nT-1)*subvoice+1:nT*subvoice,1) = repmat(sides(nT),subvoice,1);
    fallingignoredside((nT-1)*subvoice+1:nT*subvoice,1) = repmat(sides(iinv),subvoice,1);
end
%%
% Add trials which are not counterbalanced

for nT=1:nrest
    randvoice = randi(1:numel(voices));
    voicinv = mod(randvoice,2)+1;
    fallingattendedvoice(ncounterb+nT,1) = voices{randvoice};
    fallingignoredvoice(ncounterb+nT,1) = voices{voicinv};
    fallingattendedside(ncounterb+nT,1) = sides(randvoice);
    fallingignoredside(ncounterb+nT,1) = sides(voicinv);
end

%%
% Condition permutation

permcond = randperm(nNoMemoryTest);

fallingattendedvoice = fallingattendedvoice(permcond);
fallingignoredvoice =  fallingignoredvoice(permcond);
fallingattendedside = fallingattendedside(permcond);
fallingignoredside = fallingignoredside(permcond);


%%

% Building BehavTraining Trials

permattended = randperm(nbehavTraining);
permcond = randperm(nbehavTraining);

behavTrainingvoice = [repmat('f',nbehavTraining/2,1);repmat('h',nbehavTraining/2,1)];

ncond = numel(voices);
ncounterb = floor(nbehavTraining/ncond)*ncond;
nrest = nbehavTraining-ncounterb;

for nT=1:numel(voices)
    iinv = mod(nT,2)+1;
    subvoice = ncounterb/numel(voices);
    behavTrainingattendedvoice((nT-1)*subvoice+1:nT*subvoice,1) = repmat(voices{nT},subvoice,1);
    behavTrainingignoredvoice((nT-1)*subvoice+1:nT*subvoice,1) = repmat(voices{iinv},subvoice,1);
    behavTrainingattendedside((nT-1)*subvoice+1:nT*subvoice,1) = repmat(sides(nT),subvoice,1);
    behavTrainingignoredside((nT-1)*subvoice+1:nT*subvoice,1) = repmat(sides(iinv),subvoice,1);
end

%%
% Add trials which are not counterbalanced

for nT=1:nrest
    randvoice = randi(1:numel(voices));
    voicinv = mod(randvoice,2)+1;
    behavTrainingattendedvoice(ncounterb+nT,1) = voices{randvoice};
    behavTrainingignoredvoice(ncounterb+nT,1) = voices{voicinv};
    behavTrainingattendedside(ncounterb+nT,1) = sides(randvoice);
    behavTrainingignoredside(ncounterb+nT,1) = sides(voicinv);
end
%%
% Condition permutation

permcond = randperm(nbehavTraining);

behavTrainingattendedvoice = behavTrainingattendedvoice(permcond);
behavTrainingignoredvoice = behavTrainingignoredvoice(permcond);
behavTrainingattendedside = behavTrainingattendedside(permcond);
behavTrainingignoredside = behavTrainingignoredside(permcond);

%%
% Building wakeTest Trials

permattended = randperm(nwakeTest);
permcond = randperm(nwakeTest);

wakeTestvoice = [repmat('f',floor(nwakeTest/2),1);repmat('h',floor(nwakeTest/2),1)];

ncond = numel(voices);
ncounterb = floor(nwakeTest/ncond)*ncond;
nrest = nwakeTest-ncounterb;

for nT=1:numel(voices)
    iinv = mod(nT,2)+1;
    subvoice = ncounterb/numel(voices);
    wakeTestattendedvoice((nT-1)*subvoice+1:nT*subvoice,1) = repmat(voices{nT},subvoice,1);
    wakeTestignoredvoice((nT-1)*subvoice+1:nT*subvoice,1) = repmat(voices{iinv},subvoice,1);
    wakeTestattendedside((nT-1)*subvoice+1:nT*subvoice,1) = repmat(sides(nT),subvoice,1);
    wakeTestignoredside((nT-1)*subvoice+1:nT*subvoice,1) = repmat(sides(iinv),subvoice,1);
end

%%
% Add trials which are not counterbalanced

for nT=1:nrest
    randvoice = randi(1:numel(voices));
    voicinv = mod(randvoice,2)+1;
    wakeTestattendedvoice(ncounterb+nT,1) = voices{randvoice};
    wakeTestignoredvoice(ncounterb+nT,1) = voices{voicinv};
    wakeTestattendedside(ncounterb+nT,1) = sides(randvoice);
    wakeTestignoredside(ncounterb+nT,1) = sides(voicinv);
end
%%
% Condition permutation

permcond = randperm(nwakeTest);

wakeTestattendedvoice = wakeTestattendedvoice(permcond);
wakeTestignoredvoice = wakeTestignoredvoice(permcond);
wakeTestattendedside = wakeTestattendedside(permcond);
wakeTestignoredside = wakeTestignoredside(permcond);

%%
% Building Wake Trials

permattended = randperm(nwake);
permcond = randperm(nwake);

wakevoice = [repmat('f',nwake/2,1);repmat('h',nwake/2,1)];

ncond = numel(voices);
ncounterb = floor(nwake/ncond)*ncond;
nrest = nwake-ncounterb;

for nT=1:numel(voices)
    iinv = mod(nT,2)+1;
    subvoice = ncounterb/numel(voices);
    wakeattendedvoice((nT-1)*subvoice+1:nT*subvoice,1) = repmat(voices{nT},subvoice,1);
    wakeignoredvoice((nT-1)*subvoice+1:nT*subvoice,1) = repmat(voices{iinv},subvoice,1);
    wakeattendedside((nT-1)*subvoice+1:nT*subvoice,1) = repmat(sides(nT),subvoice,1);
    wakeignoredside((nT-1)*subvoice+1:nT*subvoice,1) = repmat(sides(iinv),subvoice,1);
end

%%
% Add trials which are not counterbalanced

for nT=1:nrest
    randvoice = randi(1:numel(voices));
    voicinv = mod(randvoice,2)+1;
    wakeattendedvoice(ncounterb+nT,1) = voices{randvoice};
    wakeignoredvoice(ncounterb+nT,1) = voices{voicinv};
    wakeattendedside(ncounterb+nT,1) = sides(randvoice);
    wakeignoredside(ncounterb+nT,1) = sides(voicinv);
end


%%
% Condition permutation

permcond = randperm(nwake);

wakeattendedvoice = wakeattendedvoice(permcond);
wakeignoredvoice = wakeignoredvoice(permcond);
wakeattendedside = wakeattendedside(permcond);
wakeignoredside = wakeignoredside(permcond);
%%
% Building Sleep N2 Trials

ncond = numel(voices);
ncounterb = floor(nN2/ncond)*ncond;
nrest = nN2-ncounterb;

sleepN2attendedvoice = [repmat('f',nN2/2,1);repmat('h',nN2/2,1)];

for nT=1:numel(voices)
    iinv = mod(nT,2)+1;
    subvoice = ncounterb/numel(voices);
    sleepN2attendedvoice((nT-1)*subvoice+1:nT*subvoice,1) = repmat(voices{nT},subvoice,1);
    sleepN2ignoredvoice((nT-1)*subvoice+1:nT*subvoice,1) = repmat(voices{iinv},subvoice,1);
    sleepN2attendedside((nT-1)*subvoice+1:nT*subvoice,1) = repmat(sides(nT),subvoice,1);
    sleepN2ignoredside((nT-1)*subvoice+1:nT*subvoice,1) = repmat(sides(iinv),subvoice,1);
end

%%
% Add trials which are not counterbalanced

for nT=1:nrest
    randvoice = randi(1:numel(voices));
    voicinv = mod(randvoice,2)+1;
    sleepN2attendedvoice(ncounterb+nT,1) = voices{randvoice};
    sleepN2ignoredvoice(ncounterb+nT,1) = voices{voicinv};
    sleepN2attendedside(ncounterb+nT,1) =  sides(randvoice);
    sleepN2ignoredside(ncounterb+nT,1) =  sides(voicinv);
end

%%
% Conditions permutation

permcond = randperm(nN2);

sleepN2attendedvoice = sleepN2attendedvoice(permcond);
sleepN2ignoredvoice = sleepN2ignoredvoice(permcond);
sleepN2attendedside = sleepN2attendedside(permcond);
sleepN2ignoredside = sleepN2ignoredside(permcond);

%%
% Building Sleep N3 Trials

ncond = numel(voices);
ncounterb = floor(nN3/ncond)*ncond;
nrest = nN3-ncounterb;

sleepN3attendedvoice = [repmat('f',nN3/2,1);repmat('h',nN3/2,1)];

for nT=1:numel(voices)
    iinv = mod(nT,2)+1;
    subvoice = ncounterb/numel(voices);
    sleepN3attendedvoice((nT-1)*subvoice+1:nT*subvoice,1) = repmat(voices{nT},subvoice,1);
    sleepN3ignoredvoice((nT-1)*subvoice+1:nT*subvoice,1) = repmat(voices{iinv},subvoice,1);
    sleepN3attendedside((nT-1)*subvoice+1:nT*subvoice,1) = repmat(sides(nT),subvoice,1);
    sleepN3ignoredside((nT-1)*subvoice+1:nT*subvoice,1) = repmat(sides(iinv),subvoice,1);
end

%%
% Add trials which are not counterbalanced

for nT=1:nrest
    randvoice = randi(1:numel(voices));
    voicinv = mod(randvoice,2)+1;
    sleepN3attendedvoice(ncounterb+nT,1) = voices{randvoice};
    sleepN3ignoredvoice(ncounterb+nT,1) = voices{voicinv};
    sleepN3attendedside(ncounterb+nT,1) = sides(randvoice);
    sleepN3ignoredside(ncounterb+nT,1) = sides(voicinv);
end

%%
% Conditions permutation

permcond = randperm(nN3);

sleepN3attendedvoice = sleepN3attendedvoice(permcond);
sleepN3ignoredvoice = sleepN3ignoredvoice(permcond);
sleepN3attendedside = sleepN3attendedside(permcond);
sleepN3ignoredside = sleepN3attendedside(permcond);
%%
totTrials = nEEGTraining+2*nTraining+nwake+nNoMemoryTest+nN2+nN3;
SummaryLearning = cell(totTrials,12);
SessionStructureLearning = cell(totTrials,1);

%%
%Building for EEG Train
TrainingEEGtrial = latin_square_EEGtraining(gr).order;

for T=1:nEEGTraining/2;
    stimchosen = TrainingEEGtrial(T);
    for i=1:2
        nT=(T-1)*2+i;
        Tcode(nT) = stimchosen;
        Texcerpt(nT) = i;
        Ttype{nT} = TypesEEGT{stimchosen*2};
        Tvoice(nT) = VoiceEEGT(stimchosen*2);
        attendedlist = struct('TypeCode', Ttype{nT},'Code', Tcode(nT),'Freq',Texcerpt(nT),'Order',NaN,'Name',NaN,'Side','Mono','Voice',Tvoice(nT),'Presented',0);
        ignoredlist = struct('TypeCode', Ttype{nT},'Code', Tcode(nT),'Freq',Texcerpt(nT),'Order',NaN,'Name',NaN,'Side','Mono','Voice',Tvoice(nT),'Presented',0);
        TrialStruct = struct('TrialType','EEGTrainingTrial','Condition','Mono','Number',nT,'FileName',sprintf('%s%s.mat',Ttype{nT},num2str(Texcerpt(nT))),'attended',attendedlist,'ignored',ignoredlist);
        SessionStructureLearning{nT} = TrialStruct;
        SummaryLearning{nT,1} = nT;
        SummaryLearning(nT,2:12) = {'TrainingEEG',Tcode(nT),'Mono',Tvoice(nT),NaN,NaN,NaN,NaN,NaN,NaN,0};
    end
end
%%
% Building of the structure experimental phase by exp phase
nTraining = nbehavTraining +nwakeTest;

len_list = length(lat_sqr_behavTraining(gr).attended);
permTrials = randperm(len_list);
behavTrainingattended = lat_sqr_behavTraining(gr).attended(permTrials(1:nTraining));
behavTrainingignored = lat_sqr_behavTraining(gr).ignored(permTrials(1:nTraining));
acode = behavTrainingattended;
icode = behavTrainingignored;


sentencesnames = cell(nTraining,2);
sentencestypes = cell(nTraining,2);
sentencescodes = cell(nTraining,2);
sentencesfreq = [Freq(acode),Freq(icode)];

for nT=1:nbehavTraining
    
    sentencescodes(nT,1:2) = {acode(nT),icode(nT)};
    sentencesnames(nT,1:2) = {['a',NamesF{acode(nT)}],['i',NamesF{icode(nT)}]};
    sentencestypes(nT,1:2) = {Categories(acode(nT)),Categories(icode(nT))};
    
    attendedlist = struct('TypeCode',sentencestypes{nT,1},'Code',sentencescodes{nT,1},'Freq',sentencesfreq(nT,1),'Order',NaN,'Name',sentencesnames{nT,1},'Side',behavTrainingattendedside(nT),'Voice',behavTrainingattendedvoice(nT),'Presented',0);
    ignoredlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'Freq',sentencesfreq(nT,2),'Order',NaN,'Name',sentencesnames{nT,2},'Side',behavTrainingignoredside(nT),'Voice',behavTrainingignoredvoice(nT),'Presented',0);
    TrialStruct = struct('TrialType','BehavTraining','Condition','Dicho','Number',nT,'FileName',sprintf('%s%s%s%s.mat',behavTrainingattendedside(nT),behavTrainingattendedvoice(nT),num2str(sentencescodes{nT,1}),tech),'attended',attendedlist,'ignored',ignoredlist);
    nTot = nEEGTraining+(nT-1)*2+1;
    SessionStructureLearning{nTot} = TrialStruct;
    SummaryLearning{nTot,1} = nTot;
    SummaryLearning(nTot,2:12) = {'BehavTraining',acode(nT),behavTrainingattendedside(nT),behavTrainingattendedvoice(nT),NaN,NaN,icode(nT),behavTrainingignoredside(nT),behavTrainingignoredvoice(nT),NaN,0};
    
    %dummy, to be replaced by a test trial
    TrialStruct = struct('TrialType','dummy','Condition','Dicho','Number',nT,'FileName',sprintf('%s%s%s%s.mat',behavTrainingattendedside(nT),behavTrainingattendedvoice(nT),num2str(sentencescodes{nT,1}),tech),'attended',attendedlist,'ignored',ignoredlist);
    SessionStructureLearning{nTot+1} = TrialStruct;
    SummaryLearning{nTot+1,1} = nTot;
    SummaryLearning(nTot+1,2:12) = {'dummy',acode(nT),behavTrainingattendedside(nT),behavTrainingattendedvoice(nT),NaN,NaN,icode(nT),behavTrainingignoredside(nT),behavTrainingignoredvoice(nT),NaN,0};
end

for nTwT=1:nwakeTest
    nT=nbehavTraining+nTwT;
    sentencescodes(nTwT,1:2) = {acode(nT),icode(nT)};
    sentencesnames(nTwT,1:2) = {['a',NamesF{acode(nT)}],['i',NamesF{icode(nT)}]};
    sentencestypes(nTwT,1:2) = {Categories(acode(nT)),Categories(icode(nT))};
    
    attendedlist = struct('TypeCode',sentencestypes{nTwT,1},'Code',sentencescodes{nTwT,1},'Freq',sentencesfreq(nTwT,1),'Order',NaN,'Name',sentencesnames{nTwT,1},'Side',wakeTestattendedside(nTwT),'Voice',wakeTestattendedvoice(nTwT),'Presented',0);
    ignoredlist = struct('TypeCode',sentencestypes{nTwT,2},'Code',sentencescodes{nTwT,2},'Freq',sentencesfreq(nTwT,2),'Order',NaN,'Name',sentencesnames{nTwT,2},'Side',wakeTestignoredside(nTwT),'Voice',wakeTestignoredvoice(nTwT),'Presented',0);
    TrialStruct = struct('TrialType','BehavTraining','Condition','Dicho','Number',nT,'FileName',sprintf('%s%s%s%s.mat',wakeTestattendedside(nTwT),wakeTestattendedvoice(nTwT),num2str(sentencescodes{nTwT,1}),tech),'attended',attendedlist,'ignored',ignoredlist);
    nTot = nEEGTraining+2*nbehavTraining+(nTwT-1)*2+1;
    SessionStructureLearning{nTot} = TrialStruct;
    SummaryLearning{nTot,1} = nTot;
    SummaryLearning(nTot,2:12) = {'BehavTraining',acode(nT),wakeTestattendedside(nTwT),wakeTestattendedvoice(nTwT),NaN,NaN,icode(nT),wakeTestignoredside(nTwT),wakeTestignoredvoice(nTwT),NaN,0};
    
    %dummy, to be replaced by a test trial
    TrialStruct = struct('TrialType','dummy','Condition','Dicho','Number',nT,'FileName',sprintf('%s%s%s%s.mat',wakeTestattendedside(nTwT),wakeTestattendedvoice(nTwT),num2str(sentencescodes{nTwT,1}),tech),'attended',attendedlist,'ignored',ignoredlist);
    SessionStructureLearning{nTot+1} = TrialStruct;
    SummaryLearning{nTot+1,1} = nTot;
    SummaryLearning(nTot+1,2:12) = {'dummy',acode(nT),wakeTestattendedside(nTwT),wakeTestattendedvoice(nTwT),NaN,NaN,icode(nT),wakeTestignoredside(nTwT),wakeTestignoredvoice(nTwT),NaN,0};
end
%%
% build sentences list and structure for each trial

permTrials = randperm(length(lat_sqr(gr).wake_ignored));
wakeattended = lat_sqr(gr).wake_attended(permTrials);
wakeignored = lat_sqr(gr).wake_ignored(permTrials);
acode = wakeattended;
icode = wakeignored;

sentencesnames = cell(nwake,2);
sentencestypes = cell(nwake,2);
sentencescodes = cell(nwake,2);
sentencesfreq = [Freq(acode),Freq(icode)];

for nT=1:nwake
    sentencescodes(nT,1:2) = {acode(nT),icode(nT)};
    sentencesnames(nT,1:2) = {['a' NamesF{acode(nT)}],['i',NamesF{icode(nT)}]};
    sentencestypes(nT,1:2) = {Categories(acode(nT)),Categories(icode(nT))};
    
    attendedlist = struct('TypeCode',sentencestypes{nT,1},'Code',sentencescodes{nT,1},'Freq',sentencesfreq(nT,1),'Order',NaN,'Name',sentencesnames{nT,1},'Side',wakeattendedside(nT),'Voice',wakeattendedvoice(nT),'Presented',0);
    ignoredlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'Freq',sentencesfreq(nT,2),'Order',NaN,'Name',sentencesnames{nT,2},'Side',wakeignoredside(nT),'Voice',wakeignoredvoice(nT),'Presented',0);
    
    nTtot = nEEGTraining+nT+2*nTraining;
    TrialStruct = struct('TrialType','WakeTrial','Condition','Dicho','Number',nTtot,'FileName',sprintf('%s%s%s%s.mat',wakeattendedside(nT),wakeattendedvoice(nT),num2str(sentencescodes{nT,1}),tech),'attended',attendedlist,'ignored',ignoredlist);
    SessionStructureLearning{nTtot} = TrialStruct;
    SummaryLearning{nTtot,1} = nTtot;
    SummaryLearning(nTtot,2:12) = {'Wake',acode(nT),wakeattendedside(nT),wakeattendedvoice(nT),NaN,NaN,icode(nT),wakeignoredside(nT),wakeignoredvoice(nT),NaN,0};

end
%%
listNoMemoryTestattended = [lat_sqr_falling(gr).attended];
listNoMemoryTestignored = [lat_sqr_falling(gr).ignored];
if nfalling>1 && nN3==0
    listNoMemoryTestattended = [listNoMemoryTestattended,lat_sqr_falling(gr).N3_attended];
    listNoMemoryTestignored = [listNoMemoryTestignored,lat_sqr_falling(gr).N3_ignored];
 end
permTrials = randperm(length(listNoMemoryTestattended));

NoMemoryTestattended = listNoMemoryTestattended(permTrials);
NoMemoryTestignored = listNoMemoryTestignored(permTrials);
acode = NoMemoryTestattended;
icode = NoMemoryTestignored;

sentencesnames = cell(nNoMemoryTest,2);
sentencestypes = cell(nNoMemoryTest,2);
sentencescodes = cell(nNoMemoryTest,2);
sentencesfreq = [Freq(acode),Freq(icode)];

for nT=1:nNoMemoryTest
    indf=find(attentionpairs==attentionpairs(acode(nT)));
    icode(nT) = setdiff(indf,acode(nT));
    sentencescodes(nT,1:2) = {acode(nT),icode(nT)};
    sentencesnames(nT,1:2) = {['a' NamesF{acode(nT)}],['i',NamesF{icode(nT)}]};
    sentencestypes(nT,1:2) = {Categories(acode(nT)),Categories(icode(nT))};
    
    attendedlist = struct('TypeCode',sentencestypes{nT,1},'Code',sentencescodes{nT,1},'Freq',sentencesfreq(nT,1),'Order',NaN,'Name',sentencesnames{nT,1},'Side',fallingattendedside(nT),'Voice',fallingattendedvoice(nT),'Presented',0);
    ignoredlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'Freq',sentencesfreq(nT,2),'Order',NaN,'Name',sentencesnames{nT,2},'Side',fallingignoredside(nT),'Voice',fallingignoredvoice(nT),'Presented',0);
    nTtot = nT+2*nTraining+nEEGTraining+nwake;
    
    TrialStruct = struct('TrialType','fallingTrial','Condition','Dicho','Number',nTtot,'FileName',sprintf('%s%s%s%s.mat',fallingattendedside(nT),fallingattendedvoice(nT),num2str(sentencescodes{nT,1}),tech),'attended',attendedlist,'ignored',ignoredlist);
    SessionStructureLearning{nTtot} = TrialStruct;
    SummaryLearning{nTtot,1} = nTtot;
    SummaryLearning(nTtot,2:12) = {'falling',acode(nT),fallingattendedside(nT),fallingattendedvoice(nT),NaN,NaN,icode(nT),fallingignoredside(nT),fallingignoredvoice(nT),NaN,0};
end


%%
% build sentences list and structure for each trial

permTrials = randperm(length(lat_sqr(gr).N2_attended));

sleepN2attended = lat_sqr(gr).N2_attended(permTrials);
sleepN2ignored = lat_sqr(gr).N2_ignored(permTrials);


acode = sleepN2attended;
icode = sleepN2ignored;
sentencesnames = cell(nN2,2);
sentencestypes = cell(nN2,2);
sentencescodes = cell(nN2,2);
sentencesfreq = [Freq(acode),Freq(icode)];


for nT=1:nN2
%     acode = sleepN2attendedvoice(nT);
%     icode = sleepN2attendedvoice(nT);

    ida = 0;
    idi = 0;

    sentencescodes(nT,1:2) = {acode(nT),icode(nT)};
    sentencesnames(nT,1:2) = {['a',NamesF{acode(nT)}],['i',NamesF{icode(nT)}]};
    sentencestypes(nT,1:2) = {Categories(acode(nT)),Categories(icode(nT))};
    
    attendedlist = struct('TypeCode',sentencestypes{nT,1},'Code',sentencescodes{nT,1},'Freq',sentencesfreq(nT,1),'Order',NaN,'Name',sentencesnames{nT,1},'Side',sleepN2attendedside(nT),'Voice',sleepN2attendedvoice(nT),'Presented',0);
    ignoredlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'Freq',sentencesfreq(nT,2),'Order',NaN,'Name',sentencesnames{nT,2},'Side',sleepN2ignoredside(nT),'Voice',sleepN2ignoredvoice(nT),'Presented',0);
    
    nTtot = nT+2*nTraining+nwake+nEEGTraining+nNoMemoryTest;
    TrialStruct = struct('TrialType','SleepN2Trial','Condition','Dicho','Number',nTtot,'FileName',sprintf('%s%s%s%s.mat',sleepN2attendedside(nT),sleepN2attendedvoice(nT),num2str(sentencescodes{nT,1}),tech),'attended',attendedlist,'ignored',ignoredlist);
    SessionStructureLearning{nTtot} = TrialStruct;
    SummaryLearning{nTtot,1} = nTtot;
    SummaryLearning(nTtot,2:12) = {'SleepN2',acode(nT),sleepN2attendedside(nT),sleepN2attendedvoice(nT),NaN,NaN,icode(nT),sleepN2ignoredside(nT),sleepN2ignoredvoice(nT),NaN,0};

end

%%
% build sentences list and structure for each trial

permTrials = randperm(length(lat_sqr(gr).N3_attended));

sleepN3attended = lat_sqr(gr).N3_attended(permTrials);
sleepN3ignored = lat_sqr(gr).N3_ignored(permTrials);

acode = sleepN3attended;
icode = sleepN3ignored;

sentencesnames = cell(nN3,2);
sentencestypes = cell(nN3,2);
sentencescodes = cell(nN3,2);
sentencesfreq = [Freq(acode),Freq(icode)];


for nT=1:nN3
    sentencescodes(nT,1:2) = {acode(nT),icode(nT)};
    sentencesnames(nT,1:2) = {['a' NamesF{acode(nT)}],['i',NamesF{icode(nT)}]};
    sentencestypes(nT,1:2) = {Categories(acode(nT)),Categories(icode(nT))};
    
    attendedlist = struct('TypeCode',sentencestypes{nT,1},'Code',sentencescodes{nT,1},'Freq',sentencesfreq(nT,1),'Order',NaN,'Name',sentencesnames{nT,1},'Side',sleepN3attendedside(nT),'Voice',sleepN3attendedvoice(nT),'Presented',0);
    ignoredlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'Freq',sentencesfreq(nT,2),'Order',NaN,'Name',sentencesnames{nT,2},'Side',sleepN3ignoredside(nT),'Voice',sleepN3ignoredvoice(nT),'Presented',0);
    
    nTtot = nT+2*nTraining+nwake+nNoMemoryTest+nEEGTraining+nN2;
    TrialStruct = struct('TrialType','SleepN3Trial','Condition','Dicho','Number',nTtot,'FileName',sprintf('%s%s%s%s.mat',sleepN3attendedside(nT),sleepN3attendedvoice(nT),num2str(sentencescodes{nT,1}),tech),'attended',attendedlist,'ignored',ignoredlist);
    SessionStructureLearning{nTtot} = TrialStruct;
    SummaryLearning{nTtot,1} = nTtot;
    SummaryLearning(nTtot,2:12) = {'SleepN3',acode(nT),sleepN3attendedside(nT),sleepN3attendedvoice(nT),NaN,NaN,icode(nT),sleepN3ignoredside(nT),sleepN3ignoredvoice(nT),NaN,0};
end

%%
filename = [stimPath, '\\SubjectSessions\\Subject', int2str(SubjectID)];
mkdir(filename);
save([filename, '\\Session'],'SessionStructureLearning')
save([filename, '\\SessionSummaryLearning.mat'],'SummaryLearning')
