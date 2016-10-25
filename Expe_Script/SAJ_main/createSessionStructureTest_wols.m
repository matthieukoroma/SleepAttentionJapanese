function [SessionStructureTest, SummaryTest] = createSessionStructureTest(SubjectID,SummaryLearning,listtests,oldStimLS)

% Create a structure containing all trials informations. It requires the
% ID of the subject and its group as well as the number of trials of the
% sessions. In these n trials, nwake specify the number of trials where the
% subject will learn, ie wake or N2. The argument nNoMemoryTest specify the number of
% trials played continuously during the falling asleep period.


%reset the rand function
rand('state',sum(100*clock));

%randomize order of presentation of tested items
order = {'p','s'}; %first and second
nwtest = {'n','o'};
%get the file for information about stimuli and which are paired acorss
%lists and attended sides in a trial

%stimPath=[pwd filesep '..' filesep 'Expe_Material'];
stimPath=['D:\LSCPData\SleepAttentionAllocation\Expe_Folder\Expe_Material'];

load([stimPath, '\\TranslationsCorrespondances3.mat']);

fid = fopen([stimPath, '\\CorrespondancesTextest.csv']);
Scaned = textscan(fid,'%s%s%s%s%s');

Codes = Scaned{1}(2:end);%reference number of the stimulus
Categories = Scaned{4}(2:end);%category of the word
NamesJ = Scaned{2}(2:end);%names in japanese
NamesF = Scaned{3}(2:end);%names in french
Freq = Scaned{5}(2:end);%freq of the french name
acrosslistpairs = Scaned{7}(2:end);
attentionpairs = Scaned{6}(2:end);

% stimToTest = latinsquare(gr).toTest;

%%get characteristics of the trials during learning to have the same voice
%%presentation
SummaryTest = cell(totTrials,12);
SessionStructureTest = cell(totTrials,1);
TrialsTypeSession = SummaryLearning{:}(2);

%% building the lists

%behavTrainingTrials
behavTrainingTrials = find(strcmp('behavTraining',TrialsTypeSession));
nbehavTraining = length(behavTrainingTrials);

%tests after awakening
if nargin<3
    %defining the lists used for the tests
    listtests = {{'sleepN2','sleepN3'},{'wake'}};
end

nbtest = length(listtests);
listtestTrials = cell(1,nbtest);

for testnr = 1:nbtest
    listtestTrials{testnr} = [];
    for cond=1:length(listtests{testnr})
        listtestTrials = [listtestTrials{testnr},strcmp(listtests{cond},TrialsTypeSession)];
    end
    listtestTrials{testnr} = unique(listtestTrials{testnr});
end
nlisttests = cellfun(@length,listtestTrials);%each trial has attended and ignored side but they are treated in parallel

% Building BehavTraining trial phase

ncond = numel(order);
ncounterb = floor(nbehavTraining/ncond)*ncond;
nrest = nbehavTraining-ncounterb;

for nT=1:numel(order)
    iinv = mod(nT,2)+1;
    suborder = ncounterb/numel(order);
    behavTrainingoldorder((nT-1)*suborder+1:nT*suborder,1) = repmat(order{nT},suborder,1);
    behavTrainingneworder((nT-1)*suborder+1:nT*suborder,1) = repmat(order{iinv},suborder,1);
end

%%
% Add trials which are not counterbalanced

for nT=1:nrest
    randorder = randi(1:numel(order));
    ordinv = mod(randorder,2)+1;
    behavTrainingoldorder(ncounterb+nT,1) = behavTrainingoldorder{randorder};
    behavTrainingneworder(ncounterb+nT,1) = behavTrainingneworder{ordinv};
end

%%
permcond = randperm(nbehavTraining);

attendedoldorder = behavTrainingoldorder(permcond);
attendedneworder = behavTrainingneworder(permcond);

attendedvoice = SummaryLearning{behavTrainingTrials}(5);
ocode = SummaryLearning{behavTrainingTrials}(3);
ncode = zeros(1,nbehavTraining);

sentencesnames = cell(nbehavTraining,3);
sentencestypes = cell(nbehavTraining,3);
sentencescodes = cell(nbehavTraining,3);
sentencesfreq = [Freq(acode);zeros(1,length(acode))];

for nT=1:nbehavTraining
    
    idoa = find(Codes == ocode(nT));
    index_acrosspairs = find(acrosslistpairs == acrosslistpairs(idoa));
    idn = setdiff(index_acrosspairs,idoa);
    ncode(nT) = Codes(idn);
    
    sentencescodes(nT,1:3) = {[3,ocode(nT)],[3,ncode(nT)],[2,ocode(nT)]};
    sentencesnames(nT,1:3) = {NamesF{idn},NamesF{idoa},NamesJ{idoa}};
    sentencestypes(nT,1:3) = {Categories{idoa},Categories{idn},Categories{idoa}};
    
    %first word
    oldlist = struct('TypeCode',sentencestypes{nT,1},'Code',sentencescodes{nT,1},'type',sentencesfreq(nT,1),'Name',sentencesnames{nT,1},'Side','Stereo','Voice',attendedvoice(nT),'tfTag',NaN,'tAmpTag',NaN,'Presented',0);
    %second word
    newlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'type',sentencesfreq(nT,2),'Name',sentencesnames{nT,2},'Side','Stereo','Voice',attendedvoice(nT),'tfTag',NaN,'tAmpTag',NaN,'Presented',0);
    %japanese word
    japlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'type',sentencesfreq(nT,2),'Name',sentencesnames{nT,2},'Side','Stereo','Voice',attendedvoice(nT),'tfTag',NaN,'tAmpTag',NaN,'Presented',0);
    
    TrialStruct = struct('TrialType','behavTrainingTrial','Condition','Dicho','Number',nT,'FileName',sprintf('%st%s%s%sj%s%s.mat',attendedvoice(nT),sentencesnames{nT,1},behavTrainingoldorder(nT),attendedvoice(nT),sentencesnames{nT,2},behavTrainingneworder(nT),attendedvoice(nT),sentencesnames{nT,3},behavTrainingoldorder(nT)),'old',oldlist,'new',newlist,'jap',japlist);
    
    nTtot = nT;
    SessionStructureTest{nTtot} = TrialStruct;
    SummaryTest{nTtot,1} = nTtot;
    SummaryTest(nTtot,2:12) = {'behavTraining',sentencesnames{nT,1},attendedvoice(nT),behavTrainingoldorder(nT),attendedvoice(nT),sentencesnames{nT,2},behavTrainingneworder(nT),attendedvoice(nT),sentencesnames{nT,3},behavTrainingoldorder(nT)};
end

%% loop on the test phases

for tnr = 1:nbtests
    % Building trials for test during the wake learning phase
    
    ntrials = nlisttests(tnr);
    
    ncond = numel(order)*numel(nwtest);
    ncounterb = floor(ntrials/ncond)*ncond;%???/2 because it will be done for old primed trials and new primed trials
    nrest = ntrials-ncounterb;
    
    for nT=1:numel(order)
        iinv = mod(nT,2)+1;
        suborder = ncounterb/numel(order);
        oldorder((nT-1)*suborder+1:nT*suborder,1) = repmat(order{nT},suborder,1);
        neworder((nT-1)*suborder+1:nT*suborder,1) = repmat(order{iinv},suborder,1);
        for j = 1:numel(nwtest)
            jinv = mod(k,2)+1;
            subnwtest = subtag/numel(nwtest);
            oldnw(1+(nT-1)*subnwtest+(j-1)*suborder:(nT-1)*subnwtest+j*suborder,1) = repmat(nwtest{j},subnwtest,1);
            newnw(1+(nT-1)*subnwtest+(j-1)*suborder:(nT-1)*subnwtest+j*suborder,1) = repmat(nwtest{jinv},subnwtest,1);
        end
    end
    
    %%
    % Add trials which are not counterbalanced
    
    for nT=1:nrest
        randorder = randi(1:numel(order));
        ordinv = mod(randorder,2)+1;
        oldorder(ncounterb+nT,1) = order{randorder};
        neworder(ncounterb+nT,1) = order{ordinv};
        
        randnw = randi(1:numel(nwtest));
        nwinv = mod(randtag,2)+1;
        waketaletag(ncounterb+nT,1) = fTags{randtag};
        wakejabtag(ncounterb+nT,1) = fTags{taginv};
    end
    
    %%
    %attended old
    permcond = randperm(ntrials);    
    aoorder1 = aoorder(permcond);%ao : attended old
    aoorder2 = anorder(permcond);%an : attended new
    
    %attended new
    permcond = randperm(ntrials);    
    anorder1 = aoorder(permcond);
    anorder2 = anorder(permcond);
    
    %ignored old
    permcond = randperm(ntrials);
    ioorder1 = ioorder(permcond);
    ioorder2 = inorder(permcond);
   
    %ignored old    
    permcond = randperm(ntrials);
    inorder1 = ioorder(permcond);
    inorder2 = inorder(permcond);
    
%     % assemble the vector to make the random permutation
%     toPermVec = [aoorder,aoorder,ioorder,ioorder;anorder]
    %%
    % Condition permutation
    
    %%
    
    %randomisation at the end 
    testTrials = listtestTrials{tnr};
    attendedvoice = SummaryLearning{testTrialsrand}(5);
    ignoredvoice = SummaryLearning{testTrialsrand}(5);
    oacode = SummaryLearning{testTrials}(3);
    oicode = SummaryLearning{testTrials}(8);
    nacode = zeros(1,ntrials);
    nicode = zeros(1,ntrials);

    sentencesnames = cell(ntrials,3);
    sentencestypes = cell(ntrials,3);
    sentencescodes = cell(ntrials,3);
    sentencesfreq = [Freq(ocode);Freq(oicode)];
    
    aocounting = 1;
    ancounting = 1;
    iocounting = 1;
    incounting = 1;
    
    for nT=1:ntrials
        
        idoa = find(Codes == oacode(nT));
        index_acrosspairs = find(acrosslistpairs == acrosslistpairs(idoa));
        idna = setdiff(index_acrosspairs,idoa);
        ncode(nT) = Codes(idna);
        
        if ismember(oacode(nT),oldStimLS)
            oorder = aoorder1(aocounting);
            norder = aoorder2(aocounting);
            aocounting = aocounting+1;
            %match the japanese word
            japorder = oorder;
            idj = idoa;
            sentencescodes(nT,1:3) = {[3,ocode(nT)],[3,ncode(nT)],[2,ocode(nT)]};
       else  
            oorder = aoorder1(ancounting);
            norder = aoorder2(ancounting);
            ancounting = ancounting+1;
            japorder = norder;
            idj = idna;
            sentencescodes(nT,1:3) = {[3,ocode(nT)],[3,ncode(nT)],[2,ncode(nT)]};
       end
        sentencesnames(nT,1:3) = {NamesF{idn},NamesF{idoa},NamesJ{idj}};
        sentencestypes(nT,1:3) = {Categories{idoa},Categories{idna},Categories{idj}};
        
        %first word
        oldlist = struct('TypeCode',sentencestypes{nT,1},'Code',sentencescodes{nT,1},'type',sentencesfreq(nT,1),'Name',sentencesnames{nT,1},'Side','Stereo','Voice',attendedvoice(nT),'tfTag',NaN,'tAmpTag',NaN,'Presented',0);
        %second word
        newlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'type',sentencesfreq(nT,2),'Name',sentencesnames{nT,2},'Side','Stereo','Voice',attendedvoice(nT),'tfTag',NaN,'tAmpTag',NaN,'Presented',0);
        %japanese word
        japlist = struct('TypeCode',sentencestypes{nT,3},'Code',sentencescodes{nT,3},'type',sentencesfreq(nT,3),'Name',sentencesnames{nT,2},'Side','Stereo','Voice',attendedvoice(nT),'tfTag',NaN,'tAmpTag',NaN,'Presented',0);
        
        TrialStruct = struct('TrialType',['Test_',num2str(nrtest), '_across'],'Condition','Dicho','Number',nT,'FileName',sprintf('%st%s%s%sj%s%s.mat',onjapendedvoice(nT),sentencesnames{nT,1},oorder(nT),attendedvoice(nT),sentencesnames{nT,2},norder(nT),attendedvoice(nT),sentencesnames{nT,3},japorder(nT)),'old',oldlist,'new',newlist,'jap',japlist);
        
        nTtot = nT;
        SessionStructureTest{nTtot} = TrialStruct;
        SummaryTest{nTtot,1} = nTtot;
        SummaryTest(nTtot,2:12) = {['Test_',num2str(nrtest), '_across'],sentencesnames{nT,1},onjapendedvoice(nT),sleepacrossoldorder(nT),onjapendedvoice(nT),sentencesnames{nT,2},sleepacrossneworder(nT),onjapendedvoice(nT),sentencesnames{nT,3},sleepacrossoldorder(nT)};
        
        idoi = find(Codes == oicode(nT));
        index_acrosspairs = find(acrosslistpairs == acrosslistpairs(idoi));
        idni = setdiff(index_acrosspairs,idoi);
        ncode(nT) = Codes(idni);
        if ismember(oacode(nT),oldStimLS)
            oorder = ioorder1(iocounting);
            norder = ioorder2(iocounting);
            aocounting = iocounting+1;
            japorder = oorder;
            idj = idoi;
            sentencescodes(nT,1:3) = {[3,ocode(nT)],[3,ncode(nT)],[2,ocode(nT)]};
        else
            oorder = inorder1(incounting);
            norder = inorder2(incounting);
            ancounting = incounting+1;
            japorder = norder;
            idj = idni;
            sentencescodes(nT,1:3) = {[3,ocode(nT)],[3,ncode(nT)],[2,ncode(nT)]};
        end
        
        sentencesnames(nT,1:3) = {NamesF{idoi},NamesF{idni},NamesJ{idj}};
        sentencestypes(nT,1:3) = {Categories{idoi},Categories{idni},Categories{idj}};
        
        %first word
        oldlist = struct('TypeCode',sentencestypes{nT,1},'Code',sentencescodes{nT,1},'type',sentencesfreq(nT,1),'Name',sentencesnames{nT,1},'Side','Stereo','Voice',ignoredvoice(nT),'tfTag',NaN,'tAmpTag',NaN,'Presented',0);
        %second word
        newlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'type',sentencesfreq(nT,2),'Name',sentencesnames{nT,2},'Side','Stereo','Voice',ignoredvoice(nT),'tfTag',NaN,'tAmpTag',NaN,'Presented',0);
        %japanese word
        japlist = struct('TypeCode',sentencestypes{nT,3},'Code',sentencescodes{nT,3},'type',sentencesfreq(nT,3),'Name',sentencesnames{nT,2},'Side','Stereo','Voice',ignoredvoice(nT),'tfTag',NaN,'tAmpTag',NaN,'Presented',0);
        
        TrialStruct = struct('TrialType','WithinTrial','Condition','Dicho','Number',nT,'FileName',sprintf('%st%s%s%sj%s%s.mat',ignoredvoice(nT),sentencesnames{nT,1},oorder(nT),ignoredvoice(nT),sentencesnames{nT,2},norder(nT),ignoredvoice(nT),sentencesnames{nT,3},japorder(nT)),'old',oldlist,'new',newlist,'jap',japlist);
        
        nTtot = nT;
        SessionStructureTest{nTtot} = TrialStruct;
        SummaryTest{nTtot,1} = nTtot;
        SummaryTest(nTtot,2:12) = {['Test_',num2str(nrtest), '_across'],sentencesnames{nT,1},onjapendedvoice(nT),oorder(nT),ignoredvoice(nT),sentencesnames{nT,2},norder(nT),ignoredvoice(nT),sentencesnames{nT,3},japorder(nT)};
    end
    %%
    % build sentences list and structure for each trial
    %within list
    
    ncond = numel(order);
    ncounterb = floor(ntrials/2/ncond)*ncond;%/2 because it will be done for old primed trials and new primed trials
    nrest = ntrials-ncounterb;
    
    for nT=1:numel(order)
        iinv = mod(nT,2)+1;
        suborder = ncounterb/numel(order);
        attendedoldorder((nT-1)*suborder+1:nT*suborder,1) = repmat(order{nT},suborder,1);
        attendedneworder((nT-1)*suborder+1:nT*suborder,1) = repmat(order{iinv},suborder,1);
    end
    
    %%
    % Add trials which are not counterbalanced
    
    for nT=1:nrest
        randorder = randi(1:numel(order));
        ordinv = mod(randorder,2)+1;
        attendedoldorder(ncounterb+nT,1) = order{randorder};
        attendedneworder(ncounterb+nT,1) = order{ordinv};
    end
    
    %%
    %attended old
    permcond = randperm(ntrials);    
    aoorder1 = aoorder(permcond);%ao : attended old
    aoorder2 = anorder(permcond);%an : attended new
        
    %ignored old
    permcond = randperm(ntrials);
    ioorder1 = ioorder(permcond);
    ioorder2 = inorder(permcond);
   
    
    for nT=1:ntrials1
        ocode = find(lettertypes==oldtypecode{nT}(1))*100+str2num(oldtypecode{nT}(2:3));
        ite = 1;
        ida = 0;
        idoi = 0;
        for j=1:numel(Codes)
            if strcmp(oldtypecode{nT},char(Codes(j,:)))
                ida = ite;
            end
            if strcmp(newtypecode{nT},char(Codes(j,:)))
                idoi = ite;
            end
            ite = ite+1;
        end
        sentencescodes(nT,1:2) = {ocode,ncode};
        sentencesnames(nT,1:2) = {Names{ida},['i',Names{idoi}]};
        sentencestypes(nT,1:2) = {Types{ida},Types{idoi}};
        
        oldlist = struct('TypeCode',oldtypecode{nT},'Code',sentencescodes{nT,1},'type',sentencestypes{nT,1},'Name',sentencesnames{nT,1},'Side',oldside(nT),'order',aoorder(nT),'tfTag',NaN,'tAmpTag',NaN);
        newlist = struct('TypeCode',newtypecode{nT},'Code',sentencescodes{nT,2},'type',sentencestypes{nT,2},'Name',sentencesnames{nT,2},'Side',newside(nT),'order',anorder(nT),'tfTag',NaN,'tAmpTag',NaN);
        TrialStruct = struct('TrialType','SleepTrial','Condition','Dicho','Number',nT+ntraining,'FileName',sprintf('%st%s%s%sj%s%s.mat',sleepN3oldside(nT),oldtypecode{nT},sleepN3oldorder(nT),oldside(nT),newtypecode{nT},sleepN3neworder(nT),sleepN3newside(nT)),'old',oldlist,'new',newlist);
        
        nTtot = nT+nBehavTraining+nSleepN2;
        SessionStructure{nTtot} = TrialStruct;
        SummaryTest{nTtot,1} = nTtot;
        SummaryTest(nTtot,2:12) = {'SleepN3',oldtypecode{nT,1},sleepN3oldside(nT),sleepN3oldorder(nT),newtypecode{nT,1},NaN,NaN,sleepN3newside(nT),sleepN3neworder(nT),NaN,NaN};
    end
    
    filename = [stimPath, '\\SubjectSessions\\Subject', int2str(SubjectID)];
    mkdir(filename);
    save([filename, '\\Session'],'SessionStructure')
    save([filename, '\\SessionSummary.mat'],'Summary')
    
end