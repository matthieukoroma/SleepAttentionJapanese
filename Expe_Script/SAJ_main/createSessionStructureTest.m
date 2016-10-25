% function [SessionStructureTestTest, SummaryTest] = createSessionStructureTest(SubjectID,SummaryLearning,listtests,oldStimLS)

% Create a structure containing all trials informations. It requires the
% ID of the subject and its group as well as the number of trials of the
% sessions. In these n trials, nwake specify the number of trials where the
% subject will learn, ie wake or N2. The argument nNoMemoryTest specify the number of
% trials played continuously during the falling asleep period.


%reset the rand function
rand('state',sum(100*clock));

%randomize order of presentation of tested items
order = [1,2]; %first and second
tech = 'LKFS';
%get the file for information about stimuli and which are paired acorss
%lists and attended sides in a trial

%stimPath=[pwd filesep '..' filesep 'Expe_Material'];
% stimPath=['D:\LSCPData\SleepAttentionAllocation\Expe_Folder\Expe_Material'];
% stimPath='C:\Data\SleepAttentionJapanese\Expe_Material';

[Codes,Categories,NamesJ,NameF,Freq,accrosslistpairs,attentionpairs] = load_stim_info(stimPath) ;
% oldStimLS = [cell2mat(SummaryLearning(:,3));cell2mat(SummaryLearning(:,8))];
%% building the lists

%behavTrainingTrials
% behavTrainingTrials = find(strcmp('behavTraining',TrialsTypeSession));
% nbehavTraining = length(behavTrainingTrials);

%tests after awakening
%defining the lists used for the tests
listtests = {{'SleepN2','SleepN3'},{'Wake'}};
merge_tests = [1,1,2];

nlisttests = 16;%cellfun(@,cellfun(@length,listtests),'uni',0);

%%get characteristics of the trials during learning to have the same voice
%%presentation
totTrials = sum(nlisttests.*1.5)+nbehavTraining+nwakeTest;
SummaryTest = cell(totTrials,12);
SessionStructureTest = cell(totTrials,1);
TrialsTypeSession = SummaryLearning(:,2);


nbtest = length(listtests);
listtestTrials = cell(max(cellfun(@length,listtests)),nbtest);

for testnr = 1:nbtest
    listtestTrials{testnr} = [];
    for cond=1:length(listtests{testnr})
        stim_file = cellfun(@(x) strcmp(x,listtests{testnr}(cond)),TrialsTypeSession);
        index = find(stim_file);
        listtestTrials{testnr,cond} = index;
    end
    listtestTrials{testnr,cond} = unique(listtestTrials{testnr,cond});
end
nlisttests = cellfun(@length,listtestTrials);%each trial has attended and ignored side but they are treated in parallel

% Building nWakeTest trial phase

ncond = numel(order);
ncounterb = floor(nwakeTest/ncond)*ncond;
nrest = nwakeTest-ncounterb;

for nT=1:numel(order)
    iinv = mod(nT,2)+1;
    suborder = ncounterb/numel(order);
    wakeTestoldorder((nT-1)*suborder+1:nT*suborder,1) = repmat(order(nT),suborder,1);
    wakeTestneworder((nT-1)*suborder+1:nT*suborder,1) = repmat(order(iinv),suborder,1);
end

%%
% Add trials which are not counterbalanced

for nT=1:nrest
    randorder = randi(1:numel(order));
    ordinv = mod(randorder,2)+1;
    wakeTestoldorder(ncounterb+nT,1) = wakeTestoldorder(randorder);
    wakeTestneworder(ncounterb+nT,1) = wakeTestneworder(ordinv);
end
%%

% Building BehavTraining trial phase

ncond = numel(order);
ncounterb = floor(nbehavTraining/ncond)*ncond;
nrest = nbehavTraining-ncounterb;

for nT=1:numel(order)
    iinv = mod(nT,2)+1;
    suborder = ncounterb/numel(order);
    behavTrainingoldorder((nT-1)*suborder+1:nT*suborder,1) = repmat(order(nT),suborder,1);
    behavTrainingneworder((nT-1)*suborder+1:nT*suborder,1) = repmat(order(iinv),suborder,1);
end

%%
% Add trials which are not counterbalanced

for nT=1:nrest
    randorder = randi(1:numel(order));
    ordinv = mod(randorder,2)+1;
    behavTrainingoldorder(ncounterb+nT,1) = behavTrainingoldorder(randorder);
    behavTrainingneworder(ncounterb+nT,1) = behavTrainingneworder(ordinv);
end
%%
stim_file = cellfun(@(x) strcmp(x,'BehavTraining'),TrialsTypeSession);
index1 = find(stim_file);
behavTrainingTrials =index1;

Presented = cell2mat(SummaryLearning(behavTrainingTrials,12));
permcond = randperm(nbehavTraining);
attendedoldorder = behavTrainingoldorder(permcond);
attendedneworder = behavTrainingneworder(permcond);


attendedvoice = cell2mat(SummaryLearning(behavTrainingTrials,5));
attendedside = cell2mat(SummaryLearning(behavTrainingTrials,4));
ignoredvoice = cell2mat(SummaryLearning(behavTrainingTrials,10));
ignoredside = cell2mat(SummaryLearning(behavTrainingTrials,9));

ocode = cell2mat(SummaryLearning(behavTrainingTrials,3));
ncode = zeros(1,nTraining);

sentencesnames = cell(nTraining,3);
sentencestypes = cell(nTraining,3);
sentencescodes = cell(nTraining,3);
sentencesfreq = [Freq(ocode),zeros(length(ocode),1)];

for nT=1:nbehavTraining
    
    idoa = find(Codes == ocode(nT));
    index_acrosspairs = find(acrosslistpairs == acrosslistpairs(idoa));
    idn = setdiff(index_acrosspairs,idoa);
    ncode(nT) = Codes(idn);
    sentencesfreq(nT,2) = ncode(nT);
    sentencescodes(nT,1:3) = {[ocode(nT)],[ncode(nT)],[ocode(nT)]};
    sentencesnames(nT,1:3) = {NamesF{idn},NamesF{idoa},NamesJ{idoa}};
    sentencestypes(nT,1:3) = {Categories(idoa),Categories(idn),Categories(idoa)};
    [dummy,first] = min([attendedoldorder(nT),attendedneworder(nT)]);second=setdiff(1:2,first);
    goodresp=[1,-1];
    
    %first word
    oldlist = struct('TypeCode',sentencestypes{nT,1},'Code',sentencescodes{nT,1},'Freq',sentencesfreq(nT,1),'Name',sentencesnames{nT,1},'Order',attendedoldorder(nT),'Side','Stereo','Voice',attendedvoice(nT),'Presented',goodresp(1));
    %second word
    newlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'Freq',sentencesfreq(nT,2),'Name',sentencesnames{nT,2},'Order',attendedneworder(nT),'Side','Stereo','Voice',attendedvoice(nT),'Presented',goodresp(2));
    %japanese word
    japlist = struct('TypeCode',sentencestypes{nT,3},'Code',sentencescodes{nT,3},'Freq',sentencesfreq(nT,1),'Name',sentencesnames{nT,3},'Order',attendedoldorder(nT),'Side','Stereo','Voice',attendedvoice(nT),'Presented',0);
    
    TrialStruct = struct('TrialType','behavTraining','Condition','attended','Number',nT,'FileName',sprintf('%s%st%so%s%s%s.mat','a',attendedvoice(nT),num2str(sentencescodes{nT,3}),num2str(sentencescodes{nT,first}),num2str(sentencescodes{nT,second}),tech),'old',oldlist,'new',newlist,'jap',japlist);
    
    nTtot = nT;
    SessionStructureTest{nTtot} = TrialStruct;
    SummaryTest{nTtot,1} = nTtot;
    SummaryTest(nTtot,2:12) = {'behavTraining',ocode(nT),attendedside(nT),attendedvoice(nT),behavTrainingoldorder(nT),ncode(nT),ignoredside(nT),'attended',behavTrainingneworder(nT),ocode(nT),Presented(nT)};
end

%%
permcond = randperm(nwakeTest);
oldwakeTestoldorder = wakeTestoldorder(permcond);
newwakeTestoldorder = wakeTestneworder(permcond);

for nTwT=1:nwakeTest
    nT=nbehavTraining+nTwT;
    
    idoa = find(Codes == ocode(nT));
    index_acrosspairs = find(acrosslistpairs == acrosslistpairs(idoa));
    idn = setdiff(index_acrosspairs,idoa);
    ncode(nT) = Codes(idn);
    sentencesfreq(nTwT,2:3) = [Freq(ncode(nT)),sentencesfreq(nTwT,1)];
    sentencescodes(nTwT,1:3) = {[ocode(nT)],[ncode(nT)],[ocode(nT)]};
    sentencesnames(nTwT,1:3) = {NamesF{idn},NamesF{idoa},NamesJ{idoa}};
    sentencestypes(nTwT,1:3) = {Categories(idoa),Categories(idn),Categories(idoa)};
    [dummy,first] = min([oldwakeTestoldorder(nTwT),newwakeTestoldorder(nTwT)]);second=setdiff(1:2,first);
    goodresp=[1,-1];
    
    %first word
    oldlist = struct('TypeCode',sentencestypes{nTwT,1},'Code',sentencescodes{nTwT,1},'Freq',sentencesfreq(nTwT,1),'Name',sentencesnames{nTwT,1},'Order',oldwakeTestoldorder(nTwT),'Side','Stereo','Voice',attendedvoice(nT),'Presented',goodresp(1));
    %second word
    newlist = struct('TypeCode',sentencestypes{nTwT,2},'Code',sentencescodes{nTwT,2},'Freq',sentencesfreq(nTwT,2),'Name',sentencesnames{nTwT,2},'Order',newwakeTestoldorder(nTwT),'Side','Stereo','Voice',attendedvoice(nT),'Presented',goodresp(2));
    %japanese word
    japlist = struct('TypeCode',sentencestypes{nTwT,3},'Code',sentencescodes{nTwT,3},'Freq',sentencesfreq(nTwT,1),'Name',sentencesnames{nTwT,3},'Order',oldwakeTestoldorder(nTwT),'Side','Stereo','Voice',attendedvoice(nT),'Presented',0);
    
    TrialStruct = struct('TrialType','behavTraining','Condition','attended','Number',nT,'FileName',sprintf('%s%st%so%s%s%s.mat','a',attendedvoice(nT),num2str(sentencescodes{nTwT,3}),num2str(sentencescodes{nTwT,first}),num2str(sentencescodes{nTwT,second}),tech),'old',oldlist,'new',newlist,'jap',japlist);
    
    nTtot = nT;
    SessionStructureTest{nTtot} = TrialStruct;
    SummaryTest{nTtot,1} = nTtot;
    SummaryTest(nTtot,2:12) = {'behavTraining',ocode(nT),attendedside(nT),attendedvoice(nT),oldwakeTestoldorder(nTwT),ncode(nT),ignoredside(nT),'attended',newwakeTestoldorder(nTwT),ocode(nT),Presented(nT)};
end
%% loop on the test phases

nlast = 0;

for tnr = 1:nbtest
    
    trialstomerge_within = [];
    trialstomerge_across = [];
    nfirst = nlast;
    
    for cond=1:length(listtests{tnr})
        % Building trials for test during the wake learning phase
        
        title_list = listtests{tnr}{cond};
        ntrials = nlisttests(tnr,cond);
        
        ncond = numel(order);
        ncounterb = floor(ntrials/4/ncond)*ncond;%/2 because it will be done for old primed trials and new primed trials
        nrest = ntrials-ncounterb;
        
        for nT=1:numel(order)
            iinv = mod(nT,2)+1;
            suborder = ncounterb/numel(order);
            aoorder((nT-1)*suborder+1:nT*suborder,1) = repmat(order(nT),suborder,1);
            anorder((nT-1)*suborder+1:nT*suborder,1) = repmat(order(iinv),suborder,1);
        end
        
        %%
        % Add trials which are not counterbalanced
        
        for nT=1:nrest
            randorder = randi(1:numel(order));
            ordinv = mod(randorder,2)+1;
            aoorder(ncounterb+nT,1) = order(randorder);
            anorder(ncounterb+nT,1) = order(ordinv);
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
        ioorder1 = aoorder(permcond);
        ioorder2 = anorder(permcond);
        
        %ignored old
        permcond = randperm(ntrials);
        inorder1 = aoorder(permcond);
        inorder2 = anorder(permcond);
        
        %     % assemble the vector to make the random permutation
        %     toPermVec = [aoorder,aoorder,ioorder,ioorder;anorder]
        %%
        % Condition permutation
        
        %%
        
        testTrials = listtestTrials{tnr,cond};
        Presented = cell2mat(SummaryLearning(testTrials,12));
        attendedvoice = cell2mat(SummaryLearning(testTrials,5));
        ignoredvoice = cell2mat(SummaryLearning(testTrials,10));
        oldattcode = cell2mat(SummaryLearning(testTrials,3));
        oldigncode = cell2mat(SummaryLearning(testTrials,8));
        nacode = zeros(1,ntrials);
        nicode = zeros(1,ntrials);
        
        randTrials = randperm(length(oldattcode));
        
        oacode = oldattcode(randTrials);
        oicode = oldigncode(randTrials);
        avoice = attendedvoice(randTrials);
        ivoice = ignoredvoice(randTrials);
        
        oavmcounting = find(avoice == 'm');
        oavfcounting = find(avoice == 'f');
        oivmcounting = find(ivoice == 'm');
        oivfcounting = find(ivoice == 'f');
        
        Presentedo=Presented(randTrials);
        
        sentencesnames = cell(ntrials,3);
        sentencestypes = cell(ntrials,3);
        sentencescodes = cell(ntrials,3);
        sentencesfreq = cell(ntrials,3);
        
        aocounting = 1;
        ancounting = 1;
        iocounting = 1;
        incounting = 1;
        japcode = zeros(ntrials,1);
        
        %randomize trials
        %         indicenT = randperm(ntrials*2);
        trialstomerge_across = [trialstomerge_across,nTraining+nlast+[1:2*ntrials]];
        
        %trials for test on attended side
        for nT=1:ntrials/2
            
            oac = oacode(oavmcounting(nT));
            idoa = find(Codes == oac);
            index_acrosspairs = find(acrosslistpairs == acrosslistpairs(idoa));
            idna = setdiff(index_acrosspairs,idoa);
            nac = Codes(idna);
            
            if nT<=ntrials/4
                oorder = aoorder1(aocounting);
                norder = aoorder2(aocounting);
                aocounting = aocounting+1;
                %match the japanese word
                japorder = oorder;
                jpc = oac;
                idj = idoa;
                sentencescodes(nT,1:3) = {[oac],[nac],[oac]};
                goodresp=[1,-1];
            else
                oorder = aoorder1(ancounting);
                norder = aoorder2(ancounting);
                %             ancounting = ancounting+1;
                japorder = norder;
                jpc = nac;
                idj = idna;
                sentencescodes(nT,1:3) = {[oac],[nac],[nac]};
                goodresp=[-1,1];
            end
            sentencesnames(nT,1:3) = {NamesF{idna},NamesF{idoa},NamesJ{idj}};
            sentencestypes(nT,1:3) = {Categories(idoa),Categories(idna),Categories(idj)};
            sentencesfreq(nT,1:3) = {Freq(oac),Freq(oicode(nT)),Freq(jpc)};
            [dummy,first] = min([oorder,norder]);second=setdiff(1:2,first);
            
            %first word
            oldlist = struct('TypeCode',sentencestypes{nT,1},'Code',sentencescodes{nT,1},'Freq',sentencesfreq(nT,1),'Name',sentencesnames{nT,1},'Order',oorder,'Side','Stereo','Voice','m','Presented',0);
            %second word
            newlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'Freq',sentencesfreq(nT,2),'Name',sentencesnames{nT,2},'Order',norder,'Side','Stereo','Voice','m','Presented',0);
            %japanese word
            japlist = struct('TypeCode',sentencestypes{nT,3},'Code',sentencescodes{nT,3},'Freq',sentencesfreq(nT,3),'Name',sentencesnames{nT,3},'Order',japorder,'Side','Stereo','Voice','m','Presented',0);
            
            TrialStruct = struct('TrialType',['Test_',title_list, '_across'],'Condition','attended','Number',nT,'FileName',sprintf('%s%st%so%s%s%s.mat','a','m',num2str(sentencescodes{nT,3}),num2str(sentencescodes{nT,first}),num2str(sentencescodes{nT,second}),tech),'old',oldlist,'new',newlist,'jap',japlist);
            
            
            nTtot = nlast+nTraining+nT;
            SessionStructureTest{nTtot} = TrialStruct;
            SummaryTest{nTtot,1} = nTtot;
            SummaryTest(nTtot,2:12) = {['Test_',title_list, '_across'],oac,'oa',avoice(nT),oorder,nac,'na','attended',norder,jpc,Presentedo(nT)};
            
            %attended f voice
            oac = oacode(oavfcounting(nT));
            idoa = find(Codes == oac);
            index_acrosspairs = find(acrosslistpairs == acrosslistpairs(idoa));
            idna = setdiff(index_acrosspairs,oac);
            nac = Codes(idna);
            
            if nT<=ntrials/4
                oorder = aoorder1(aocounting);
                norder = aoorder2(aocounting);
                aocounting = aocounting+1;
                %match the japanese word
                japorder = oorder;
                jpc = oac;
                idj = idoa;
                sentencescodes(nT,1:3) = {[oac],[nac],[oac]};
                goodresp=[1,-1];
            else
                oorder = aoorder1(ancounting);
                norder = aoorder2(ancounting);
                ancounting = ancounting+1;
                japorder = norder;
                jpc = nac;
                idj = idna;
                sentencescodes(nT,1:3) = {[oac],[nac],[nac]};
                goodresp=[-1,1];
            end
            sentencesnames(nT,1:3) = {NamesF{idna},NamesF{idoa},NamesJ{idj}};
            sentencestypes(nT,1:3) = {Categories(idoa),Categories(idna),Categories(idj)};
            sentencesfreq(nT,1:3) = {Freq(oac),Freq(nac),Freq(jpc)};
            [dummy,first] = min([oorder,norder]);second=setdiff(1:2,first);
            
            %first word
            oldlist = struct('TypeCode',sentencestypes{nT,1},'Code',sentencescodes{nT,1},'Freq',sentencesfreq(nT,1),'Name',sentencesnames{nT,1},'Order',oorder,'Side','Stereo','Voice','f','Presented',0);
            %second word
            newlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'Freq',sentencesfreq(nT,2),'Name',sentencesnames{nT,2},'Order',norder,'Side','Stereo','Voice','f','Presented',0);
            %japanese word
            japlist = struct('TypeCode',sentencestypes{nT,3},'Code',sentencescodes{nT,3},'Freq',sentencesfreq(nT,3),'Name',sentencesnames{nT,2},'Order',japcode,'Side','Stereo','Voice','f','Presented',0);
            
            TrialStruct = struct('TrialType',['Test_',title_list, '_across'],'Condition','attended','Number',nT,'FileName',sprintf('%s%st%so%s%s%s.mat','a','f',num2str(sentencescodes{nT,3}),num2str(sentencescodes{nT,first}),num2str(sentencescodes{nT,second}),tech),'old',oldlist,'new',newlist,'jap',japlist);
            
            nTtot = nlast+nTraining+ntrials/2+nT;
            SessionStructureTest{nTtot} = TrialStruct;
            SummaryTest{nTtot,1} = nTtot;
            SummaryTest(nTtot,2:12) = {['Test_',title_list, '_across'],oac,'oa','f',oorder,nac,'na','attended',norder,jpc,Presentedo(nT)};
            
            %test on ignored voice
            %ignored male voice
            ioc = oicode(oivmcounting(nT));
            idoi = find(Codes == ioc);
            index_acrosspairs = find(acrosslistpairs == acrosslistpairs(idoi));
            idni = setdiff(index_acrosspairs,idoi);
            nicode = Codes(idni);
           
            if nT<=ntrials/4
                oorder = ioorder1(iocounting);
                norder = ioorder2(iocounting);
                aocounting = iocounting+1;
                japorder = oorder;
                japcode = ioc;
                idj = idoi;
                sentencescodes(nT,1:3) = {[ioc],[nicode],[ioc]};
                goodresp=[1,-1];
            else
                oorder = inorder1(incounting);
                norder = inorder2(incounting);
                ancounting = incounting+1;
                japorder = norder;
                japcode = nicode;
                idj = idni;
                sentencescodes(nT,1:3) = {[ioc],[nicode],[nicode]};
                goodresp=[-1,1];
            end
            
            sentencesnames(nT,1:3) = {NamesF{idoi},NamesF{idni},NamesJ{idj}};
            sentencestypes(nT,1:3) = {Categories(idoi),Categories(idni),Categories(idj)};
            [dummy,first] = min([oorder,norder]);second=setdiff(1:2,first);
            %first word
            oldlist = struct('TypeCode',sentencestypes{nT,1},'Code',sentencescodes{nT,1},'Freq',sentencesfreq(nT,1),'Name',sentencesnames{nT,1},'Order',oorder,'Side','Stereo','Voice','m','Presented',0);
            %second word
            newlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'Freq',sentencesfreq(nT,2),'Name',sentencesnames{nT,2},'Order',norder,'Side','Stereo','Voice','m','Presented',0);
            %japanese word
            japlist = struct('TypeCode',sentencestypes{nT,3},'Code',sentencescodes{nT,3},'Freq',sentencesfreq(nT,3),'Name',sentencesnames{nT,2},'Order',japorder,'Side','Stereo','Voice','m','Presented',0);
            
            TrialStruct = struct('TrialType',['Test_',title_list, '_across'],'Condition','ignored','Number',nT,'FileName',sprintf('%s%st%so%s%s%s.mat','a','m',num2str(sentencescodes{nT,3}),num2str(sentencescodes{nT,first}),num2str(sentencescodes{nT,second}),tech),'old',oldlist,'new',newlist,'jap',japlist);
            
            nTtot = nlast+nTraining+ntrials+nT;
            SessionStructureTest{nTtot} = TrialStruct;
            SummaryTest{nTtot,1} = nTtot;
            SummaryTest(nTtot,2:12) = {['Test_',title_list, '_across'],ioc,'oi','m',oorder,nicode,'ni','ignored',norder,japcode,Presentedo(nT)};
            
            %ignored female voice
            ioc = oicode(oivfcounting(nT));
            idoi = find(Codes == ioc);
            index_acrosspairs = find(acrosslistpairs == acrosslistpairs(idoi));
            idni = setdiff(index_acrosspairs,idoi);
            nicodek = Codes(idni);
            
            if nT<=ntrials/4
                oorder = ioorder1(iocounting);
                norder = ioorder2(iocounting);
                aocounting = iocounting+1;
                japorder = oorder;
                japcode = ioc;
                idj = idoi;
                sentencescodes(nT,1:3) = {[ioc],[nicodek],[ioc]};
                goodresp=[1,-1];
            else
                oorder = inorder1(incounting);
                norder = inorder2(incounting);
                ancounting = incounting+1;
                japorder = norder;
                japcode = nicode;
                idj = idni;
                sentencescodes(nT,1:3) = {[ioc],[nicodek],[nicodek]};
                goodresp=[-1,1];
            end
            
            sentencesnames(nT,1:3) = {NamesF{idoi},NamesF{idni},NamesJ{idj}};
            sentencestypes(nT,1:3) = {Categories(idoi),Categories(idni),Categories(idj)};
            [dummy,first] = min([oorder,norder]);second=setdiff(1:2,first);
            %first word
            oldlist = struct('TypeCode',sentencestypes{nT,1},'Code',sentencescodes{nT,1},'Freq',sentencesfreq(nT,1),'Name',sentencesnames{nT,1},'Order',oorder,'Side','Stereo','Voice','f','Presented',0);
            %second word
            newlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'Freq',sentencesfreq(nT,2),'Name',sentencesnames{nT,2},'Order',norder,'Side','Stereo','Voice','f','Presented',0);
            %japanese word
            japlist = struct('TypeCode',sentencestypes{nT,3},'Code',sentencescodes{nT,3},'Freq',sentencesfreq(nT,3),'Name',sentencesnames{nT,2},'Order',japorder,'Side','Stereo','Voice','f','Presented',0);
            
            TrialStruct = struct('TrialType',['Test_',title_list, '_across'],'Condition','ignored','Number',nT,'FileName',sprintf('%s%st%so%s%s%s.mat','a','f',num2str(sentencescodes{nT,3}),num2str(sentencescodes{nT,first}),num2str(sentencescodes{nT,second}),tech),'old',oldlist,'new',newlist,'jap',japlist);
            
            nTtot = nlast+nTraining+ntrials/2*3+nT;
            SessionStructureTest{nTtot} = TrialStruct;
            SummaryTest{nTtot,1} = nTtot;
            SummaryTest(nTtot,2:12) = {['Test_',title_list, '_across'],ioc,'oi','f',oorder,nicode,'ni','ignored',norder,japcode,Presentedo(nT)};
            
        end
        
        %%
        % build sentences list and structure for each trial
        %within list
        
        ncond = numel(order);
        ncounterb = floor(ntrials/8/ncond)*ncond;%/2 because it will be done for old primed trials and new primed trials
        %and for the same attended voice
        nrest = ntrials-ncounterb;
        
        for nT=1:numel(order)
            iinv = mod(nT,2)+1;
            suborder = ncounterb/numel(order);
            attendedoldorder((nT-1)*suborder+1:nT*suborder,1) = repmat(order(nT),suborder,1);
            attendedneworder((nT-1)*suborder+1:nT*suborder,1) = repmat(order(iinv),suborder,1);
        end
        
        %%
        % Add trials which are not counterbalanced
        
        for nT=1:nrest
            randorder = randi(1:numel(order));
            ordinv = mod(randorder,2)+1;
            attendedoldorder(ncounterb+nT,1) = order(randorder);
            attendedneworder(ncounterb+nT,1) = order(ordinv);
        end
        
        %%
        %attended old
        
        
        permcond = randperm(ntrials);
        aovforder1 = aoorder(permcond);%ao : attended old
        aovforder2 = anorder(permcond);%an : attended new
        
        permcond = randperm(ntrials);
        aovmorder1 = aoorder(permcond);%ao : attended old
        aovmorder2 = anorder(permcond);%an : attended new
        
        permcond = randperm(ntrials);
        iovforder1 = aoorder(permcond);
        iovforder2 = anorder(permcond);
        
        %ignored old
        permcond = randperm(ntrials);
        iovmorder1 = aoorder(permcond);
        iovmorder2 = anorder(permcond);
        
        %     %building trials
        %     randTrialsa = randperm(length(testTrials));
        % %     randTrialsi = randperm(length(testTrials));
        % %
        %
        sentencesnames = cell(ntrials,3);
        sentencestypes = cell(ntrials,3);
        sentencescodes = cell(ntrials,3);
        sentencesfreq = cell(ntrials,3);
        
        
        nlast = nlast+2*ntrials;
        trialstomerge_within = [trialstomerge_within,nTraining+nlast+[1:ntrials]];
        %%
        for nT=1:ntrials/4
            
            
            %feminine attended
            oac = oacode(oavfcounting(nT));
            oa2c = oacode(oavfcounting(ntrials/4+nT));
            idoa = find(Codes == oac);
            idoa2 = find(Codes == oa2c);
            
            oorder = aovforder1(nT);
            o2order = aovforder2(nT);
            japorder = o2order;
            japcode = oa2c;
            idj = idoa2;
            
            %match the japanese word
            sentencesnames(nT,1:3) = {NamesF{idoa2},NamesF{idoa},NamesJ{idj}};
            sentencestypes(nT,1:3) = {Categories(idoa2),Categories(idoa),Categories(idj)};
            [dummy,first] = min([oorder,norder]);second=setdiff(1:2,first);
            sentencesfreq(nT,1:3) = {Freq(oa2c),Freq(oac),Freq(japcode)};
            sentencescodes(nT,1:3) = {[oa2c],[oac],[japcode]};
            goodresp=[1,-1];
            
            %first word
            oldlist = struct('TypeCode',sentencestypes{nT,1},'Code',sentencescodes{nT,1},'Freq',sentencesfreq(nT,1),'Name',sentencesnames{nT,1},'Order',o2order,'Side','Stereo','Voice','f','Presented',0);
            %second word
            newlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'Freq',sentencesfreq(nT,2),'Name',sentencesnames{nT,2},'Order',oorder,'Side','Stereo','Voice','f','Presented',0);
            %japanese word
            japlist = struct('TypeCode',sentencestypes{nT,3},'Code',sentencescodes{nT,3},'Freq',sentencesfreq(nT,3),'Name',sentencesnames{nT,3},'Order',japorder,'Side','Stereo','Voice','f','Presented',0);
            
            TrialStruct = struct('TrialType',['Test_',title_list, '_within'],'Condition','attended','Number',nT,'FileName',sprintf('%s%st%so%s%s%s.mat','w','f',num2str(sentencescodes{nT,3}),num2str(sentencescodes{nT,first}),num2str(sentencescodes{nT,second}),tech),'old',oldlist,'new',newlist,'jap',japlist);
            
            nTtot = nlast+nTraining+nT;
            SessionStructureTest{nTtot} = TrialStruct;
            SummaryTest{nTtot,1} = nTtot;
            SummaryTest(nTtot,2:12) = {['Test_',title_list, '_within'],oa2c,'oac','f',o2order,oac,'oa2c','attended',oorder,oa2c,Presentedo(nT)};
            
            
            %masculine attended
            
            oac = oacode(oavmcounting(nT));
            oa2c = oacode(oavmcounting(ntrials/4+nT));
            idoa = find(Codes == oac);
            idoa2 = find(Codes == oa2c);
            
            oorder = aovmorder1(nT);
            o2order = aovmorder2(nT);
            %match the japanese word
            
            japorder = o2order;
            japcode = oa2c;
            idj = idoa2;
            sentencescodes(nT,1:3) = {[oa2c],[oac],[oa2c]};
            goodresp=[1,-1];
            
            sentencesnames(nT,1:3) = {NamesF{idoa2},NamesF{idoa},NamesJ{idj}};
            sentencestypes(nT,1:3) = {Categories(idoa2),Categories(idoa),Categories(idj)};
            [dummy,first] = min([o2order,oorder]);second=setdiff(1:2,first);
            sentencesfreq(nT,1:3) = {Freq(oa2c),Freq(oac),Freq(japcode)};
            
            %first word
            newlist = struct('TypeCode',sentencestypes{nT,1},'Code',sentencescodes{nT,1},'Freq',sentencesfreq(nT,1),'Name',sentencesnames{nT,1},'Order',o2order,'Side','Stereo','Voice','m','Presented',0);
            %second word
            oldlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'Freq',sentencesfreq(nT,2),'Name',sentencesnames{nT,2},'Order',oorder,'Side','Stereo','Voice','m','Presented',0);
            %japanese word
            japlist = struct('TypeCode',sentencestypes{nT,3},'Code',sentencescodes{nT,3},'Freq',sentencesfreq(nT,3),'Name',sentencesnames{nT,3},'Order',japorder,'Side','Stereo','Voice','m','Presented',0);
            
            TrialStruct = struct('TrialType',['Test_',title_list, '_within'],'Condition','attended','Number',nT,'FileName',sprintf('%s%st%so%s%s%s.mat','w','m',num2str(sentencescodes{nT,3}),num2str(sentencescodes{nT,first}),num2str(sentencescodes{nT,second}),tech),'old',oldlist,'new',newlist,'jap',japlist);
            
            nTtot = nlast+nTraining+ntrials/4+nT;
            SessionStructureTest{nTtot} = TrialStruct;
            SummaryTest{nTtot,1} = nTtot;
            SummaryTest(nTtot,2:12) = {['Test_',title_list, '_within'],oa2c,'oac','m',o2order,oac,'oa2c','attended',oorder,oa2c,Presentedo(nT)};
            
            
            %ignored feminineS
            oic = oicode(oivfcounting(nT));
            oi2c = oicode(oivfcounting(ntrials/4+nT));
            idoi = find(Codes == oic);
            idoi2 = find(Codes == oi2c);
            
            oorder = iovforder1(nT);
            o2order = iovforder2(nT);
            %match the japanese word
            japorder = o2order;
            
            japcode = oi2c;
            idj = idoi2;
            sentencescodes(nT,1:3) = {[oi2c],[oic],[oi2c]};
            goodresp=[1,-1];
            
            sentencesnames(nT,1:3) = {NamesF{idoa2},NamesF{idoa},NamesJ{idj}};
            sentencestypes(nT,1:3) = {Categories(idoa2),Categories(idoa),Categories(idj)};
            [dummy,first] = min([o2order,oorder]);second=setdiff(1:2,first);
            sentencesfreq(nT,1:3) = {Freq(oi2c),Freq(oic),Freq(japcode)};
            
            %first word
            newlist = struct('TypeCode',sentencestypes{nT,1},'Code',sentencescodes{nT,1},'Freq',sentencesfreq(nT,1),'Name',sentencesnames{nT,1},'Order',o2order,'Side','Stereo','Voice','f','Presented',0);
            %second word
            oldlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'Freq',sentencesfreq(nT,2),'Name',sentencesnames{nT,2},'Order',oorder,'Side','Stereo','Voice','f','Presented',0);
            %japanese word
            japlist = struct('TypeCode',sentencestypes{nT,3},'Code',sentencescodes{nT,3},'Freq',sentencesfreq(nT,3),'Name',sentencesnames{nT,3},'Order',japorder,'Side','Stereo','Voice','f','Presented',0);
            
            TrialStruct = struct('TrialType',['Test_',title_list, '_within'],'Condition','ignored','Number',nT,'FileName',sprintf('%s%st%so%s%s%s.mat','w','f',num2str(sentencescodes{nT,3}),num2str(sentencescodes{nT,first}),num2str(sentencescodes{nT,second}),tech),'old',oldlist,'new',newlist,'jap',japlist);
            
            nTtot = nlast+nTraining+ntrials/2+nT;
            SessionStructureTest{nTtot} = TrialStruct;
            SummaryTest{nTtot,1} = nTtot;
            SummaryTest(nTtot,2:12) = {['Test_',title_list, '_within'],oi2c,'oic','f',o2order,oic,'oi2c','ignored',oorder,oi2c,Presentedo(nT)};
            
            % ignored masculine
            oic = oicode(oivmcounting(nT));
            oi2c = oicode(oivmcounting(ntrials/4+nT));
            idoi = find(Codes == oic);
            idoi2 = find(Codes == oi2c);
            
            oorder = iovmorder1(nT);
            o2order = iovmorder2(nT);
            %match the japanese word
            japorder = o2order;
            
            idj = idoi2;
            sentencescodes(nT,1:3) = {[oi2c],[oic],[oi2c]};
            [dummy,first] = min([o2order,oorder]);
            second=setdiff(1:2,first);
            
            sentencesnames(nT,1:3) = {NamesF{idoi2},NamesF{idoi},NamesJ{idj}};
            sentencestypes(nT,1:3) = {Categories(idoi2),Categories(idoi),Categories(idj)};
            goodresp=[1,-1];
            
            %first word
            oldlist = struct('TypeCode',sentencestypes{nT,1},'Code',sentencescodes{nT,1},'Freq',sentencesfreq(nT,1),'Name',sentencesnames{nT,1},'Order',o2order,'Side','Stereo','Voice','m','Presented',0);
            %second word
            newlist = struct('TypeCode',sentencestypes{nT,2},'Code',sentencescodes{nT,2},'Freq',sentencesfreq(nT,2),'Name',sentencesnames{nT,2},'Order',oorder,'Side','Stereo','Voice','m','Presented',0);
            %japanese word
            japlist = struct('TypeCode',sentencestypes{nT,3},'Code',sentencescodes{nT,3},'Freq',sentencesfreq(nT,3),'Name',sentencesnames{nT,3},'Order',japorder,'Side','Stereo','Voice','m','Presented',0);
            
            TrialStruct = struct('TrialType',['Test_',title_list, '_within'],'Condition','ignored','Number',nT,'FileName',sprintf('%s%st%so%s%s%s.mat','w','m',num2str(sentencescodes{nT,3}),num2str(sentencescodes{nT,first}),num2str(sentencescodes{nT,second}),tech),'old',oldlist,'new',newlist,'jap',japlist);
            
            nTtot = nlast+nTraining+3*ntrials/4+nT;
            SessionStructureTest{nTtot} = TrialStruct;
            SummaryTest{nTtot,1} = nTtot;
            SummaryTest(nTtot,2:12) = {['Test_',title_list, '_within'],oi2c,'oic','m',o2order,oic,'oi2c','ignored',oorder,oi2c,Presentedo(nT)};
            
        end
        
        nlast = nlast+ntrials;
     end
       
        %%
        %merge and shuffle trials independantly for within and across test
        
        %randperm the trials
        rpacross = randperm(length(trialstomerge_across));
        rpwithin = randperm(length(trialstomerge_within));
        
        %store trials
        SummaryTestacross = SummaryTest(trialstomerge_across(rpacross),:);
        SessionStructureTestacross = SessionStructureTest(trialstomerge_across(rpacross));
        SummaryTestwithin = SummaryTest(trialstomerge_within(rpwithin),:);
        SessionStructureTestwithin = SessionStructureTest(trialstomerge_within(rpwithin));
        
        %across
        SummaryTest(nTraining+nfirst+[1:length(rpacross)],:) = SummaryTestacross;
        SessionStructureTest(nTraining+nfirst+[1:length(rpacross)],1) = SessionStructureTestacross;
        for trialsnr=1:length(rpwithin)
            SummaryTest{trialsnr,1} = nTotWake+nTotSleep+nfirst+trialsnr;
            SessionStructureTest{trialsnr}.number = nTotWake+nTotSleep+nfirst+trialsnr;
        end
        
        title_test =  cellfun(@cell2mat,listtests(tnr),'uni',0);
        Memory_Test_Transition((tnr-1)*2+1).trial = nTotWake+nTotSleep+nfirst+1;
        Memory_Test_Transition((tnr-1)*2+1).name = ['Test_',title_test{:}, '_across'];
        
        %within
        SummaryTest(nTraining+nfirst+length(rpacross)+[1:length(rpwithin)],:)=SummaryTestwithin;
        SessionStructureTest(nTraining+nfirst+length(rpacross)+[1:length(rpwithin)],1) = SessionStructureTestwithin;
        for trialsnr = 1:length(rpwithin)
            SummaryTest{trialsnr,1} = nTotWake+nTotSleep+length(rpacross)+nfirst+trialsnr;
            SessionStructureTest{trialsnr}.number = nTotWake+nTotSleep+length(rpacross)+nfirst+trialsnr;
        end
        
        Memory_Test_Transition((tnr-1)*2+2).trial = nTotWake+nTotSleep+nfirst+length(rpacross)+1;
        Memory_Test_Transition((tnr-1)*2+2).name = ['Test_',title_test{:}, '_within'];
    
end


filename = [stimPath, '\\SubjectSessions\\Subject', int2str(SubjectID)];
mkdir(filename);
save([filename, '\\SessionTest'],'SessionStructureTest')
save([filename, '\\SummaryTest.mat'],'SummaryTest')

nlast = length(SummaryTest);
