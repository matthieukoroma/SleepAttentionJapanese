%merge SessionStructure and SessionStructureTest

%% put behaveTraining Test trials in the global structure
SessionStructureLearning_final(nEEGTraining+2:2:nEEGTraining+nbehavTraining*2) = SessionStructureTest(1:nbehavTraining);
SummaryLearning_final(nEEGTraining+2:2:nEEGTraining+nbehavTraining*2,:) = SummaryTest(1:nbehavTraining,:);

%% Put learning trial tests

%take the nwakeTest trial of BehavTraining and inject them at random position into wake
stems = floor(nwake/nwakeTest);
test_position = stems:stems:nwake;
test_items_position = randi(stems,[1,nwakeTest]);

% make ordering vector
vectorwake = [1:nwake;repmat(2,1,nwake)];
if nwake>0;
    vectorwakeTestitems = [[[0,test_position(1:end-1)]+test_items_position];repmat(1,1,length(test_items_position))];
    %arrange stems at max
    ind=find(test_items_position==stems);
    tp1=test_position(1:end);
    vectorwakeTestitems(:,ind)=[tp1(ind)-1;repmat(3,1,length(ind))];
    
    vectorwakeTestTrials = [[test_position];repmat(4,1,length(test_items_position))];
end
vectorTest = [vectorwake,vectorwakeTestitems,vectorwakeTestTrials];
sortvectorTest=sortrows(vectorTest',1);
%rearrange vectors
vectorItemsind = find(ismember(sortvectorTest(:,2),[1,3]));
vectorTestind = find(sortvectorTest(:,2)==4);
vectorwakeind = find(sortvectorTest(:,2)==2);

%%
%Take the items for SessionLearning
nTr = nEEGTraining+2*nbehavTraining;
SSwakeitems = SessionStructureLearning(nTr+[1:2:2*nwakeTest]);
SLwakeitems = SummaryLearning(nTr+[1:2:2*nwakeTest],:);
SSwake = SessionStructureLearning(nTr+2*nwakeTest+[1:nwake],:);
SLwake = SummaryLearning(nTr+2*nwakeTest+[1:nwake],:);

SSwaketest = SessionStructureTest(nbehavTraining+[1:nwakeTest]);
SLwaketest = SummaryTest(nbehavTraining+[1:nwakeTest],:);


%% replace
SessionStructureLearning_final=SessionStructureLearning;
SummaryLearning_final=SummaryLearning;
for i=1:nbehavTraining*2
    if bitget(i,1)
        SessionStructureLearning_final{i}.TrialType = 'BehavTrainingItems';
        SummaryLearning_final{i,2} = 'BehavTrainingItems';
    else
        SessionStructureLearning_final{i}.TrialType = 'BehavTrainingTest';
        SummaryLearning_final{i,2} = 'BehavTrainingTest';
    end
end

SessionStructureLearning_final(nTr+vectorItemsind) = SSwakeitems;
SummaryLearning_final(nTr+vectorItemsind,:) = SLwakeitems;
for i=1:length(vectorItemsind);
    SessionStructureLearning_final{nTr+vectorItemsind(i)}.TrialType = 'WakeTestItems';
    SummaryLearning_final{nTr+vectorItemsind(i),2} = 'WakeTestItems';
end

SessionStructureLearning_final(nTr+vectorTestind) = SSwaketest;
SummaryLearning_final(nTr+vectorTestind,:) = SLwaketest;
for i=1:length(vectorTestind);
    SessionStructureLearning_final{nTr+vectorTestind(i)}.TrialType = 'WakeTest';
    SummaryLearning_final{nTr+vectorTestind(i),2} = 'WakeTestTrials';
end

SessionStructureLearning_final(nTr+vectorwakeind) = SSwake;
SummaryLearning_final(nTr+vectorwakeind,:) = SLwake;
for i=1:length(vectorwakeind);
    SessionStructureLearning_final{nTr+vectorwakeind(i)}.TrialType = 'Wake';
    SummaryLearning_final{nTr+vectorwakeind(i),2} = 'Wake';
end
%% concatenate both structures
SessionStructureExp = [SessionStructureLearning_final;SessionStructureTest(nbehavTraining+nwakeTest+1:end)];
SummaryExp = [SummaryLearning_final;SummaryTest(nbehavTraining+nwakeTest+1:end,:)];

%%
listN= {'SleepN2','SleepN3','falling'};
for i=1:length(listN)
    if strcmp(listN{i},'SleepN2') && nN2>0
        stim_file = cellfun(@(x) strcmp(x,'SleepN2'),TrialsTypeSession);
        index1 = find(stim_file);
        listnames.SleepN2.list = index1;
        listnames.SleepN2.active = index1(1);
    elseif strcmp(listN{i},'SleepN3') && nN3>0
        stim_file = cellfun(@(x) strcmp(x,'SleepN3'),TrialsTypeSession);
        index1 = find(stim_file);
        listnames.SleepN3.list = index1;
        listnames.SleepN3.active = index1(1);
    elseif strcmp(listN{i},'falling') && nfalling>0
        stim_file = cellfun(@(x) strcmp(x,'falling'),TrialsTypeSession);
        index1 = find(stim_file);
        listnames.falling.list = index1;
        listnames.falling.active = index1(1);
    end
end

%%
TestTrialsAwake = [nEEGTraining+[2:2:2*nbehavTraining],nTr+vectorTestind(:)'];