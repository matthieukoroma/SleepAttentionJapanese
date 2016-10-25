%merge SessionStructure and SessionStructureTest

%% put behaveTraining Test trials in the global structure
SessionStructureLearning(nEEGTraining+2:2:nEEGTraining+nbehavTraining*2) = SessionStructureTest(1:nbehavTraining);
SummaryLearning(nEEGTraining+2:2:nEEGTraining+nbehavTraining*2,:) = SummaryTest(1:nbehavTraining,:);

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
ind=find(test_items_position==stems+1);
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
for i=1:nbehavTraining*2
    if bitget(i,1)
        SessionStructureLearning{i}.TrialType = 'BehavTrainingItems';
        SummaryLearning{i,2} = 'BehavTrainingItems';
    else
        SessionStructureLearning{i}.TrialType = 'BehavTrainingTest';
        SummaryLearning{i,2} = 'BehavTrainingTest';
    end
end

SessionStructureLearning(nTr+vectorItemsind) = SSwakeitems;
SummaryLearning(nTr+vectorItemsind,:) = SLwakeitems;
for i=1:length(vectorItemsind);
    SessionStructureLearning{nTr+vectorItemsind(i)}.TrialType = 'WakeTestItems';
    SummaryLearning{nTr+vectorItemsind(i),2} = 'WakeTestItems';
end

SessionStructureLearning(nTr+vectorTestind) = SSwaketest;
SummaryLearning(nTr+vectorTestind,:) = SLwaketest;
for i=1:length(vectorTestind);
    SessionStructureLearning{nTr+vectorTestind(i)}.TrialType = 'WakeTest';
    SummaryLearning{nTr+vectorTestind(i),2} = 'WakeTestTrials';
end

SessionStructureLearning(nTr+vectorwakeind) = SSwake;
SummaryLearning(nTr+vectorwakeind,:) = SLwake;
for i=1:length(vectorwakeind);
    SessionStructureLearning{nTr+vectorwakeind(i)}.TrialType = 'Wake';
    SummaryLearning{nTr+vectorwakeind(i),2} = 'Wake';
end
%% concatenate both structures
SessionStructureExp = [SessionStructureLearning;SessionStructureTest(nbehavTraining+nwakeTest+1:end)];
SummaryExp = [SummaryLearning;SummaryTest(nbehavTraining+nwakeTest+1:end,:)];

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