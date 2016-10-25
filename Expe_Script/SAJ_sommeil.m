%%%% Sleep Attention Japanese; Sleep version 1.0
%%%% Original: Sleep Attention Allocation
%%%% Thomas Andrillon - 21/01/2011
%%%% New Thomas Andrillon - 10/03/2014
%%%% Modified version Guillaume Legendre - 06/02/2015
%%%% Modified for Sleep Attention Japanese - 28/03/2015


%%
% mycomp = 1;
% if mycomp
%     AcquisitionFlag=0; % for test on computer
% end

%%



%%
%%%%%%                  Matlab Initialization
clear all
close all

mycomp=0;


cd('C:\Data\SleepAttentionJapanese\Expe_Folder\')
rand('state',sum(100*clock));
warning('off','MATLAB:dispatcher:InexactMatch'); % Swith off the case-sensitive warnings
source_path = genpath(fullfile(pwd,'../'));
% addpath(source_path)
SyncLag=[];
logstring = [];
maxdiffsynchro=30; % increase if synchronisation fail or too long
maxloopsynchro=500;
%dataPath = 'D:\\LSCPData\\SleepAttentionAllocation\\Expe_Folder\\Expe_Data'
%stimPath = 'D:\\LSCPData\\SleepAttentionAllocation\\Expe_Folder\\Expe_Material'
% dataPath=[pwd filesep '..' filesep 'Expe_Data'];
% stimPath=[pwd filesep '..' filesep 'Expe_Material'];
dataPath=[pwd filesep 'Expe_Data'];
stimPath=[pwd filesep 'Expe_Material'];

expeName = 'SAJ';

Category = {'Wo','NW'};
HandSide = {'RIGHT' , 'LEFT'};

if mycomp
    nChannels = 2; %on computer
else 
    nChannels=4; %Number of channels
end
nDevice=30; % AudioFire12:DeviceID=30
runMode=0;

nbtestkeyh = msgbox(sprintf('Did you take care of rebooting Matlab between this subject and the previous one ?\n (Matlab can run out of memory between two SessionStructures. Can cause the stop of the script in the middle of the experiment)'),'Warning : Rebooting','warn');
uiwait(nbtestkeyh);

%%
%%%%%%%%%%              Params manip

% nwake = 16;
% nbehavTraining = 10;
% nwakeTest = 4;
% nN2 = 0;%16;
% nN3 = 0;%16;
% nEEGTraining = 0;

load('C:\Data\SleepAttentionJapanese\Expe_Folder\Expe_Material\TTLson');
%%
%%%%%%%%                Set Home or Hotel Dieu
mycomp = 0;
if mycomp
    AcquisitionFlag=0; % for test on computer
end
AcquisitionFlag=1;
%%
%%%%%%                  Initialization PsychToolBox
if mycomp==0
    try
        AssertOpenGL;
        GetSecs; WaitSecs(ceil(rand(1)*1000)/1000); % time initialisation
        KbName('KeyNamesWindows');
        InitializePsychSound(1);
    catch
        PsychPortAudio('Close');
        %     psychrethrow(psychlasterror);
        AssertOpenGL;
        GetSecs; WaitSecs(ceil(rand(1)*1000)/1000); % time initialisation
        keyTouch=KbName('KeyNamesWindows');
        InitializePsychSound(1);
    end
end
fs=44100;

%%
%%%%%%                  Check keybord

name_key = {'FAR LEFT','MID LEFT','MID RIGHT','FAR RIGHT'};
nbtestkey = 2;
keyUsed = 'dfjk';
AnswerNames = ['S','F','H','K'];
AnswerKeys = [68,70,74,75];
% AnswerKeys = [40:2:46];
% AnswerNumbers = [1:4];
KbName('UnifyKeyNames');
%
% for keynr=1:length(name_key)
%     for testkey=1:nbtestkey
%         nbtestkeyh= msgbox(['Check keyboard is working... Ask subject for one ', name_key{keynr},' press'],'Info','none');
%         pause(0.1);
%         checkStart=GetSecs;
%         [KBKey dummyPressed keyboardKey]= get_SubjectResponse_overPeriod_KB(checkStart,10,1);
%         KBKey
%             while KBKey~=AnswerKeys(keynr)
%             close(nbtestkeyh);
%             nbtestkeyh = msgbox(['Wrong response... Ask subject again for one ', name_key{keynr},' press'],'Info','none');
%             pause(0.1);
%             checkStart=GetSecs;
%             [KBKey dummyPressed keyboardKey] = get_SubjectResponse_overPeriod_KB(checkStart,10,1);
%             KBKey
%         end
%         close(nbtestkeyh);
%     end
% end

% %%
% %%%%%%         x gdfdgtred        SessionStructure parameters
answer = inputdlg({'Code','Group','Initials','Age','Gender','Acquisition','Volume'},'Subject Info',1,{'000','XX','JD','33','M','0','0.6'});
Subject.ID=(answer{1});
Subject.Group=(answer{2});
Subject.Initials=(answer{3});
Subject.Age=(answer{4});
Subject.Gender=(answer{5});
Subject.AcquisitionFlag=str2num(answer{6}); AcquisitionFlag=Subject.AcquisitionFlag;
Subject.volumeSound=(answer{7});
volumeSound=str2num(Subject.volumeSound);

alreadyExistingSubject = 0;
keepAsking = 0;
newSub = 0;

FormatID='XXX';
% if (exist([dataPath filesep expeName Subject.ID])~=0 && strcmp(Subject.ID,'000')~=1)
%     keepAsking = 1;
% end

% Subject.ID = '129';
% Subject.Group='1';
% Subject.Initials='DD';
% Subject.Age=23;
% Subject.Gender='F';
% Subject.AcquisitionFlag=str2num('0'); AcquisitionFlag=Subject.AcquisitionFlag;
% Subject.volumeSound='0.6';
% volumeSound=str2num(Subject.volumeSound)

%%

startingTrial = 1;

while isempty(Subject.ID) || length(Subject.ID)~=length(FormatID) || keepAsking
    if keepAsking==1
        nbtestkeyh = msgbox(sprintf('Subject ID %s already exist',Subject.ID),'Already existing Subject ID','error');
        uiwait(nbtestkeyh);
        
        subjectPanel = figure();
        subjP1 = uicontrol('String','Erase Subject file','Position', [20 240 100 40],'Callback',sprintf('keepAsking=0;\nuiresume(gcbf);\n'));
        subjP2 = uicontrol('String',sprintf('Use the already existing Subject file'),'Position', [140 240 100 40],'Callback',sprintf('keepAsking=0;\nalreadyExistingSubject = 1;\nuiresume(gcbf);\n'));
        subjP3 = uicontrol('String',sprintf('Enter a new Subject ID'),'Position', [260 240 100 40],'Callback',sprintf('keepAsking=0;\nnewSub = 1;\nuiresume(gcbf);\n'));
        subjP4 = uicontrol('Style','text','String',sprintf('What do you want to do ?'),'Position', [20 300 200 40]);
        uiwait(subjectPanel);
        close(subjectPanel)
        
        if alreadyExistingSubject == 1
            sameInfos = 1;
            previousData = load([dataPath filesep expeName Subject.ID filesep 'Result_' expeName  Subject.ID]);
            corresSubj = [  strcmp(previousData.Subject.ID,Subject.ID),strcmp(previousData.Subject.Group,Subject.Group),...
                strcmp(previousData.Subject.Initials,Subject.Initials),strcmp(previousData.Subject.Age,Subject.Age),...
                strcmp(previousData.Subject.Gender,Subject.Gender),previousData.Subject.volumeSound==Subject.volumeSound,...
                previousData.Subject.AcquisitionFlag==Subject.AcquisitionFlag];
            if sum(corresSubj)~=numel(corresSubj)
                incorrectPanel = figure();
                iPanel1 = uicontrol('Style','text','String',sprintf('Data from the previous SessionStructure do not correspond to the new subject data provided.\nPlease, check the subject data.'),'Position', [20 300 200 40]);
                iPanel2 = uitable(incorrectPanel,'Data',{   'ID',Subject.ID,previousData.Subject.ID;...
                    'Group',Subject.Group,previousData.Subject.Group;...
                    'Initials',Subject.Initials,previousData.Subject.Initials;...
                    'Age',Subject.Age,previousData.Subject.Age;...
                    'Gender',Subject.Gender,previousData.Subject.Gender;...
                    'volumeSound',Subject.volumeSound,previousData.Subject.volumeSound;...
                    'AcquisitionFlag',Subject.AcquisitionFlag,previousData.Subject.AcquisitionFlag});
                iPanel3 = uicontrol('String','Continue','Position', [300 20 100 40],'Callback',sprintf('uiresume(gcbf);\n'));
                sizeiTable = get(iPanel2,'Extent');
                set(iPanel2,'Position',sizeiTable);
                uiwait(incorrectPanel);
                close(incorrectPanel);
                sameInfos = 0;
                alreadyExistingSubject = 0;
            end
            if sameInfos
                stoppedEvent = numel(previousData.Timings.stopTime);
                stoppedTrial = previousData.TrialsCaracs(end).Trials;
                inputTrial = inputdlg({sprintf('Recordings are present until trial %g.\nAt which trial do you want to restart the experiment ?',stoppedTrial)},'Starting Trial',1,{int2str(stoppedTrial)});
                startingTrial = str2num(inputTrial{1});
            else
                newSub = 1;
            end
        end
    end
    if isempty(Subject.ID) || length(Subject.ID)~=length(FormatID) || newSub
        nbtestkeyh = msgbox(sprintf('Please enter a new Subject ID (%s)',FormatID),'INVALID SUBJECT ID','error');
        uiwait(nbtestkeyh);
        answer = inputdlg({'Code','Group','Initials','Age','Gender','Acquisition'},'Subject Info',1,{'000','XX','JD','33','M','0'});
        Subject.ID=(answer{1});
        Subject.Group=(answer{2});
        Subject.Initials=(answer{3});
        Subject.Age=(answer{4});
        Subject.Gender=(answer{5});
        Subject.volumeSound=(answer{7});
        Subject.AcquisitionFlag=str2num(answer{6}); AcquisitionFlag=Subject.AcquisitionFlag;
    end
    if (exist([dataPath filesep expeName Subject.ID])~=0 && strcmp(Subject.ID,'000')~=1 && newSub)
        keepAsking = 1;
        newSub = 0;
    end
end



%%
%%%%%%                  Prepare NetSation Sync
if AcquisitionFlag==1
    fprintf('#################### \n #################### \n #################### \n Prepare NetSation \n #################### \n #################### \n #################### \n',1);
    % Initializing communications with the NETSTATION
    [status, error] = NetStationPC_Debugged('connect','10.0.0.42');
    if (status ~= 0)
        nbtestkeyh = msgbox('COMMUNICATION WITH THE NETSTATION IMPOSSIBLE','Error','error');
        uiwait(nbtestkeyh);
    else
        nbtestkeyh = msgbox('NetStation is alive. Press ok to start the recordings','Info','none');
        uiwait(nbtestkeyh);
    end
    
    % Synchronisation between PC and MAC
    fprintf('\n Synchronisation');
    d=GetSecs*1000; e=d+100;
    tmp=0;
    while (e-d>maxdiffsynchro) && (tmp<maxloopsynchro)
        %         [status, error] = NetStationPC_Debugged('attention');
        d=GetSecs*1000;
        [status, error] = NetStationPC_Debugged('synchronize',d);
        e=GetSecs*1000;
        tmp=tmp+1;
        diffsynchro=e-d;
    end
    if tmp==maxloopsynchro
        nbtestkeyh = msgbox('Synchronisation FAILED. Press ok to abort','Error','error');
        uiwait(nbtestkeyh);
        [status, error] = NetStationPC_Debugged('stoprecording');
        %         [status, error] = NetStationPC_Debugged('disconnect');
        return
    end
    logstring = sprintf('%s\nDifference in synchronisation: %g at begining\n',logstring,diffsynchro);
    
    % Start recording on the NETSTATION
    [status, error] = NetStationPC_Debugged('startrecording');
    recstart=GetSecs*1000;
    logstring = sprintf('%s\nStart of Recording: %g\n',logstring,recstart);
    nbtestkeyh = msgbox('Check NetStation is recording and press ok. Program will pause for 10 seconds.','Info','none');
    uiwait(nbtestkeyh);
    WaitSecs(10); % TOUJOURS LAISSER AU MOINS 3000ms APRES UN START RECORD (Augmant\E9 \E0 10sec: Guillaume 05022015)
    % $$Event marking beginning of SessionStructure and info
    [status, error] = NetStationPC_Debugged('event','SESS',GetSecs,0.001,'Expe',1,'Subj',Subject.ID);
    
    nbtestkeyh = msgbox('NetStation is ready and recording. Press ok to resume.','Info','none');
    uiwait(nbtestkeyh);
end


%%
%%%%%%                  Initialization of PsychPortAudio

if mycomp==0
    
    pahandle = PsychPortAudio('Open',nDevice, 1, 2, fs, nChannels);
    PsychPortAudio('Volume',pahandle,volumeSound);
    PsychPortAudio('UseSchedule', pahandle, 1);
    PsychPortAudio('RunMode', pahandle, 0);
    [beepTrigger    ,sr] = MakeBeep(2000, 0.08, fs);
    doubleBeep = [MakeBeep(600,0.08,fs),zeros(1,44100*0.5),MakeBeep(600,0.08,fs)];
    endTrialBeep = [doubleBeep; doubleBeep; zeros(1,numel(doubleBeep)); zeros(1,numel(doubleBeep))];
    beepKey = [MakeBeep(800,0.08,fs);MakeBeep(800,0.08,fs);zeros(1,0.08*44100+1);zeros(1,0.08*44100+1)];
    beepSpace = [MakeBeep(1000,0.08,fs);MakeBeep(1000,0.08,fs);zeros(1,0.08*44100+1);zeros(1,0.08*44100+1)];
    
    for nrms=1:2
        endTrialBeep(nrms,:)=endTrialBeep(nrms,:)/rms(endTrialBeep(nrms,:))*0.07;
        beepKey(nrms,:)=beepKey(nrms,:)/rms(beepKey(nrms,:))*0.07;
        beepSpace(nrms,:)=beepSpace(nrms,:)/rms(beepSpace(nrms,:))*0.07;
    end
    
    endTrialBeep = endTrialBeep(1:nChannels,:);
    beepKey = beepKey(1:nChannels,:);
    beepSpace = beepSpace(1:nChannels,:);
    
end
%%
%%%%%%                  Initialization loop on trials


sampleSize = 10.0;

qY = 200;
rX = 100;
rY = [300,400,500,600];
instrY = 700;

clc;

countNoResponse=0;
LastResponse=GetSecs;
expeStart=GetSecs;

% mem = memory;
running = 1;
exitK = 1;
falling = 0;
awake = 1;
N2 = 0;
N3 = 0;
memorytest = 0;

%intialize variable
nEEGTraining =4;
nwake = 16;
nbehavTraining = 10;%<=12
nfalling = 16;
nN3 = 16;
nN2 = 16;
nwakeTest = 4;
nTestWake = 4;
nTotWake = nEEGTraining+2*nbehavTraining+2*nwakeTest+nwake+nfalling;
nTotSleep = nN2 + nN3;
nNoMemoryTest = nfalling;

fixInterval=4;
jitterInterval=2;
fixIntervalmemorytest = 2;


% WaitSecs(9);

if alreadyExistingSubject
    SessionStructureExp = previousData.SessionStructureExp;
else
    %     SessionStructureExp = createSessionStructureFinal_Sleepv2(str2num(Subject.ID),str2num(Subject.Group),nbehavTraining,nwake,nN2,nN3,nNoMemoryTest);
    gr = str2num(Subject.Group);
    SubjectID = str2num(Subject.ID);
    createSessionStructureFinal_Sleepv2
    createSessionStructureTest
    MergeSessionStructure
    TestTrials = [nEEGTraining+2:2:nbehavTraining,nTr+vectorTestind'];
    
    %arrange SessionStructure
end
nTotTrials = length(SummaryExp);
endTrialJitter = fixInterval+jitterInterval*floor(rand(nTotTrials*10,1)*1000)/1000;

%%
%%%%%%                  Opening of the control panel


% Control Panel
ctrlPanel = figure();
h1 = uicontrol('Position', [20 20 80 40], 'String', 'Continue', 'Callback', sprintf('uiresume(ctrlPanel)\nset(h1,''Enable'',''off'');\ndrawnow;\nset(h1,''Enable'',''on'')\n'));
h3 = uicontrol('String','Skip to next sound','Position', [320 20 100 40],'Callback',sprintf('PsychPortAudio(''Stop'', pahandle, 2);\nset(h3,''Enable'',''off'');\ndrawnow;\nset(h3,''Enable'',''on'')\n'));
h22 = uicontrol('String',sprintf('Replay\nSound'),'Position', [240 20 60 40],'Callback',sprintf('rep = 1;\nPsychPortAudio(''Stop'', pahandle, 2);\nset(h22,''Enable'',''off'');\ndrawnow;\nset(h22,''Enable'',''on'')\n'));
h4 = uicontrol('Style','text','String','What','Position', [180 160 60 20]);
h5 = uicontrol('Style','text','String','Inactive','Position', [80 200 50 20]);
h6 = uicontrol('String','Exit','Position', [440 20 50 40], 'Callback','exitK = 0;');
h7 = uicontrol('Style','text','String','AudioPort:','Position', [20 200 50 20]);
h8 = uicontrol('Style','text','String','Currently playing:','Position', [20 160 100 20]);
h9 = uicontrol('Style','text','String','Trial','Position', [130 160 40 20]);
h10 = uicontrol('Style','text','String','','Position', [130 100 150 20]);
h11 = uicontrol('Style','text','String','','Position', [370 80 80 20]);
h12 = uicontrol('String','Check Memory','Position', [370 100 80 20],'Callback',sprintf('mem=memory;\nset(h11,''String'',[int2str(mem.MemAvailableAllArrays/1000000),'' MB'']);'));
h17 = uicontrol('Style','text','String','Trial Name','Position', [190 130 150 20]);
h18 = uicontrol('Style','text','String','attName','Position', [130 130 50 20]);
h19 = uicontrol('Style','text','String','igName','Position', [70 130 50 20]);
h20 = uicontrol('Style','text','String','attvoice','Position', [310 130 50 20]);
h21 = uicontrol('Style','text','String','attside','Position', [370 130 50 20]);
h25 = uicontrol('Style','text','String', 'Wake Phase', 'Position', [120 80 100 20]);
if nN2>0
    h23 = uicontrol('Position', [120 40 100 20], 'String', 'N2 Phase', 'Callback', sprintf('currtrial=listnames.SleepN2.active+1;Changelist = 1;\nset(h25,''String'',''N2 Phase'');\ndrawnow;\nset(h25,''Enable'',''on'')\n'));
end;if nN3>0;h26 = uicontrol('Position', [120 20 100 20], 'String', 'N3 Phase', 'Callback', sprintf('currtrial=listnames.SleepN3.active+1;Changelist = 1;\nset(h25,''String'',''N3 Phase'');\ndrawnow;\nset(h25,''Enable'',''on'')\n'));
end;if nfalling>0; h24 = uicontrol('Position', [120 60 100 20], 'String', 'Falling asleep Phase', 'Callback', sprintf('currtrial=listnames.falling.active+1;Changelist = 1;\nset(h25,''String'',''Falling asleep Phase'');\ndrawnow;\nset(h25,''Enable'',''on'')\n'));
end;h27 = uicontrol('Position', [120 100 100 20], 'String', 'Memory Test Phase', 'Callback', sprintf('currtrial=nTotWake+nTotSleep+1;Changelist = 1;\nset(h25,''String'',''Memory Test Phase'');\ndrawnow;\nset(h25,''Enable'',''on'')\n'));

uiwait(ctrlPanel);

%%
%%%%%%                  Prepare sounds

if mycomp==0
    expeStart = GetSecs();
    
    handleStatus = PsychPortAudio('GetStatus',pahandle);
    
    keyBuffer = PsychPortAudio('CreateBuffer',pahandle,beepKey);
    spaceBuffer = PsychPortAudio('CreateBuffer',pahandle,beepSpace);
    endTrialBuffer = PsychPortAudio('CreateBuffer',pahandle,endTrialBeep);
    
    startingBeep = zeros(1,2.08*44100+1);
    for i = 1:3
        startingBeep((i-1)*44100+1:(i-1)*44100+0.08*44100+1) = MakeBeep(400,0.08,44100);
    end
    startingBeep = [startingBeep;startingBeep;zeros(1,numel(startingBeep));zeros(1,numel(startingBeep))];
    
    for nrms=1:2
        startingBeep(nrms,:)=startingBeep(nrms,:)/rms(startingBeep(nrms,:))*0.07;
    end
    
    startingBeep = startingBeep(1:nChannels,:);
    startBuffer = PsychPortAudio('CreateBuffer',pahandle,startingBeep);
end

%%

trialType = cell(0);

trialattended = cell(0);
trialaCode = cell(0);
trialaVoice = cell(0);
trialaSide = cell(0);
trialaTag = cell(0);
trialaAmp = cell(0);
trialaPresented = cell(0);

trialignored = cell(0);
trialiOrder = cell(0);
trialiVoice = cell(0);
trialiSide = cell(0);
trialiTag = cell(0);
trialiAmp = cell(0);
trialiPresented = cell(0);

startTimes = [];
stopTimes = [];
stopFunctionTime = [];

playedTrial = cell(0);

if alreadyExistingSubject
    startTimes = cell2mat(previousData.Timings.startTime(1:stoppedEvent));
    stopTimes = cell2mat(previousData.Timings.stopTime(1:stoppedEvent));
    stopFunctionTime = cell2mat(previousData.Timings.stopFunctionTime(1:stoppedEvent,:));
    
    
    for i=1:stoppedEvent
        playedTrial{i,1} = previousData.TrialsCaracs(i).Trials;
        
        trialType{i,1} = previousData.TrialsCaracs(i).TrialType;
        
        trialattended{i,1} = previousData.TrialsCaracs(i).attended;
        trialaCode{i,1} = previousData.TrialsCaracs(i).aCode;
        trialaOrder{i,1} = previousData.TrialsCaracs(i).aOrder;s
        trialaVoice{i,1} = previousData.TrialsCaracs(i).aVoice;
        trialaSide{i,1} = previousData.TrialsCaracs(i).aSide;
        trialaCategories{i,1} = previousData.TrialsCaracs(i).aCategories;
        trialaName{i,1} = previousData.TrialsCaracs(i).aName;
        trialaFreq{i,1} = previousData.TrialsCaracs(i).aFreq;
        %         trialaTag{i,1} = previousData.TrialsCaracs(i).aTag;
        %         trialaAmp{i,1} = previousData.TrialsCaracs(i).aAmp;
        trialaPresented{i,1} = previousData.TrialsCaracs(i).aPresented;
        
        trialignored{i,1} = previousData.TrialsCaracs(i).ignored;
        trialiCode{i,1} = previousData.TrialsCaracs(i).iCode;
        trialiOrder{i,1} = previousData.TrialsCaracs(i).iOrder;
        trialiVoice{i,1} = previousData.TrialsCaracs(i).iVoice;
        trialiSide{i,1} = previousData.TrialsCaracs(i).iSide;
        trialiCategories{i,1} = previousData.TrialsCaracs(i).iCategories;
        trialiName{i,1} = previousData.TrialsCaracs(i).iName;
        trialiFreq{i,1} = previousData.TrialsCaracs(i).iFreq;
        %         trialiTag{i,1} = previousData.TrialsCaracs(i).iTag;
        %         trialiAmp{i,1} = previousData.TrialsCaracs(i).iAmp;
        trialiPresented{i,1} = previousData.TrialsCaracs(i).iPresented;
    end
end


%%
falling=0;
currevent = 1;
currtrial = 1;
awake = 1;
memorytest = 1;
Changelist = 0;
playedTrial{1,1} = 1;
countingTrial = 1;
if alreadyExistingSubject
    currtrial = startingTrial;
    currevent = stoppedEvent+1;
end
% sleepTrials = nEEGTraining+2*nbehavTraining+2*nwakeTest+nwake+[1:nN2+nN3];

%%

% tryf
while running&&exitK
    if mycomp==1
        nChannels=4;
    end
    % Load the first trial, delete the unused buffer and replace it contents by the trial
%     try
        j;
        if mycomp==0
            load([stimPath, '\\testStims\\', SessionStructureExp{currtrial}.FileName]);
        else
            load([stimPath, '/testStims/', SessionStructureExp{currtrial}.FileName]);
        end 
        son = Stim.Mat;
        
        fs = Stim.Caracs.fsample;
        
        if currtrial <= nEEGTraining
            training = 1;
        else
            training = 0;
        end
         Changelist = 0;

        if currtrial > nTotWake-nfalling && currtrial <= nTotWake
            memorytest = 0;falling=1;
            set(h25,'String','Falling asleep phase');
            listnames.falling.active = currtrial;
        elseif currtrial > nTotWake && currtrial <= nTotWake+nN2
            memorytest = 0;
            listnames.SleepN2.active = currtrial;
            set(h25,'String','N2 Phase');
        elseif currtrial >= nTotWake+nN2 && currtrial <= nTotWake+nTotSleep
            memorytest = 0;
            listnames.SleepN3.active = currtrial;
            set(h25,'String','N3 Phase');
        end
        if ismember(currtrial,TestTrialsAwake) || currtrial > nTotWake+nTotSleep
            memorytest = 1;
        else memorytest = 0;
        end
        
        if currtrial==nEEGTraining+1 && nbehavTraining>0
            nbtestkeyh = msgbox(sprintf('Going to the behavorial Training phase ?\n'),'Warning : Memory','warn');
            uiwait(nbtestkeyh);
        elseif currtrial==nEEGTraining+nbehavTraining*2+1
            nbtestkeyh = msgbox(sprintf('Going to the Wake phase ?\n'),'Warning : Memory','warn');
            uiwait(nbtestkeyh);
        elseif currtrial == nTotWake-nfalling+1
            nbtestkeyh = msgbox(sprintf('Going to the falling asleep phase ?\n'),'Warning : Memory','warn');
            uiwait(nbtestkeyh);
        elseif currtrial == nTotWake+nTotSleep+1
            endTrialJitter = fixIntervalmemorytest+jitterInterval*floor(rand(nTotTrials*10,1)*1000)/1000;
            memorytest=1;
%    
%             
%             %verify the keys
%             for keynr=1:length(name_key)
%                 nbtestkeyh= msgbox(['Check keyboard is working... Ask subject for one ', name_key{keynr},' press'],'Info','none');
%                 pause(0.1);
%                 checkStart=GetSecs;
%                 [KBKey dummyPressed keyboardKey]= get_SubjectResponse_overPeriod_KB(checkStart,10,1);
%                 KBKey
%                 while KBKey~=AnswerKeys(keynr)
%                     close(nbtestkeyh);
%                     nbtestkeyh = msgbox(['Wrong response... Ask subject again for one ', name_key{keynr},' press'],'Info','none');
%                     pause(0.1);
%                     checkStart=GetSecs;
%                     [KBKey dummyPressed keyboardKey] = get_SubjectResponse_overPeriod_KB(checkStart,10,1);
%                     KBKey
%                 end
%                 close(nbtestkeyh);
%             end
        end
        if ismember(currtrial,nTotWake+nTotSleep+[Memory_Test_Transition.trial])
            TTname = find([Memory_Test_Transition.trial]==currtrial);
            set(h25,'String','Memory_Test Phase');
            nbtestkeyh = msgbox(sprintf('Going to the Memory Test ?\n'),'Warning : Memory','warn');
            uiwait(nbtestkeyh);      
        end
        if currtrial < nEEGTraining && AcquisitionFlag==1
            
            nextStimCarac = struct( 'FileName',SessionStructureExp{currtrial}.FileName,...
                'TrialType',SessionStructureExp{currtrial}.TrialType,...
                'aCode',SessionStructureExp{currtrial}.attended.Code,'aFreq',SessionStructureExp{currtrial}.attended.Freq,'aOrder',SessionStructureExp{currtrial}.attended.Order,...
                'aName',SessionStructureExp{currtrial}.attended.Name,'aCategories',SessionStructureExp{currtrial}.attended.TypeCode,...
                'aVoice',SessionStructureExp{currtrial}.attended.Voice,'aSide',SessionStructureExp{currtrial}.attended.Side,...
                'aPresented',SessionStructureExp{currtrial}.attended.Presented,'iCode',NaN,'iFreq',NaN,'iType',NaN,'iOrder',NaN,...
                'iName',NaN,'iCategories',NaN,'iVoice',NaN,'iSide',NaN,'jCode',NaN,'iPresented',SessionStructureExp{currtrial}.ignored.Presented);
        end
        

        if memorytest
            nextStimCarac = struct('FileName',SessionStructureExp{currtrial}.FileName,...
                'TrialType',SessionStructureExp{currtrial}.TrialType,...
                'aCode',SessionStructureExp{currtrial}.old.Code,'aFreq',SessionStructureExp{currtrial}.old.Freq,...
                'aOrder',SessionStructureExp{currtrial}.old.Order,...
                'aName',SessionStructureExp{currtrial}.old.Name,'aCategories',SessionStructureExp{currtrial}.old.TypeCode,...
                'aVoice',SessionStructureExp{currtrial}.old.Voice,'aSide',SessionStructureExp{currtrial}.old.Side,...
                'aPresented',SessionStructureExp{currtrial}.old.Presented,...
                'iCode',SessionStructureExp{currtrial}.new.Code,'iFreq',SessionStructureExp{currtrial}.new.Freq,...
                'iOrder',SessionStructureExp{currtrial}.new.Order,'iName',SessionStructureExp{currtrial}.new.Name,...
                'iCategories',SessionStructureExp{currtrial}.new.TypeCode,'iVoice',SessionStructureExp{currtrial}.new.Voice,...
                'iSide',SessionStructureExp{currtrial}.new.Side,'jCode',SessionStructureExp{currtrial}.jap.Code,'iPresented',SessionStructureExp{currtrial}.new.Presented);
        else
            SessionStructureExp{currtrial}.attended.Presented = SessionStructureExp{currtrial}.attended.Presented+1;
            SessionStructureExp{currtrial}.ignored.Presented = SessionStructureExp{currtrial}.attended.Presented+1;
            nextStimCarac = struct('FileName',SessionStructureExp{currtrial}.FileName,...
                'TrialType',SessionStructureExp{currtrial}.TrialType,...
                'aCode',SessionStructureExp{currtrial}.attended.Code,'aFreq',SessionStructureExp{currtrial}.attended.Freq,'aOrder',NaN,...
                'aName',SessionStructureExp{currtrial}.attended.Name,'aCategories',SessionStructureExp{currtrial}.attended.TypeCode,...
                'aVoice',SessionStructureExp{currtrial}.attended.Voice,'aSide',SessionStructureExp{currtrial}.attended.Side,...
                'aPresented',SessionStructureExp{currtrial}.attended.Presented,...
                'iCode',SessionStructureExp{currtrial}.ignored.Code,'iFreq',SessionStructureExp{currtrial}.ignored.Freq,...
                'iOrder',NaN,'iName',SessionStructureExp{currtrial}.ignored.Name,...
                'iCategories',SessionStructureExp{currtrial}.ignored.TypeCode,'iVoice',SessionStructureExp{currtrial}.ignored.Voice,...
                'iSide',SessionStructureExp{currtrial}.ignored.Side,'jCode',SessionStructureExp{currtrial}.attended.Code,'iPresented',SessionStructureExp{currtrial}.ignored.Presented);
        end
        
        
        son = transpose(son);
%         spaz_son = SpatializeSound(son(1:2,:),fs,0.0003);
%          son(1:2,:) =  spaz_son(1:length(son),1:2)';
         
        son(3,:)=1000000*son(3,:);
        clear('Stim');
        
        if mycomp==0
            if mod(currtrial+1-startingTrial,2)==1
                if currtrial>startingTrial
                    
                    PsychPortAudio('DeleteBuffer',trialbuffer1);
                    clear('trialbuffer1');
                    if exist('trialbuffer2')==0
                        warning('fix me!')
                    trialbuffer1 = PsychPortAudio('CreateBuffer', pahandle, son);
                    currentbuffer = trialbuffer1;
                    else
                    currentbuffer = trialbuffer2;
                    end
                end
                if mycomp==0
                    trialbuffer1 = PsychPortAudio('CreateBuffer', pahandle, son);
                end
                nextbuffer = trialbuffer1;
            elseif mod(currtrial+1-startingTrial,2)==0
                if currtrial>startingTrial+1
                    if exist('trialbuffer2')==0
                    else
                        PsychPortAudio('DeleteBuffer',trialbuffer2);
                    end
                end
                currentbuffer = trialbuffer1;
                trialbuffer2 = PsychPortAudio('CreateBuffer', pahandle, son);
                nextbuffer = trialbuffer2;
                
                clear('son','trial');
            end
        end
        set(h9,'String',sprintf('Trial %d',currtrial-1));
        if mycomp==0
            handleStatus = PsychPortAudio('GetStatus',pahandle);
        end
        
        if currtrial==startingTrial
            looping = 0;
        else
            looping = 1;
        end
        
        rep=0;
        
        audioStopped = 0;
        clock = 0;
        ite = 1;
        down = 0;
        % Wait until the trial end
        if mycomp==0
            while looping
%                 figure(ctrlPanel);
                handleStatus = PsychPortAudio('GetStatus',pahandle);
                if ~handleStatus.Active
                    if ~audioStopped
                        stopTimes(currevent,1) = GetSecs();
                        if AcquisitionFlag==1
                            [status, error] = NetStationPC_Debugged('Event','endT',stopTimes(currevent,1),0.001,'Subj',str2num(Subject.ID),'Trai',StimCarac.TrialType,'nTri',countingTrial,...
                                'aCod',StimCarac.aCode,'aVoi',StimCarac.aVoice,'aSid',StimCarac.aSide,'aOrd',StimCarac.aOrder,...%'aTag',StimCarac.afTag,'aAmp',StimCarac.aAmpTag*100,...
                                'iCod',StimCarac.iCode,'iVoi',StimCarac.iVoice,'iSid',StimCarac.iSide,'iOrd',StimCarac.iOrder);%'iTag',StimCarac.ifTag,'iAmp',StimCarac.iAmpTag*100);
                        end
                        audioStopped = 1;
                        set(h10,'String','Jitter');
                    end

                    PsychPortAudio('Stop',pahandle,1);
                    set(h5,'String','Inactive');
                    figure(ctrlPanel);
                end
                looping = ~audioStopped;
            end
        end
        set(h10,'String','');
%         figure(ctrlPanel);
        
        % Jitter
        if rep==1 || ...
                (currtrial-1==2*nbehavTraining+2*nwakeTest+nwake+nEEGTraining && rep==0)...
                (currtrial==1+2*nbehavTraining+2*nwakeTest+nfalling+nN2+nN3+nwake+nEEGTraining && rep==0)
            % Replay button has been pressed or subject has finished the
            % Wake phase
            set(h10,'String','Waiting for you to click on Continue');
            pausing = GetSecs;
            if AcquisitionFlag==1
                [status, error] = NetStationPC_Debugged('Event','Paus',pausing,0.001,'Subj',str2num(Subject.ID),'nTri',0);
            end
            figure(ctrlPanel);
            uiwait(ctrlPanel);
            ContinueTime = GetSecs;
            if AcquisitionFlag==1
                [status, error] = NetStationPC_Debugged('Event','Cont',ContinueTime,0.001,'Subj',str2num(Subject.ID),'nTri',0);
            end
            set(h10,'String','');
        elseif currtrial~=startingTrial 
            while GetSecs-stopTimes(currevent,1)<endTrialJitter(countingTrial-1)
            end
        end
        
        checkStatus=PsychPortAudio('GetStatus',pahandle);
        while checkStatus.Active
            checkStatus=PsychPortAudio('GetStatus',pahandle);
        end
        if currtrial~=startingTrial
            currevent = currtrial;
        end
        
        
        if mycomp==0
            % Create a new schedule and add to it the new trial
            PsychPortAudio('UseSchedule', pahandle, 2);
        end
        if rep==0
            if mycomp==0
                PsychPortAudio('AddToSchedule', pahandle, nextbuffer, 1, 0, [], 1);
            end
            if exitK==1
                playedTrial{countingTrial,1} = currtrial;
            end
        elseif rep==1
            if mycomp==0
                PsychPortAudio('AddToSchedule', pahandle, currentbuffer, 1, 0, [], 1);
            end
            if exitK==1
                playedTrial{countingTrial,1} = currtrial-1;
            end
        end
        currtrial=currtrial+1;    

        if mycomp==0
           
            % Wait until any sound store in schedule stop

            PsychPortAudio('Stop', pahandle , 1);
        end
        
       

        % Play the new trial
        if exitK==1
%             sound(son(1:2,:)',fs)
            if mycomp==0
                PsychPortAudio('Start', pahandle , 1 , 0 , 0, [], 0);
            end
            startTimes(currevent,1) = GetSecs();
            StimCarac = nextStimCarac;
            
            %%%%%% Sending Event Information to NetSation
            if AcquisitionFlag==1
                [status, error] = NetStationPC_Debugged('Event','newT',startTimes(currevent,1),0.001,'Subj',str2num(Subject.ID),'Trai',training,'nTri',countingTrial,...
                    'aCod',StimCarac.aCode,'aVoi',StimCarac.aVoice,'aSid',StimCarac.aSide,'aOrd',StimCarac.aOrder,...
                    'iCod',StimCarac.iCode);
            end
        end
     
        fprintf('New event: newT %g, Subj %s nTri %g \n\t\attended %d\n\t\ignored %d\n',num2str(startTimes(currevent)),Subject.ID,playedTrial{countingTrial,1},...
            num2str(StimCarac.aCode),...%StimCarac.afTag,StimCarac.aAmpTag,...dd
            num2str(StimCarac.iCode))%StimCarac.ifTag,StimCarac.iAmpTag)
        
        trialType{countingTrial,1} = StimCarac.TrialType;
        trialaCode{countingTrial,1} = StimCarac.aCode;
        trialaVoice{countingTrial,1} = StimCarac.aVoice;
        trialaSide{countingTrial,1} = StimCarac.aSide;
        trialaCategories{countingTrial,1} = StimCarac.aCategories;
        trialaName{countingTrial,1} = StimCarac.aName;
        trialaFreq{countingTrial,1} = StimCarac.aFreq;
        trialaPresented{countingTrial,1} = StimCarac.aPresented;
        trialjCode{countingTrial,1} = StimCarac.jCode;
        trialaOrder{countingTrial,1} = StimCarac.aOrder;
        
        trialiCode{countingTrial,1} = StimCarac.iCode;
        trialiVoice{countingTrial,1} = StimCarac.iVoice;
        trialiSide{countingTrial,1} = StimCarac.iSide;
        trialiCategories{countingTrial,1} = StimCarac.iCategories;
        trialiName{countingTrial,1} = StimCarac.iName;
        trialiFreq{countingTrial,1} = StimCarac.iFreq;
        trialiPresented{countingTrial,1} = StimCarac.iPresented;
        trialiOrder{countingTrial,1} = StimCarac.iOrder;
        
        
        % Update of the control panel
        set(h5,'String','Playing');
        set(h4,'String','Trial');
        set(h4,'String',num2str(currevent));
        set(h17,'String',StimCarac.TrialType);
        set(h18,'String',StimCarac.aName);
        set(h19,'String',StimCarac.iName);
        set(h20,'String',StimCarac.aVoice);
        set(h21,'String',StimCarac.aSide);
        figure(ctrlPanel);
        
        
        %wait until the trial has finished
        if mycomp==0
            audioStopped=0;
            while looping
                figure(ctrlPanel);
                handleStatus = PsychPortAudio('GetStatus',pahandle);
                
                if ~handleStatus.Active
                    if ~audioStopped
                        stopTimes(currevent,1) = GetSecs();
                        if AcquisitionFlag==1
                            [status, error] = NetStationPC_Debugged('Event','endT',stopTimes(currevent,1),0.001,'Subj',str2num(Subject.ID),'Trai',training,'nTri',countingTrial,...
                                'aCod',StimCarac.aCode,'aVoi',StimCarac.aVoice,'aSid',StimCarac.aSide,'aOrd',StimCarac.aOrder,...%'aTag',StimCarac.afTag,'aAmp',StimCarac.aAmpTag*100,...
                                'iCod',StimCarac.iCode,'iVoi',StimCarac.iVoice,'iSid',StimCarac.iSide,'iOrd',StimCarac.iOrder);%'iTag',StimCarac.jfTag,'iAmp',StimCarac.iAmpTag*100);
                        end
                        audioStopped = 1;
                        set(h10,'String','Jitter');
                    end
                    PsychPortAudio('Stop',pahandle,1);
                    set(h5,'String','Inactive');
                    figure(ctrlPanel);
                end
                looping = ~audioStopped;
            end
        end
           %%
                
        if ismember(currtrial,TestTrialsAwake) || currtrial>nTotWake+nTotSleep
            
            waitingContinue=1;
            fprintf('press the key')
            checkStart=GetSecs;
            while waitingContinue
                [down,timingKey,kCode]=KbCheck;
                if down==1 && sum(kCode)==1
                    keyboardKey = find(kCode);
                    waitingContinue = 0;
                    if AcquisitionFlag==1
                        [status, error] = NetStationPC_Debugged('Event','spcP',timingKey,0.001,'Subj',str2num(Subject.ID),'nTri',countingTrial);
                    end
                end
            end
                        
            disp(keyboardKey)
            disp(timingKey-checkStart)
            %     Response.StimCarac = SessionStructureExp{currtrial}.ignored.Name,'iCategories',SessionStructureExp{currtrial}.ignored.TypeCode,'iVoice';
            ResponseType{countingTrial,1} = keyboardKey;
            ResponseRT{countingTrial,1} = timingKey;%temps réponse - timing du son
            Responseevent{countingTrial,1} = countingTrial;
            ResponseType{countingTrial,1} = StimCarac.TrialType;
            ResponseoCode{countingTrial,1} = StimCarac.aCode;
            ResponsenCode{countingTrial,1} = StimCarac.iCode;
            ResponsejCode{countingTrial,1} = StimCarac.jCode;
            ResponseVoice{countingTrial,1} = StimCarac.aVoice;
            ResponseoOrder{countingTrial,1} = StimCarac.aOrder;
            ResponseiOrder{countingTrial,1} = StimCarac.iOrder;
            ResponseaCategories{countingTrial,1} = StimCarac.aCategories;
            ResponseaFreq{countingTrial,1} = StimCarac.aFreq;
            ResponseiCategories{countingTrial,1} = StimCarac.iCategories;
            ResponseaPresented{countingTrial,1} = StimCarac.aPresented;
            ResponseiFreq{countingTrial,1} = StimCarac.iFreq;
            ResponseiPresented{countingTrial,1} = StimCarac.iPresented;
            
        else 
         ResponseType{countingTrial,1} = nan;
            ResponseRT{countingTrial,1} = nan;%temps réponse - timing du son
            Responseevent{countingTrial,1} = nan;
            ResponseType{countingTrial,1} = nan;
            ResponseoCode{countingTrial,1} = nan;
            ResponsenCode{countingTrial,1} = nan;
            ResponsejCode{countingTrial,1} = nan;
            ResponseVoice{countingTrial,1} = nan;
            ResponseoOrder{countingTrial,1} = nan;
            ResponseiOrder{countingTrial,1} = nan;
            ResponseaCategories{countingTrial,1} = nan;
            ResponseaFreq{countingTrial,1} = nan;
            ResponseiCategories{countingTrial,1} = nan;
            ResponseaPresented{countingTrial,1} = nan;
            ResponseiFreq{countingTrial,1} = nan;
            ResponseiPresented{countingTrial,1} = nan;
        end
        
        
        if training
            if mycomp==0
                % End trial sound
                PsychPortAudio('UseSchedule', pahandle, 2);
                PsychPortAudio('AddToSchedule', pahandle, endTrialBuffer, 1, 0, [], 1);
                PsychPortAudio('Start', pahandle , 1 , 0 , 0, [], 0);
            end
            [beep,samplingRate] = MakeBeepAP(1000,0.01,fs);
            sound(beep*0.07)
            waitingContinue = 1;
            set(h10,'String','Waiting for the subject to press any key');
            figure(ctrlPanel);
            
            fprintf('PRESS A KEY TO CONTINUE!\n')
            while waitingContinue
                [down,secs,kCode]=KbCheck;
                if down==1 && sum(kCode)==1 %&& sum(ismember(find(kCode),AnswerKeys))>1
                    waitingContinue = 0;
                    if AcquisitionFlag==1
                        [status, error] = NetStationPC_Debugged('Event','spcP',secs,0.001,'Subj',str2num(Subject.ID),'nTri',0);
                    end
                end
            end
            set(h10,'String','');
            timer = GetSecs;
            while GetSecs-timer<1
            end
        end
        %%
        if Changelist==0
            %loop on the same list for falling, N2 or N3
           if currtrial == nTotWake
                randpermfalling = randperm(length(listnames.falling.list));
                SummaryExp(listnames.falling.list,:) = SummaryExp(listnames.falling.list(randpermfalling),:);
                SessionStructureExp(listnames.falling.list,:) = SessionStructureExp(listnames.falling.list(randpermfalling),:);
                currtrial=currtrial-nfalling+1;
                listnames.falling.active =currtrial;
                currevent=currtrial;
            elseif currtrial == nTotWake+nN2
                randpermSleepN2 = randperm(length(listnames.SleepN2.list));
                SummaryExp(listnames.SleepN2.list,:) = SummaryExp(listnames.SleepN2.list(randpermSleepN2),:);
                SessionStructureExp(listnames.SleepN2.list,:) = SessionStructureExp(listnames.SleepN2.list(randpermSleepN2),:);
                currtrial=currtrial-nN2+1;
                listnames.SleepN2.active =currtrial;
                 currevent=currtrial;
           elseif currtrial == nTotWake+nN2+nN3
                randpermSleepN3 = randperm(length(listnames.SleepN3.list));
                SummaryExp(listnames.SleepN3.list,:) = SummaryExp(listnames.SleepN3.list(randpermSleepN3),:);
                SessionStructureExp(listnames.SleepN3.list,:) = SessionStructureExp(listnames.SleepN3.list(randpermSleepN3),:);
                currtrial=currtrial-nN3+1;
                listnames.SleepN3.active =currtrial;
                 currevent=currtrial;
            end
        end
%          if rep == 1 %replay
%                 currtrial=currtrial;
%          else:end
         
        if currtrial>nTotTrials
            running = 0;
        end
         tic;
        while  toc <= 5
        end
            
        clock = 0;
        if exitK==1
            looping = 1;
            % Small jitter so that the AudioPort update its informations
            while GetSecs-startTimes(currevent,1)<1
            end
        else
            looping = 0;
        end
        
        
        % Wait until the last question is answered
        audioStopped = 0;
        clock = 0;
        ite = 1;
        down = 0;
        % Wait until the trial end
        if mycomp==0
            while looping
                figure(ctrlPanel);
                
                handleStatus = PsychPortAudio('GetStatus',pahandle);
                if ~handleStatus.Active
                    if ~audioStopped
                        stopTimes(currevent,1) = GetSecs();
                        if AcquisitionFlag==1
                            [status, error] = NetStationPC_Debugged('Event','endT',stopTimes(currevent,1),0.001,'Subj',str2num(Subject.ID),'Trai',training,'nTri',countingTrial,...
                                'aCod',StimCarac.aCode,'aVoi',StimCarac.aVoice,'aSid',StimCarac.aSide,'aOrd',StimCarac.aOrder,...%'iCategories',SessionStructureExp{currtrial}.ignored.TypeCode,'iVoice'.aCode,...%'aTag',StimCarac.afTag,'aAmp',StimCarac.aAmpTag*100,...
                                'iCod',StimCarac.iCode,'iVoi',StimCarac.iVoice,'iSid',StimCarac.iSide,'iOrd',StimCarac.iOrder);%'iTag',StimCarac.jfTag,'iAmp',StimCarac.iAmpTag*100);
                        end
                        audioStopped = 1;
                        set(h10,'String','Jitter');
                    end

                    PsychPortAudio('Stop',pahandle,1);
                    set(h5,'String','Inactive');
                    figure(ctrlPanel);
                    
                    looping = ~audioStopped;
%                 else
%                     if ~audioStopped
%                         stopTimes(currevent,1) = GetSecs();
%                         audioStopped = 1;
%                         set(h10,'String','Jitter');
%                     end
%                     set(h5,'String','Inactive');
%                     figure(ctrlPanel);
%                     looping = ~audioStopped;
                end
            end
        end

                
%     catch
%         fprintf('Une erreur est survenue durant le script, les donn\E9es ont \E9t\E9 sauvegard\E9es, vous pourrez reprendre l\E0 o\F9 le sujet s''est arr\EAt\E9. Tentez de trouver la faille.')
%     end
    countingTrial = countingTrial+1;
    if mycomp==0

        PsychPortAudio('Stop', pahandle, 2);
    end
    if exitK==1
        stopFunctionTime(currevent,1) = GetSecs();
    end
    
    if numel(startTimes)>numel(stopTimes)
        stopTimes(currevent,1) = GetSecs;
        stopFunctionTime(currevent,1) = GetSecs;
    end
    
    Timings.expeStart = expeStart;
    Timings.startTime = num2cell(startTimes);
    Timings.stopTime = num2cell(stopTimes);
    Timings.stopFunctionTime = num2cell(stopFunctionTime);
end
%%
TrialsCaracs = struct(  'Trials',playedTrial,...
    'TrialType',trialType,...
    'aCode',trialaCode,'aOrder',trialaOrder,'aVoice',trialaVoice,'aSide',trialaSide,'aName',trialaName,'aCategories',trialaCategories,'aFreq',trialaFreq,'aPresented',trialaPresented,...%'aTag',trialaTag,'aAmp',trialaAmp,...
    'iCode',trialiCode,'iOrder',trialiOrder,'iVoice',trialiVoice,'iSide',trialiSide,'iName',trialiName,'iCategories',trialiCategories,'iFreq',trialiFreq,'jCode',trialjCode,'iPresented',trialiPresented);%'iTag',trialiTag,'iAmp',trialiAmp);

ResponsesCaracs = struct('Type',ResponseType,'RT',ResponseRT,...
    'event',Responseevent,'TrialType',trialType,...
    'aCode',ResponseoCode,'aOrder',ResponseoOrder,'Voice',ResponseVoice,'oCategories',ResponseaCategories,'oFreq',ResponseaFreq,'oPresented',ResponseaPresented,...%'aTag',trialaTag,'aAmp',trialaAmp,...
    'iCode',ResponsenCode,'nOrder',ResponseiOrder,'nCategories',ResponseiCategories,'nFreq',ResponseiFreq,'jCode',ResponsejCode,'nPresented',ResponseiPresented);%'iTag',trialiTag,'iAmp',trialiAmp);

if mycomp==0
    PsychPortAudio('Close',pahandle);
end
close(ctrlPanel);

%%
%%%%%%                  End Recordings
if AcquisitionFlag==1
    fprintf('STOPPING NET STATION\n'); fprintf(' \n')
    
    [status, error] = NetStationPC_Debugged('event','END!',GetSecs,0.001,'Subj',Subject.ID);
    WaitSecs(10);  %to avoid nul data when filter is used on the EEG data
    [status, error] = NetStationPC_Debugged('stoprecording');
end
Timings.EndExpe=GetSecs;

%%
%%%%%%                  Save BehavioralResults and Close
fprintf('Saving data\n'); fprintf(' \n')
savePath=[dataPath filesep expeName Subject.ID];
mkdir(savePath)
saveName=[savePath filesep 'Result_' expeName  Subject.ID];
save(saveName,'Subject','Timings','TrialsCaracs','ResponsesCaracs','SessionStructureExp')
