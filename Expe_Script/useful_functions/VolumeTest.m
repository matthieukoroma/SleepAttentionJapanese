
clear all

fs = 44100;

nChannels=4; %Number of channels
nDevice=30; % AudioFire12:DeviceID=30
runMode=0;

initialVolume = 0.4;

rand('state',sum(100*clock));

InitializePsychSound(1);

pahandle = PsychPortAudio('Open',nDevice, 1, 2, fs, nChannels);
PsychPortAudio('Volume',pahandle,initialVolume);
PsychPortAudio('UseSchedule', pahandle, 1);
PsychPortAudio('RunMode', pahandle, 0);

soundPath = 'C:\Data\SleepLexicalDecision2\Expe_Folder\Expe_Material\secondversion';
soundDir = dir(soundPath);
sounds = {soundDir.name};

[finalSound,sr] = createWordStringTest(5,soundPath,4);

soundBuffer = PsychPortAudio('CreateBuffer',pahandle,finalSound);

% Control Panel
ctrlPanel = figure();
h1 = uicontrol('Position', [20 20 80 40], 'String', 'Sound', 'Callback',...
    sprintf([   'PsychPortAudio(''UseSchedule'', pahandle, 2);\n',...
                'PsychPortAudio(''AddToSchedule'', pahandle, soundBuffer, 1, 0, [], 1);\n',...
                'PsychPortAudio(''Start'', pahandle , 1 , 0 , 0, [], 0);\n'...
                'PsychPortAudio(''Stop'', pahandle , 1);\n']));
h2 = uicontrol('Position', [120 20 80 40], 'String', 'New Sound', 'Callback',...
    sprintf([   '[finalSound,sr] = createWordStringTest(5,soundPath,4);\n',...
                'soundBuffer = PsychPortAudio(''CreateBuffer'',pahandle,finalSound);']));
h5 = uicontrol('Position', [300 100 80 40],'Style','text','String',num2str(initialVolume));
h3 = uicontrol('Style','slider','Position', [220 20 20 100],'Value',initialVolume,'Min',0,'Max',1.0,'SliderStep',[0.01,0.1],...
    'Callback',sprintf([    'initialVolume = get(h3,''Value'');\n'...
                            'PsychPortAudio(''Volume'',pahandle,initialVolume);\n',...
                            'set(h5,''String'',num2str(initialVolume));']));
uiwait(ctrlPanel);

PsychPortAudio('Close',pahandle);
