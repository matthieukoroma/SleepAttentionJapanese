%call a sound to be played for testing

fs=44100;
snrWN=0.4;
raspiness=2;
for i=1:1
    namefile=['home/lscp00/Codes/Projects/attention japanese/SleepAttentionJapanese/Expe_Folder/Expe_Material/testStims/Lf',num2str(i),'LKFS.mat'];
    namedest=['/home/lscp00/Codes/Projects/attention japanese/SleepAttentionJapanese/Expe_Folder/Expe_Material/testStimswav/Lf',num2str(i),'LKFS.wav'];
    load(namefile)
    son=Stim.Mat(:,1:2);
    %     whitenoise1=addAWGN(son(:,1),snrWN);
    %     whitenoise2=addAWGN(son(:,1),snrWN);
    %     son2=son+[whitenoise1,whitenoise2];
    audiowrite(namedest,son(:,1),fs)
    [status, result] = system (['Praat.exe raspiness.praat' namedest])
%     sound(son,fs)
    b=input('r')
end