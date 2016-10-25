for i=1:10
    for voice=['f','m'];
    load(['/home/lscp00/Codes/Projects/attention_japanese/SleepAttentionJapanese/Expe_Folder/Expe_Material/testStims/L',voice,num2str(i),'LKFS.mat'])
    son=Stim.Mat(:,1:2);
    audiowrite(['L',voice,num2str(i),'LKFSrasp.wav'],son,44100)
    end
end