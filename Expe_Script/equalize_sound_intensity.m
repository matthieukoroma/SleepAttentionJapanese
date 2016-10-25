% dirname = 'home\\lscp00\\Matthieu\\SleepAttentionJapanese\\Expe_Folder\\Expe_Material';
dirname = '/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/';
langue={'f','j'};
voices = {'f','m'};

[Codes,Categories,NamesJ,NamesF,Freq,accrosslistpairs,attentionpairs] = load_stim_info(dirname);
%%
% sr =   Stim.Caracs.fsample;
equa_tech ={'LKFS'};%'STL';% {};%
carrying_sentences = {'suivez_la_femme','suivez_lhomme','veut_dire'};
lkfs = -30;
%%
for ind=205;%1:length(Codes)
ind
for voice=voices(2)
   

           %length(Codes)
        for techn=equa_tech
            namf=[dirname,'stimsTrialsVsola',techn{:}];
%                         namf=[namf{:}];
            cd(namf)
            
            name_write= [voice,num2str(Codes(ind)),'.wav'];
            name_write= [name_write{:}];
            
            %     load(name_write);
            
            name_load= [voice,num2str(Codes(ind)),'.mat'];
            name_load= [name_load{:}];
            load(name_load)
            name_load2= [voice,num2str(Codes(ind)),'.wav'];
            name_load2= [name_load2{:}];
            
            %     cd([dirname,'sonsTrialsVsola/'])
            %     s=wavread(name_load2);
            %     length(s)
            s = Stim.Mat;
            
            Stim.Caracs.File = name_write;
            Stim.Caracs.Voice = voice;
            Stim.Caracs.Type = Categories(ind);
            Stim.Caracs.NamesF = NamesF{ind};
            Stim.Caracs.NamesJ = NamesJ{ind};
            Stim.Caracs.Codes = Codes(ind);
            Stim.Caracs.fsample = sr;
            
            
            %         techn = tech{:};
            if strcmp(techn,'LKFS')
                
                factor = normalise_loudness(s(:,1), sr, lkfs);
                Stim.Caracs.LKFSfactor = factor;
                Stim.Caracs.LKFSval = lkfs;
                Stim.Mat(:,1:2)=s(:,1:2)*factor;
            elseif strcmp(techn,'RMS')
                Stim.Mat(:,1) = 0.07*s(:,1)./rms(s(:,1));
                Stim.Mat(:,2)=Stim.Mat(:,1);
                %             Stim.Mat(:,1:2)=Stim.Mat(:,1:2);
            elseif strcmp(techn,'STL')
                %get info for sonic normalization
                %             res = Loudness_TimeVaryingSound_Moore(s(:,1), sr, 'mic',0);
                %             Stim.Caracs.STLmax = res.STLmax;
                Stim.Mat(:,1)=s(:,1)/STLmax;
                Stim.Mat(:,2)=Stim.Mat(:,1);
            end
            
            namf = [dirname,'stimsTrialsVsola',techn{:}];
%             namf=[namf{:}];
            cd(namf)
            save(name_load,'Stim')
            %         namf= [dirname,'sonsTrials',techn];
            %         cd(namf)
            %         wavwrite(Stim.Mat(:,1:4),sr,name_write);
            
            
        end
    end
end
%{
for lang=langue
    
    
    %             try
    %
    %                     cd([dirname,'stimsTest/'])
    %         name_load= [voice,num2str(Codes(ind)),'.mat'];
    %         name_load= [name_load{:}];
    %         name_write= [voice,num2str(Codes(ind)),'.wav'];
    %         name_write= [name_write{:}];
    %
    
    cd([dirname,'stimsTest/'])
    name_load= [voice,lang,num2str(Codes(ind)),'.mat'];
    name_load= [name_load{:}];
    name_write= [voice,lang,num2str(Codes(ind)),'.wav'];
    name_write= [name_write{:}];
    
    load(name_load);
    s = Stim.Mat;
    %
    %
    %         cd([dirname,'sonsTestRMS/'])
    %         % rms normalization
    %         wavwrite(s,sr,name_write);
    %            cd([dirname,'stimsTestRMS/'])
    
    Stim.Caracs.File = name_write;
    Stim.Caracs.Voice = voice;
    Stim.Caracs.Langue = lang;
    Stim.Caracs.Type = Categories(ind);
    Stim.Caracs.NamesF = NamesF{ind};
    Stim.Caracs.NamesJ = NamesJ{ind};
    Stim.Caracs.Codes = Codes(ind);
    Stim.Caracs.fsample = sr;
    
    for tech=equa_tech
        techn = tech{:};
        
        if strcmp(techn,'LKFS')
            factor = normalise_loudness(s(:,2), sr, lkfs);
            Stim.Caracs.LKFSfactor = factor;
            Stim.Caracs.LKFSval = lkfs;
            Stim.Mat(:,1:2)=s(:,1:2)*factor;
        elseif strcmp(techn,'RMS')
            Stim.Mat(:,1)=s(:,1)/rms(s(:,1));
            Stim.Mat(:,2)=Stim.Mat(:,1);
            Stim.Mat(:,1:2)=Stim.Mat(:,1:2)*0.07;
        elseif strcmp(techn,'STL')
            %get info for sonic normalization
            res = Loudness_TimeVaryingSound_Moore(s(:,1), sr, 'mic',0);
            Stim.Caracs.STLmax = res.STLmax;
            Stim.Mat(:,1)=s(:,1)/res.STLmax*0.07;
            Stim.Mat(:,2)=Stim.Mat(:,1);
        end
        
        cd([dirname,'stimsTest',techn])
        save(name_load,'Stim')
        
        %                 cd([dirname,'sonsTrials',techn])
        %                 wavwrite(s,sr,name_write);
    end
    
    %             cd([dirname,'stimsTests',techn])
    %             save(name_load,'Stim')
    
    %             cd([dirname,'sonsTrials',techn])
    %             wavwrite(s,sr,name_write);
    
    
    %         Stim.Caracs.compression = tsmfactor(ind,voice);
    
    %             %get info for sonic normalization
    %             res = Loudness_TimeVaryingSound_Moore(s(:,1), sr, 'mic',0);
    %             Stim.Caracs.MooreTVS = res;
    %             s(:,1)=s(:,1)/res.STLmax/0.07;
    %             s(:,2)=s(:,1);
    %             cd([dirname,'sonsTestSTL/'])
    %             wavwrite(s,sr,name_write);
    %             cd([dirname,'stimsTestSTL/'])
    %             save(name_load,'Stim')
    %                 pitch = res.ERB ???
    %             catch
    %                 continue
                end
end
%{
        cd([dirname,'stimsCarrying/'])
        for cs=carrying_sentences
            name_load= [voice,cs,'.mat'];
            name_load= [name_load{:}];
            name_write= [voice,cs,'.wav'];
            name_write= [name_write{:}];
            
            load(name_load);
            s = Stim.Mat;
            
            Stim.Caracs.File = name_write;
            Stim.Caracs.Voice = voice;
            Stim.Caracs.fsample = sr;
            
            for tech=equa_tech
                techn = tech{:};
                
                if strcmp(techn,'LKFS')
                    factor = normalise_loudness(s(:,2), sr, lkfs);
                    Stim.Caracs.LKFSfactor = factor;
                    Stim.Caracs.LKFSval = lkfs;
                    Stim.Mat(:,1:2)=s(:,1:2)*factor;
                elseif strcmp(techn,'RMS')
                    Stim.Mat(:,1)=s(:,1)/rms(s(:,1));
                    Stim.Mat(:,2)=Stim.Mat(:,1);
                    Stim.Mat(:,1:2)=Stim.Mat(:,1:2)*0.07;
                elseif strcmp(techn,'STL')
                    %get info for sonic normalization
                    res = Loudness_TimeVaryingSound_Moore(s(:,1), sr, 'mic',0);
                    Stim.Caracs.STLmax = res.STLmax;
                    Stim.Mat(:,1)=s(:,1)/res.STLmax/0.07;
                    Stim.Mat(:,2)=Stim.Mat(:,1);
                end
                
                if strcmp(techn,'veut_dire')
                    cd([dirname,'stimsTest',techn])
                    save(name_load,'Stim')
                    
                    cd([dirname,'sonsTest',techn])
                    wavwrite(s,sr,name_write);
                else
                    cd([dirname,'stimsTrials',techn])
                    save(name_load,'Stim')
                    
                    cd([dirname,'sonsTrials',techn])
                    wavwrite(s,sr,name_write);
                end
            end
%}
end
%}

% STLmaxmean = mean(STLmaxs(:));
% for voice=1:length(voices)
%     for ind=1:length(Codes)
%         [s,sr] = wavread([dirname,'/sonsTrials/',voice,num2str(codes(ind)),'.wav']);
%         s=s*(STLmaxmean/STLmaxs(voice,ind))*0.07;
%         wavwrite(s,sr,[dirname,'/SonsTrialsLoudness/',voice,num2str(Codes(stim)),'.wav']);
%     end
% end

%