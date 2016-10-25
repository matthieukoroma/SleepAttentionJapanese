% dirname = 'home\\lscp00\\Matthieu\\SleepAttentionJapanese\\Expe_Folder\\Expe_Material';
dirname = '/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/';
langue={'f','j'};
voices = {'f','m'};

[Codes,Categories,NamesJ,NamesF,Freq,accrosslistpairs,attentionpairs] = load_stim_info(dirname);
%%
sr =  44100;% Stim.Caracs.fsample;
equa_tech ={'RMS'};%'RMS'};%'STL';% {'LFKS'};%
carrying_sentences = {'suivez_la_femme','suivez_lhomme','veut_dire'};
lkfs = -30;

load('/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/sons/TTLson');
%%

for ind=1:length(Codes);
    ind
    for voice=voices       %length(Codes)
        
%         cd([dirname,'stimsTrialsVsola/'])
%         
%         name_write= [voice,num2str(Codes(ind)),'.wav'];
%         name_write= [name_write{:}];
%         
%         %     load(name_write);
%         
%         name_load= [voice,num2str(Codes(ind)),'.mat'];
%         name_load= [name_load{:}];
%         load(name_load)
%         name_load2= [voice,num2str(Codes(ind)),'.wav'];
%         name_load2= [name_load2{:}];
%         
%         %     cd([dirname,'sonsTrialsVsola/'])
%         %     s=wavread(name_load2);
%         %     length(s)
%         s = Stim.Mat;
%         
%         Stim.Caracs.File = name_write;
%         Stim.Caracs.Voice = voice;
%         Stim.Caracs.Type = Categories(ind);
%         Stim.Caracs.NamesF = NamesF{ind};
%         Stim.Caracs.NamesJ = NamesJ{ind};
%         Stim.Caracs.Codes = Codes(ind);
%         Stim.Caracs.fsample = sr;
%         %put TTLsound on segments
%         
        %%
%         tic
%         %find silences to add triggers
%         vec0 = s(:,1)==0;
%         lens = length(s);
%         seg = sr*0.1;
%         nb_seg = floor(lens/seg);
%         plages0 = reshape(s(1:nb_seg*seg,1),nb_seg,seg);
%         sil = find(sum(plages0));
%         if length(sil)~=4
%             keyboard
%         end
%         for i=1:length(sil)
%             TTLi=find(vec0(sil(i):end),0,'first');
%             s(sil(i)+TTLi+[1:length(TTLson)],4)=TTLson;
%             TTLi=find(vec0(sil(i):-1:1),0,'first');
%             s(sil(i)-TTLi+[1:length(TTLson)],4)=TTLson;
%         end
%         toc
%         A=findSilences(s(:,1),sr,'fixed',0,0.2);
%         if length(A)<4
%             keyboard
%         end
%         A=A(:)';
%         for i=A
%            s(i+[0:length(TTLson)-1],4)=TTLson;
%         end
%         if ~ismember(length(s),A)
%             s(end-length(TTLson)+1:end,4) = TTLson;
%         end
%         Stim.Mat=s;
%         save(name_load,'Stim')  
                Stim.Mat=[];
        %%
        for techn=equa_tech
            %         techn = tech{:};
            if strcmp(techn,'LKFS')
                factor = normalise_loudness(s(:,1)', sr, lkfs);
                Stim.Caracs.LKFSfactor = factor;
                Stim.Caracs.LKFSval = lkfs;
                Stim.Mat(:,1:2)=s(:,1:2)*factor;
            elseif strcmp(techn,'RMS')
                Stim.Mat(:,1)=s(:,1)./rms(s(:,1))*0.07;
                Stim.Mat(:,2)=Stim.Mat(:,1);
                %             Stim.Mat(:,1:2)=Stim.Mat(:,1:2);
            elseif strcmp(techn,'STL')
                %get info for sonic normalization
                %             res = Loudness_TimeVaryingSound_Moore(s(:,1), sr, 'mic',0);
                %             Stim.Caracs.STLmax = res.STLmax;
                Stim.Mat(:,1)=s(:,1)/STLmax;
                Stim.Mat(:,2)=Stim.Mat(:,1);
            end
            
%             dirfile = [dirname,'stimsTrials',techn];
%             cd([dirfile{:}])
%             save(name_load,'Stim')
            namf= [dirname,'sonsTrials',techn];
            namf=[namf{:}];
                    cd(namf)
                    wavwrite(Stim.Mat(:,1:2),sr,name_write);
            
            
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
                    factor = normalise_loudness(s(:,1), sr, lkfs);
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
                
                Stim.Mat(:,4)=Stim.Mat(:,3);
                namcd=[dirname,'stimsTest',techn];
                cd([namcd])
               
                save(name_load,'Stim')
                
                 namcd = [dirname,'sonsTest',techn];
                cd([namcd])
                %                 cd([dirname,'sonsTrials',techn])
                  wavwrite(s(:,1:2),sr,name_write);
                
           
%                 save(name_load,'Stim')
            end
            %                 cd([dirname,'sonsTrials',techn])
            %                 wavwrite(Stim.Mat(:,2),sr,name_write);
            
        end
        
        
        
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
        %             end
        
        %}
        
    end
end
%}
% 
% num=['','1','2','3','4'];
% cd([dirname,'stimsCarrying/'])
% for voice=voices
%     for cs=carrying_sentences
%         ok=1;
%         while ok==1;
%             numname = num(nu);
%             nu=nu+1;
%             try
%                 name_load= [voice,cs,numname,'.mat'];
%                 name_load= [name_load{:}];
%                 name_write= [voice,cs,numname,'.wav'];
%                 name_write= [name_write{:}];
%                 
%                 load(name_load);
%                 s = Stim.Mat;
%                 
%                 Stim.Caracs.File = name_write;
%                 Stim.Caracs.Voice = voice;
%                 Stim.Caracs.fsample = sr;
%                 
%                 for tech=equa_tech
%                     techn = tech{:};
%                     
%                     if strcmp(techn,'LKFS')
%                         factor = normalise_loudness(s(:,2), sr, lkfs);
%                         Stim.Caracs.LKFSfactor = factor;
%                         Stim.Caracs.LKFSval = lkfs;
%                         Stim.Mat(:,1:2)=s(:,1:2)*factor;
%                     elseif strcmp(techn,'RMS')
%                         Stim.Mat(:,1)=s(:,1)/rms(s(:,1));
%                         Stim.Mat(:,2)=Stim.Mat(:,1);
%                         Stim.Mat(:,1:2)=Stim.Mat(:,1:2)*0.07;
%                     elseif strcmp(techn,'STL')
%                         %get info for sonic normalization
%                         res = Loudness_TimeVaryingSound_Moore(s(:,1), sr, 'mic',0);
%                         Stim.Caracs.STLmax = res.STLmax;
%                         Stim.Mat(:,1)=s(:,1)/res.STLmax/0.07;
%                         Stim.Mat(:,2)=Stim.Mat(:,1);
%                     end
%                     
%                     if strcmp(techn,'veut_dire')
%                         cd([dirname,'stimsTest',techn])
%                         save(name_load,'Stim')
%                         
%                         cd([dirname,'sonsTest',techn])
%                         wavwrite(s,sr,name_write);
%                     else
%                         cd([dirname,'stimsTrials',techn])
%                         save(name_load,'Stim')
%                         
%                         cd([dirname,'sonsTrials',techn])
%                         wavwrite(s,sr,name_write);
%                     end
%                     
%                 end
%             catch
%                 ok=0;
%             end
%         end
%     end
% end
%     % STLmaxmean = mean(STLmaxs(:));
%     % for voice=1:length(voices)
%     %     for ind=1:length(Codes)
%     %         [s,sr] = wavread([dirname,'/sonsTrials/',voice,num2str(codes(ind)),'.wav']);
%     %         s=s*(STLmaxmean/STLmaxs(voice,ind))*0.07;
%     %         wavwrite(s,sr,[dirname,'/SonsTrialsLoudness/',voice,num2str(Codes(stim)),'.wav']);
%     %     end
%     % end
%     
%     %