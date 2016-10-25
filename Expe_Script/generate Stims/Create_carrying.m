dirname = '/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/';

voices = {'f','m'};

sentec = {'suivez_la_femme','suivez_lhomme'};%,'veut_dire'};%'veut_dire'};
load('/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/sons/TTLson');
dirsons = [dirname,'sonsCarrying'];
dirstim= [dirname,'stimsCarrying'];
[Codes,Categories,NamesJ,NamesF,Freq,accrosslistpairs,attentionpairs] = load_stim_info(dirname) ;
sr =  44100;% Stim.Caracs.fsample;
techn ={'RMS','LKFS'};%
lkfs = -30;
Stim=[];
%[2,5,6,4]
%%
combi_suiv = perms(1:3);

for ind=sentec(2)
    for voice=voices(2)
        for numc=4;%1:size(combi_suiv,1)
            stim_sentence=[];
            for num=1:size(combi_suiv,2)
                
                dname = [dirname,'sons/',voice,'/carrying_sentences/'];
                dname = [dname{:}];
                cd(dname)
                
                [s,sr] = wavread([ind{:},num2str(combi_suiv(numc,num)),'.wav']);
                
                sTTL = [TTLson;zeros(length(s)-length(TTLson),1)];
                if num==size(combi_suiv,2)
                    silence = zeros(sr*0.3,1);
                else
                    silence = zeros(sr*0.3,1);
                end
                silTTL = [TTLson;zeros(length(silence)-length(TTLson),1)];
                
                
                stim_sentence=[stim_sentence;[s(:),s(:),zeros(length(s),1),sTTL];[silence,silence,silence,silTTL]];
                
            end
            stim_seg = [1,length(stim_sentence)];
            stim_sentence(1:length(TTLson),3)=TTLson;
            %             cd(dirsons);
            name_write = [voice,ind];%
            name_write = [name_write{:}];
            %             wavwrite(stim_sentence,sr,[name_write,'.wav']);
            
            for tech=techn
                clear Stim
                s=stim_sentence;
                cd([dirstim,tech{:}]);
                if strcmp(tech,'LKFS')
                    factor = normalise_loudness(s(:,1)', sr, lkfs);
                    Stim.Caracs.LKFSfactor = factor;
                    Stim.Caracs.LKFSval = lkfs;
                    Stim.Mat(:,1:2)=s(:,1:2)*factor;
                elseif strcmp(tech,'RMS')
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
                
                Stim.Caracs.File = name_write;
                Stim.Caracs.Voice = voice{:};
                Stim.Caracs.Type = ind{:};
                Stim.Caracs.order = combi_suiv(numc,:);
                Stim.Caracs.NamesF = ind{:};
                Stim.Caracs.NamesJ = ind{:};
                Stim.Caracs.Code = num;
                Stim.Caracs.stim_seg = stim_seg;
                Stim.Caracs.fsample = sr;
                Stim.Mat = stim_sentence;
                save(name_write,'Stim')
                
                %                 cd([dirname,'stimsCarrying/'])
                %
                %                 name_write= [voice,ind,'.mat'];
                %                 name_write= [name_write{:}];
                %
                %             %     load(name_write);
                %         s =  stim_sentence;
                %         %             name_load= [voice,sentec{ind_p(ind)},'.mat'];
                %         %             name_load= [name_load{:}];
                %         for techn=equa_tech
                %             %         techn = tech{:};
                %             if strcmp(techn,'LKFS')
                %                 factor = normalise_loudness(s(:,2), sr, lkfs);
                %                 Stim.Caracs.LKFSfactor = factor;
                %                 Stim.Caracs.LKFSval = lkfs;
                %                 Stim.Mat(:,1:2)=s(:,1:2)*factor;
                %             elseif strcmp(techn,'RMS')
                %                 Stim.Mat(:,1)=s(:,1)./rms(s(:,1))*0.07;
                %                 Stim.Mat(:,2)=Stim.Mat(:,1);
                %                 %             Stim.Mat(:,1:2)=Stim.Mat(:,1:2);
                %             elseif strcmp(techn,'STL')
                %                 %get info for sonic normalization
                %                 %             res = Loudness_TimeVaryingSound_Moore(s(:,1), sr, 'mic',0);
                %                 %             Stim.Caracs.STLmax = res.STLmax;
                %                 Stim.Mat(:,1)=s(:,1)/STLmax;
                %                 Stim.Mat(:,2)=Stim.Mat(:,1);
                %             end
                
                %                             cd([dirname,'stimstrials',techn])
                %                             save(name_load,'Stim')
                %                             namf=[dirname,'stimsTest',techn];
                %                             namf=[namf{:}];
                %                             cd(namf)
                %                 save(name_write,'Stim')
                %             namf= [dirname,'sonsCarrying'];
                %                         namf=[namf{:R}];
                %             cd(namf)
                %             wavwrite(Stim.Mat(:,1:4),sr,name_write);
                
                
            end
        end
    end
end
%% equalize length vsola
%S
P = [round(44100/1000*7),round(44100/1000*15)];%the longest likely pitch period
for tech=techn
    ind_p=[1,2];
    len=[];
    for ind = 1:length(ind_p)
        for voice=1:length(voices)
%             for num=1:3
                %             namewave = [dirname,'stimsTrials/',voices{voice},num2str(Codes(ind)),'.wav'];
                %             [s{ind,voice},sr] = wavread(namewave);
                cd([dirstim,tech{:}])
                name_write = [voices(voice),sentec{ind_p(ind)}];%
                name_write = [name_write{:}];
                load([name_write,'.mat'])
                bounds{ind,voice} = Stim.Caracs.stim_seg;
                len(ind,voice)=sum(diff(bounds{ind,voice},1,2));
%             end
        end
    end
    lenmean=mean(len(:));
    %     tsmfactor1 = mean([len(1,1),len(2,2)]);
    %     tsmfactor2 = mean([len(2,1),len(1,2)]);
    %     tsmfactor = zeros(2);
    %
    %     tsmfactor1(1,1) = tsmfactor1;
    %     tsmfactor1(2,2) = tsmfactor1;
    %     tsmfactor2(2,1) = tsmfactor2;
    %     tsmfactor2(1,2) = tsmfactor2;
    tsmfactor = lenmean./len;
    
    Stim=[];
    %%
    P = [round(44100/1000*6),round(44100/1000*13)];%the longest likely pitch period
    for voice=1:length(voices)
        for ind = 1:length(ind_p)
            name_write = [voices(voice),sentec{ind_p(ind)}];%
            name_write = [name_write{:}];
            %         cd('/home/lscp00/Matthieu/SleepAttentionJapanese/ExpeMaterial/stimsCarrying/')
            cd([dirstim,tech{:}])
            clear Stim
            load([name_write,'.mat'])
            stim_m = Stim.Mat(:,1);
%             stimF =[];
%             for i=1:size(bounds{ind,voice},1)
%                 stimF = [stimF;stim_m(bounds{ind,voice}(i,1):min([length(stim_m),bounds{ind,voice}(i,2)]),1)];
%             end
%             keyboard
            svesola = vsola(stim_m, tsmfactor(ind,voice), P(voice));
            %             length(svesola)t
%             s=[];stim_seg=[];s=[];
%             %insert silence
%             trials_indices=zeros(size(bounds{ind,voice},1),2);
%             for i=1:size(bounds{ind,voice},1)
%                 trials_indices(i,:) = [max([floor(bounds{ind,voice}(i,1)*tsmfactor(ind,voice)),1]),min([length(svesola),ceil(bounds{ind,voice}(i,2)*tsmfactor(ind,voice))])];
%                 so=svesola(trials_indices(i,1):trials_indices(i,2));
%                 
                sTTL = [TTLson;zeros(length(svesola)-length(TTLson),1)];
%                 if i<5
%                     silence = zeros(sr*0.15,1);
%                     silTTL = [TTLson;zeros(length(silence)-length(TTLson),1)];
%                 else silence = zeros(length(TTLson),1);silTTL=TTLson;
%                 end
%                 s=[s;[so(:),so(:),sTTL,sTTL];[silence,silence,silTTL,silTTL]];
%             end
            %             length(s)
            Stim.Caracs.File = name_write;
            Stim.Caracs.Voice = voice;
            Stim.Caracs.Codes = sentec{ind_p(ind)};
            Stim.Caracs.stim_seg = trials_indices;
            Stim.Caracs.fsample = sr;
            Stim.Caracs.compression = tsmfactor(ind,voice);
            s= [svesola(:),svesola(:),sTTL(:),sTTL(:)];
            Stim.Mat = s;
            
            namf=[dirname,'stimsTrials',tech{:}];
            %             namf=[namf{:}];
            cd(namf)
            save(name_write,'Stim')
            cd('/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/sonsTrialsVsola/')
            wavwrite(s,sr,name_write);
        end
    end
end
%}
%% equalize intensity
ind_p=[1,2];
for ind = 1:length(ind_p)
    for voice=voices       %length(Codes)
        
        namf=[dirname,'stimsTrials',tech{:}];
        %             namf=[namf{:}];
        cd(namf)
        
        name_write= [voice,sentec{ind_p(ind)},'.wav'];
        name_write= [name_write{:}];
        
        %     load(name_write);
        
        name_load= [voice,sentec{ind_p(ind)},'.mat'];
        name_load= [name_load{:}];
        load(name_load)
        name_load2= [voice,sentec{ind_p(ind)},'.wav'];
        name_load2= [name_load2{:}];
        
        %     cd([dirname,'sonsTrialsVsola/'])
        %     s=wavread(name_load2);
        %     length(s)
        s = Stim.Mat;
        
        for techn=tech
            %         techn = tech{:};
            if strcmp(techn,'LKFS')
                factor = normalise_loudness(s(:,2), sr, lkfs);
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
            
            %                 cd([dirname,'stimstrials',techn])
            %                 save(name_load,'Stim')
            namf=[dirname,'stimsTrials',techn];
            namf=[namf{:}];
            cd(namf)
            save(name_load,'Stim')
            %         namf= [dirname,'sonsTrials',techn];
            %         namf=[namf{:}];
            %         cd(namf)
            %             wavwrite(Stim.Mat(:,1:4),sr,name_write);
            
            
        end
    end
end
   %{     
        for lang=langue
            
            try
                %             try
                %
                %                     cd([dirname,'stimsTest/'])
                %         name_load= [voice,sentec{ind},'.mat'];
                %         name_load= [name_load{:}];
                %         name_write= [voice,sentec{ind},'.wav'];
                %         name_write= [name_write{:}];
                %
                cd([dirname,'stimsTest/'])
                name_load= [voice,lang,sentec{ind},'.mat'];
                name_load= [name_load{:}];
                name_write= [voice,lang,sentec{ind},'.wav'];
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
                        Stim.Mat(:,1)=s(:,1)/res.STLmax/0.07;
                        Stim.Mat(:,2)=Stim.Mat(:,1);
                    end
                    
                    cd([dirname,'stimsTrials',techn])
                    save(name_load,'Stim')
                    
                    cd([dirname,'sonsTrials',techn])
                    wavwrite(s,sr,name_write);
                end
                
                cd([dirname,'stimsTrials',techn])
                save(name_load,'Stim')
                
%                 cd([dirname,'sonsTrials',techn])
%                 wavwrite(s,sr,name_write);
                
            end
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
            %             end
        end
        
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
                    Stim.Mat(:,1)=s(:,1)/res.STLmax*0.07;
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
        end
        %}


