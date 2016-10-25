% dirname = 'C:\\Data\\SleepAttentionAllocation\\Expe_Folder\\Expe_Material';
dirname = '/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/';

lettertypes = 'emrtu';
voices = {'f','m'};
sr=44100;
[Codes,Categories,NamesJ,NamesF,Freq,accrosslistpairs,attentionpairs] = load_stim_info(dirname) ;
sr =  44100;% Stim.Caracs.fsample;
techn ={'LKFS','RMS'};%
carrying_sentences = {'suivez_la_femme','suivez_lhomme','veut_dire'};
lkfs = -30;
%%

% files_names
%%
% %generate a random matrix of possible ordering of the sentences
% %ordered as they appear in files_names
list_carrying = {'mot','traduit','traduction','utilise','équivalent',};
cd(dirname)
load('perms_carrying')
% %
% perms_carrying = perms(1:length(list_carrying));
% perms_carrying = repmat(perms_carrying,ceil(length(Codes)/length(perms_carrying)),1);
% perms_carrying = perms_carrying(randperm(length(perms_carrying)),:);
% %
% length_stims = zeros(length(Codes),2);

%%
load('/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/sons/TTLson');
dirsons = [dirname,'sonsTrials'];
dirstim= [dirname,'stimsTrials'];
Stim=[];

%%

for stimA = 1:length(Codes);%73;%[30,51,54,72,102,189,228,232,278];
    stim=Codes(stimA)
    
    for voice = voices
        name_directory = [dirname,'sons/',voice, '/sentences/'];
        name_directory = [ name_directory{:}];
        %         if strcmp(NamesF{stim},'bouche')
        %             keyboard
        %         end
        
        cd(name_directory)
        directory = dir('*');
        files_names = {directory.name};
        
        stim_file = cellfun(@(x) strfind(x,NamesF{stim})>0,files_names,'uni',0);
        index1 = find(cellfun(@isempty,stim_file)==0);
        
        stim_file2 = cellfun(@(x) strfind(x,NamesJ{stim})>0,files_names,'uni',0);
        index2 = find(cellfun(@isempty,stim_file2)==0);
        indexes = index1(ismember(index1(:),index2(:),'rows'));
        indexes_M=zeros(1,length(indexes));
        
        if length(indexes)>5
            keyboard
        end
        %get the order of files in indexes like in list_carrying
        for i=1:length(list_carrying)
            place_list_carrying = cellfun(@(x) strfind(files_names{indexes(i)},[x,'_'])>0,list_carrying,'uni',0);
            indexes_M(i)=find(cellfun(@isempty,place_list_carrying)==0);
        end
        [a,order_list]=sort(indexes_M);
        indexes = indexes(order_list);
        
        perms_carrying_s = perms_carrying(stim,:);
        redo=1;
        %avoid "on traduit" and "la traduction" at the same position
        while redo
            paircode = find(attentionpairs==attentionpairs(stim));
            perms_carrying_pair = perms_carrying(paircode,:);
            perms_carrying_pair_test=perms_carrying_pair;
            
            %2 and 3 indexes in list_carrying of
            perms_carrying_pair_test(perms_carrying_pair_test==2)=3;
            if prod(diff(perms_carrying_pair_test))==0
                %switch so as to renew the possibility of finding a good
                %association
                indexToPermWith = stim+randi(size(perms_carrying,1)-stim,1);
                perms_carrying(paircode(2),:) = perms_carrying(indexToPermWith,:);
                perms_carrying(indexToPermWith,:) = perms_carrying_pair(2,:);
                redo=1;
            else redo=0;
            end
        end
        list_carrying_s = list_carrying(perms_carrying(stim,:));
        %     carrying_s1 = stim_file(strfind(@(x) strfind(x,list_carrying_s{1}),files_names(stim_file)));
        
        %%
        
        
        stim_seg=zeros(length(list_carrying),2);
        %add silence 150ms : 0.15*1/sampling rate
        stim_sentence = [];
        %         TTL = [];
        dname = [dirname,'sons/',voice,'/sentences/'];
        dname = [dname{:}];
        cd(dname)
        
        for tech=techn
            for i=1:length(list_carrying)
                
                fnam = files_names{indexes(perms_carrying(stim,i))};
                
                [s,sr] = wavread(fnam);
                %         wavwrite(stim_sentence,sr,[name_write,'.wav']);
                
                if strcmp(tech,'LKFS')
                    
                    factor = normalise_loudness(s(:,1)', sr, lkfs);
                    Stim.Caracs.LKFSfactor = factor;
                    Stim.Caracs.LKFSval = lkfs;
                    s=s*factor;
                elseif strcmp(tech,'RMS')
                    
                    s=s(:,1)./rms(s(:,1))*0.07;
                    
                    %             Stim.Mat(:,1:2)=Stim.Mat(:,1:2);
                elseif strcmp(techn,'STL')
                    %get info for sonic normalization
                    %             res = Loudness_TimeVaryingSound_Moore(s(:,1), sr, 'mic',0);
                    %             Stim.Caracs.STLmax = res.STLmax;
                    Stim.Mat(:,1)=s(:,1)/STLmax;
                    Stim.Mat(:,2)=Stim.Mat(:,1);
                end
                
                
                %         cd(dirsons)
                %         wavwrite(Stim.Mat(:,1:2),sr,name_write)
            end
        end
        sTTL = [TTLson;zeros(length(s)-length(TTLson),1)];
        if i<5
            silence = zeros(sr*0.3,1);
            silTTL = [TTLson;zeros(length(silence)-length(TTLson),1)];
        else silence = zeros(length(TTLson),1);silTTL=TTLson;
        end
        stim_seg(i,:) = [length(stim_sentence)+1,length(stim_sentence)+length(s)+1];
        stim_sentence=[stim_sentence;[s(:),s(:),sTTL,zeros(length(s),1)];[silence,silence,silTTL,silence]];
    end
    %                 name_write = [dirname,'/SonsTrial/',];
    name_write = [voice,num2str(Codes(stim))];%
    name_write = [name_write{:}];
    Stim.Caracs.File = name_write;
    Stim.Caracs.Voice = voice{:};
    Stim.Caracs.Type = Categories(stim);
    Stim.Caracs.order = perms_carrying(stim,:);
    Stim.Caracs.NamesF = NamesF{stim};
    Stim.Caracs.NamesJ = NamesJ{stim};
    Stim.Caracs.Code = Codes(stim);
    Stim.Caracs.stim_seg = stim_seg;
    Stim.Caracs.fsample = sr;
         
    %             cd(dirsons);
    cd([dirstim,tech{:}]);

       Stim.Mat = stim_sentence;
    save(name_write,'Stim')
end
cd(dirname)
save('perms_carrying','perms_carrying')
