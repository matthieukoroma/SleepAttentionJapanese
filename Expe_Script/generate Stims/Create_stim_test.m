% dirname = 'C:\\Data\\SleepAttentionAllocation\\Expe_Folder\\Expe_Material';
dirname = '/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/';

lettertypes = 'emrtu';
voices = {'f','m'};
jf={'f','j'};

[Codes,Categories,NamesJ,NamesF,Freq,accrosslistpairs,attentionpairs] = load_stim_info(dirname) ;

%%
load('/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/sons/TTLson');
dirsons = [dirname,'sonsTest'];
dirstim= [dirname,'stimsTest'];
Stim=[];

%%
for voice = voices
    %%
    for j_f=jf
        if strcmp(j_f,'j')
            dname = [dirname,'sons/',voice,'/stim_test_jap/'];
            Nam=NamesJ;
        elseif strcmp(j_f,'f')
            dname = [dirname,'sons/',voice,'/stim_test_fr/'];
            Nam=NamesF;
        end
        dname = [dname{:}];
        cd(dname)
        
        directory = dir('*');
        files_names = {directory.name};
        %         files_names(1:2)=[];
        
        for stimA =1:length(Codes)
            %%
            stim=Codes(stimA);
            
            stim_file = cellfun(@(x) strfind(x,[Nam{stimA} '.wav'])>0,files_names,'uni',0);
            index = find(cellfun(@isempty,stim_file)==0);
            [mm,minc] = min(cellfun(@length,files_names(index),'uni',1));
            index=index(minc);
            
            %                 stim_seg=zeros(length(list_carrying),2);
            %add silence 150ms : 0.15*1/sampling rate
            stim_sentence = [];
            cd(dname)
            fnam = files_names{index(1)};
            [s,sr] = wavread(fnam);
            sTTL = [TTLson;zeros(length(s)-length(TTLson),1)];
            silence = zeros(length(TTLson),1);silTTL=TTLson;
            stimT=[[s(:),s(:),sTTL,zeros(length(s),2)];[silence,silence,silTTL,silTTL]];
            
            cd(dirsons);
            name_write = [voice,j_f,num2str(Codes(stim))];%
            name_write = [name_write{:}];
            %                 wavwrite(stimT,sr,[name_write,'.wav']);
            
            cd(dirstim);
            Stim.Caracs.File = name_write;
            Stim.Caracs.Voice = voice{:};
            Stim.Caracs.Type = Categories(stim);
            %                 Stim.Caracs.order = perms_carrying(stim,:);
            Stim.Caracs.NamesF = NamesF{stim};
            Stim.Caracs.NamesJ = NamesJ{stim};
            Stim.Caracs.Code = Codes(stim);
            Stim.Caracs.stim_seg = length(s);
            Stim.Caracs.fsample = sr;
            Stim.Mat = stimT;
            save(name_write,'Stim')
        end
        %}
        
    end
end
