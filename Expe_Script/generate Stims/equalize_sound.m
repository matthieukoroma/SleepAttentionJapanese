%equalize voices length across attended pairs with vsola
%equalize in intensity
dirname = '/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/';
[Codes,Categories,NamesJ,NamesF,Freq,accrosslistpairs,attentionpairs] = load_stim_info(dirname);

voices = {'f','m'};
sr =  44100;% Stim.Caracs.fsample;
techn ={'LKFS'};%'RMS'};%'STL';% {'LFKS'};%
carrying_sentences = {'suivez_la_femme','suivez_lhomme','veut_dire'};
lkfs = -30;

pairs = unique(attentionpairs);
P = [round(44100*8/1000),round(44100*16/1000)];%the longest likely pitch period
s = cell(2,length(Codes));
load('/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/sons/TTLson');

% pairs=unique(attentionpairs);
%%
for tech=techn
    for stim=114;%1:length(pairs)
        
        %get pairs
        ind_p = find(attentionpairs==pairs(stim));
        
        s = cell(length(ind_p),length(voices));
        len = zeros(length(ind_p),length(voices));
        
        %%
        
        for ind = 1:length(ind_p)
            for voice=1:length(voices)
                %             namewave = [dirname,'stimsTrials/',voices{voice},num2str(Codes(ind)),'.wav'];
                %             [s{ind,voice},sr] = wavread(namewave);
                nam = ['/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/stimsTrials', tech{:}]
                cd(nam)
                name_write = [voices(voice),num2str(Codes(ind_p(ind)))];%
                name_write = [name_write{:}];
                load([name_write,'.mat'])
                bounds{ind,voice} = max(Stim.Caracs.stim_seg(:))
                len(ind,voice)=bounds{ind,voice};
            end
        end
        
        %     lenmean=mean(len(:));
        %     tsmfactor = lenmean./len;
        
        lenmin = (mean(len(:)));
        tsmfactor = lenmin./len;
        %     tsmfactor1 = mean([len(1,1),len(2,2)]);
        %     tsmfactor2 = mean([len(2,1),len(1,2)]);
        %     tsmfactor = zeros(2);
        %
        %     tsmfactor1(1,1) = tsmfactor1;
        %     tsmfactor1(2,2) = tsmfactor1;
        %     tsmfactor2(2,1) = tsmfactor2;
        %     tsmfactor2(1,2) = tsmfactor2;
        
        
        %apply vsola
        %%
        for voice=1:length(voices)
            for ind = 1:length(ind_p)
                name_write = [voices(voice),num2str(Codes(ind_p(ind)))];%
                name_write = [name_write{:}];
                nam = ['/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/stimsTrials', tech{:}]
                
                cd(nam)
                 clear Stim
                load([name_write,'.mat'])
                stim_m = Stim.Mat(:,1);
                %              stimF =[];
                %             for i=1:size(bounds{ind,voice},1)
                %                 stimF = [stimF;stim_m(bounds{ind,voice}(i,1):min([length(stim_m),bounds{ind,voice}(i,2)]),1)];
                %             end
                if tsmfactor(ind,voice)==1
                    s=stim_m;
                else
                    s = vsola(stim_m, tsmfactor(ind,voice), P(voice));
                end
                %             length(svesola)
                %             s=[];stim_seg=[];s=[];
                %insert silence
                %             trials_indices=zeros(size(bounds{ind,voice},1),2);
                %             for i=1:size(bounds{ind,voice},1)
                %                 trials_indices(i,:) = [max([floor(bounds{ind,voice}(i,1)*tsmfactor(ind,voice)),1]),min([length(svesola),ceil(bounds{ind,voice}(i,2)*tsmfactor(ind,voice))])];
                %                 so=svesola(trials_indices(i,1):trials_indices(i,2));
                %
                %                 sTTL = [TTLson;zeros(length(so)-length(TTLson),1)];
                %                 if i<5
                %                     silence = zeros(sr*0.3,1);
                %                     silTTL = [TTLson;zeros(length(silence)-length(TTLson),1)];
                %                 else silence = zeros(length(TTLson),1);silTTL=TTLson;
                %                 end
                %                 s=[s;[so(:),so(:),sTTL,zeros(length(so),1)];[silence,silence,silTTL,silence]];
                %             end
                %             length(s)
                
                sTTL = [TTLson;zeros(length(s)-length(TTLson),1)];silence = zeros(length(TTLson),1);
                
                Stim.Mat = [[s(:),s(:),sTTL,sTTL];[silence,silence,TTLson,TTLson]];
                Stim.Caracs.File = name_write;
                Stim.Caracs.Voice = voice;
                Stim.Caracs.Type = name_write;
                Stim.Caracs.order = perms_carrying(ind_p(ind),:);
                Stim.Caracs.NamesF = NamesF{ind_p(ind)};
                Stim.Caracs.NamesJ = NamesJ{ind_p(ind)};
                Stim.Caracs.Codes = Codes(ind_p(ind));
                Stim.Caracs.stim_seg = length(s);%trials_indices;
                Stim.Caracs.fsample = sr;
                Stim.Caracs.compression = tsmfactor(ind,voice);
                %             length(s)
                %             Stim.Mat = s;
                nam = ['/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/stimsTrialsVsola', tech{:}]
                 cd(nam)
                 save(name_write,'Stim')
%                 %             cd('/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/sonsTrialsVsola/')
%                 wavwrite(s,sr,name_write);
            end
        end
    end
end
