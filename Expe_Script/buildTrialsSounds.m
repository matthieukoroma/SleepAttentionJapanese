% build_Trial
dirname = '/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/';
dirtest =  '/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/';
dirsons = '/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/Ex';
[Codes,Categories,NamesJ,NameF,Freq,acrosslistpairs,attentionpairs] = load_stim_info(dirname);

voices = {'f','m'};
cues = {'suivez_la_femme','suivez_lhomme'};
sides = {'L','R'};
sr = 44100;
equa_tech = {'RMS'};%,'STL'};%{'LFKS'};%
% pairsa = unique(attentionpairs);
pairsa=unique(attentionpairs);
pairst = unique(acrosslistpairs);
load('/home/lscp00/Matthieu/SleepAttentionJapanese/Expe_Material/sons/TTLson');

%%
for stim=1:1;%length(pairsa)
    
    %     try
    stim
    for tech=equa_tech
        
        ind_p = find(attentionpairs==attentionpairs(stim));
        side=voices(randperm(length(voices)));
        stims=cell(length(ind_p),length(voices));
        %cue = cues{randi(2)};
        silence_S=zeros(0.5*sr,4);
        silence_cue = cell(1,2);
        techn=[tech{:}];
        
        
        %%
        %         for cue=1:length(cues)
        cue=1;
        for voice=1:length(voices)
            for ind = 1:length(ind_p)
                namewave = [dirname,'stimsTrials',tech,'/',voices{voice},num2str(Codes(ind)),'.mat'];
%                 [s{ind,voice},sr] = wavread(namewave);
                cdname = [dirname,'stimsTrials',tech];
                cdname = [cdname{:}];
                cd(cdname)
                name_load = [voices(voice),num2str(Codes(ind_p(ind)))];%
                name_load = [name_load{:}];
                load(name_load)
                stims{ind,voice}=[Stim.Mat(:,1),Stim.Mat(:,4)];
            end
            
            %             dnam = [voices(voice),cues{cue}];
            %             dnam = [dnam{:}];
            %             load(dnam)
            %             stimCue=Stim.Mat;
            %             silence_cue{voice}(1:length(stimCue),:) = stimCue(:,[1,3]);
        end
        %         len = diff(cellfun(@length,silence_cue));
        %         if len>0
        %             silence_cue{2}=silence_cue{2}(1:length(silence_cue{1}),:);
        %         elseif len<0
        %             silence_cue{1}=silence_cue{1}(1:length(silence_cue{2}),:);
        %         end
        
        maxl=max(max([cellfun(@length,stims,'uni',1)]));
        stimS = zeros(maxl,2);
        stimS(1:length(stims{1,1}),[1,3]) = stims{1,1};
        stimS(1:length(stims{2,2}),[2,4]) = stims{2,2};
        silenceS = [];
%         silenceS(:,[1,3,2,4]) = [silence_cue{1},silence_cue{2}];
        
        %%
        if cue ==1
            side = sides{1};
        else
            side = sides{2};
        end
        Stim.Mat = [silenceS;stimS];
        
        
        cd(dirsons);
        name_write = [side,num2str(Codes(ind_p(cue))),tech];%
        %             name_write = ['trial',tech,'.wav'];%
        name_write = [name_write{:}];
        wavwrite([stimS],sr,name_write);
        Stim.Mat = [stimS];
        cd(dirtest)
        %
                save(name_write,'Stim')
        
        %%
        if cue ==1
            side = sides{2};
        else
            side = sides{1};
        end
        
        maxl=max(max([cellfun(@length,stims,'uni',1)]));
        stimS = zeros(maxl,2);
        stimS(1:length(stims{1,1}),[2,4]) = stims{1,1};
        stimS(1:length(stims{2,2}),[1,3]) = stims{2,2};
        %         silenceS(:,[1,3,2,4]) = [silence_cue{2},silence_cue{1}];
        cd(dirsons);
        name_write = [side,num2str(Codes(ind_p(cue))),tech];%
        %             name_write = ['trial',tech,'.wav'];%
        name_write = [name_write{:}];
                Stim.Mat = [stimS];
        %
         
        wavwrite(stimS,sr,name_write);
        cd(dirtest)
        %
                save(name_write,'Stim')
        %%
        if cue ==1
            side = sides{1};
            stimL=ind_p(2);
        else
            side = sides{1};
            stimL=ind_p(1);
        end
        
        %         silenceS(:,[1,3,2,4]) = [silence_cue{2},silence_cue{1}];
        maxl=max(max([cellfun(@length,stims,'uni',1)]));
        stimS = zeros(maxl,2);
        stimS(1:length(stims{1,2}),[1,3]) = stims{1,2};
        stimS(1:length(stims{2,1}),[2,4]) = stims{2,1};
        
        cd(dirsons);
        name_write = [side,num2str(Codes(stimL)),tech];%
        %             name_write = ['trial',tech,'.wav'];%
        name_write = [name_write{:}];
        % %         Stim.Mat = [silenceS;stimS];
        Stim.Mat = [stimS];
        
        wavwrite([stimS],sr,name_write);
        cd(dirtest)
        %
                save(name_write,'Stim')
        %
        if cue ==1
            side = sides{1};
            stimL=ind_p(2);
        else
            side = sides{2};
            stimL=ind_p(1);
        end
%         silenceS(:,[1,3,2,4]) = [silence_cue{1},silence_cue{2}];
        stimS(1:length(stims{1,2}),[2,4]) = stims{1,2};
        stimS(1:length(stims{2,1}),[1,3]) = stims{2,1};
        
        cd(dirsons);
        name_write = [side,num2str(Codes(stimL)),tech];%
        %             name_write = [cue,'trial2',tech,'.wav'];%
        name_write = [name_write{:}];
        Stim.Mat = [stimS];
         cd(dirtest)
        %
        save(name_write,'Stim')
        wavwrite(stimS,sr,name_write);
        
    end
end
        
        %{
        %-----------------------------
        
        %% test_trials
        
        %     try
        cdname = [dirname,'stimsTest',techn];
        %         cdname = [cdname{:}];
        cd(cdname)
        
        
        if acrosslistpairs(stim)==0
            continue
        end
        
        
        ind_p = find(acrosslistpairs==acrosslistpairs(stim));
        
        
        for tg=1:length(ind_p)
            target = ind_p(tg);
            %         side=voices(randperm(length(voices)));
            %side=side(1);
            for side = voices
                for ord=1:2
                    %%
                    if ord==1
                        order = ind_p([1,2]);
                    else
                        order=ind_p([2,1]);
                    end
                    stimSIL=[];
                    
                    silence_cue=zeros(0.5*sr,2);
                    
                    [beep,samplingRate] = MakeBeepAP(450,0.1,sr);
                    TTLs= zeros(length(beep),1);
                    TTLs(1:length(TTLson))=TTLson;
                    stimSIL=[repmat(beep(:),1,2),repmat(TTLs,1,2)];
                    for ind = 1:length(ind_p)
                        cdname = [dirname,'stimsTest',techn];
                        %         cdname = [cdname{:}];
                        cd(cdname)
                        %             namewave = [dirname,'stimsTrials/',voices{voice},num2str(Codes(ind)),'.wav'];
                        %             [s{ind,voice},sr] = wavread(namewave);
                        name_load = [side,'j',num2str(Codes(target)),'.mat'];%
                        name_load = [name_load{:}];
                        load(name_load)
                        stimSIL=[stimSIL;zeros(0.4*sr,4);Stim.Mat];
                        
                        name_load = [side,'veut_dire1','.mat'];%
                        name_load = [name_load{:}];
                        load(name_load)
                        stimSIL = [stimSIL;zeros(0.15*sr,4);Stim.Mat];
                        
                        name_load = [side,'f',num2str(Codes(order(ind))),'.mat'];%
                        name_load = [name_load{:}];
                        load(name_load)
                        stimSIL = [stimSIL;zeros(0.5*sr,4);Stim.Mat];
                    end
                    
                    [beep2,samplingRate2] = MakeBeepAP(530,0.1,sr);
                    stimSIL=[stimSIL;zeros(0.5*sr,4);[repmat(beep2(:),1,2),repmat(TTLs,1,2)]];
                    
                    cd(dirsons);
                    name_write = ['a',side,'t',num2str(Codes(target)),'o',num2str(Codes(order(1))),num2str(Codes(order(2))),tech];%
                    name_write = [name_write{:}];
                    cdname = [dirname,'testStims2'];
                    %                     cdname = [cdname{:}];
                    cd(cdname)
                    Stim.Mat = stimSIL;
                    save(name_write,'Stim')
                    %                 name_write = [side,'test',tech,'.wav'];%
                    %                     name_write = [name_write{:}];
                    %                     wavwrite(stimSIL,sr,name_write);
                end
            end
            
        end
        %     end
        %end
    end
end

        %}
        
%% within list
% cdname = [dirname,'stimsTest',tech];
% % cdname = [cdname{:}];
% cd(cdname)


%         if acrosslistpairs(stim)==0
%             continue
%         end

%         ind_p = find(acrosslistpairs==acrosslistpairs(stim));
%
for i=5:size(list_sav,2) %list_test
    
    dyads = nchoosek(list_sav(:,i),2);
    for in=1:size(dyads,1)
        ind_p = dyads(in,:);
        
        for tech=equa_tech
            for tg=1:size(ind_p)
                target = ind_p(tg);
                %         side=voices(randperm(length(voices)));
                %side=side(1);
                for side = voices
                    for ord=1:2
                        if ord==1
                            order = ind_p([1,2]);
                        else
                            order=ind_p([2,1]);
                        end
                        stimSIL=[];
                        
                        silence_cue=zeros(0.5*sr,2);
                        
                        [beep,samplingRate] = MakeBeepAP(450,0.1,sr);
                        TTLs= zeros(length(beep),1);
                        TTLs(1:length(TTLson))=TTLson;
                        stimSIL=[repmat(beep(:),1,2),TTLs,TTLs];
                        for ind = 1:length(ind_p)
                            %             namewave = [dirname,'stimsTrials/',voices{voice},num2str(Codes(ind)),'.wav'];
                            %             [s{ind,voice},sr] = wavread(namewave);
                            cdname = [dirname,'stimsTest',tech];
                                    cdname = [cdname{:}];
                            cd(cdname)
                            name_load = [side,'j',num2str(Codes(target)),'.mat'];%
                            name_load = [name_load{:}];
                            load(name_load)
                            stimSIL=[stimSIL;zeros(0.4*sr,4);Stim.Mat];
                            
                            name_load = [side,'veut_dire1','.mat'];%
                            name_load = [name_load{:}];
                            load(name_load)
                            stimSIL = [stimSIL;zeros(0.15*sr,4);Stim.Mat];
                            
                            name_load = [side,'f',num2str(Codes(order(ind))),'.mat'];%
                            name_load = [name_load{:}];
                            load(name_load)
                            stimSIL = [stimSIL;zeros(0.5*sr,4);Stim.Mat];
                        end
                        
                        [beep2,samplingRate2] = MakeBeepAP(530,0.1,sr);
                        stimSIL=[stimSIL;zeros(0.5*sr,4);[repmat(beep2(:),1,2),repmat(TTLs,1,2)]];
                        
                       
                        name_write = ['w',side,'t',num2str(Codes(target)),'o',num2str(Codes(order(1))),num2str(Codes(order(2))),tech];%
                        %                 name_write = [side,'test',tech,'.wav'];%
                        name_write = [name_write{:}];
                        cdname = [dirname,'testStims2'];
                        %                     cdname = [cdname{:}];
                        cd(cdname)
                        Stim.Mat = stimSIL;
                        save(name_write,'Stim')
                         cd(dirsons);
                         wavwrite(stimSIL,sr,name_write);
                    end
                end
                
            end
        end
    end
end
%}