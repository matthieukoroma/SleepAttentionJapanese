%buildmatFinal

dirname = 'C:\\Data\\SleepAttentionAllocation\\Expe_Folder\\Expe_Material';
types = 'bcmw';

cat = {'Learning','Test'};
cattag = {'m','f'};
nb_syll = [2,3];
soundlen = cell(2,4);
SoundStruct = struct();


for k=1:2
    
    names = cell(20,4);
    len = cell(20,4);
    catpath = [dirname,'\\',char(cat(k))];
    for j=1:4
        
        if k==1
            lenmat = ones(20,2);
        else
            lenmat = ones(2,20);
        end

        directory = dir([catpath,'\\',char(cattag(k)),types(j),'*f.wav']);
        files = {directory.name};
        
        ite = 1;
        for i=1:numel(files)
            filename = char(files(i));
            names{i,j} = filename(2:4);
            [audio,fs] = wavread([catpath,'\\',filename]);
            audiolen = numel(audio)/fs;
            len{i,j} = audiolen;
            if k==1
                lenmat(i,1)=audiolen;
            else
                lenmat(2,i)=-audiolen;
            end
            clear('audio','fs');
            ite = ite+1;
        end
        soundlen{k,j}=lenmat;
        clear('filecar');
    end
    if k==1
        SoundStruct.Tale=names;
        SoundStruct.Talelen=len;
        
    else
        SoundStruct.Jab=names;
        SoundStruct.Jablen=len;
    end
end


%%
talelenmat = cell2mat(SoundStruct.Talelen);
jablenmat = cell2mat(SoundStruct.Jablen);
newTaleTimes = cell(20,4);
newJabTimes = cell(20,4);
newTale = cell(20,4);
newJab = cell(20,4);
Stories = cell(40,2);
tTimes = cell(40,2);
jTimes = cell(40,2);
for i=1:4
    talelist = SoundStruct.Tale(:,i);
    talelenlist = talelenmat(:,i);
    jablenlist = jablenmat(:,i);
    
    [ordtimes,I] = sort(talelenlist);
    
    newTaleTimes(:,i) = num2cell(ordtimes);
    pairTaleTimes = transpose(reshape(newTaleTimes(:,i),2,numel(newTaleTimes(:,i))/2));
    
    newJabTimes(:,i) = num2cell(jablenlist(I));
    pairJabTimes = transpose(reshape(newJabTimes(:,i),2,numel(newJabTimes(:,i))/2));
    
    newTale(:,i) = talelist(I);
    pairTales = transpose(reshape(newTale(:,i),2,numel(newTale(:,i))/2));
    
    Stories(10*(i-1)+1:10*i,:) = pairTales;
    tTimes(10*(i-1)+1:10*i,:) = pairTaleTimes;
    jTimes(10*(i-1)+1:10*i,:) = pairJabTimes;
end
Groupe={struct('Tale',Stories(:,1),'tLength',tTimes(:,1),'Jab',Stories(:,2),'jLength',jTimes(:,2)),struct('Tale',Stories(:,2),'tLength',tTimes(:,2),'Jab',Stories(:,1),'jLength',jTimes(:,1))};

save([dirname,'\\TaleCorrespondances.mat'],'Groupe');



    