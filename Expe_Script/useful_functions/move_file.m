%%

directory = dir('*');
files_names = {directory.name};
files_names(1:2)=[];
for i=1:length(files_names)
%     newf=[files_names{i},'.wav'];
    movefile(files_names{i},[files_names{i},'.wav'])
end