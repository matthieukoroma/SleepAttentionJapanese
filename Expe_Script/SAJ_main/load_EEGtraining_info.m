function [Codes,Types,Excerpt,Voice] = load_EEGtraining_info(stimPath)

cd(stimPath)
load('EEGtraining_info');

Codes = cell2mat(EEGtraining_info(:,1));%Scaned{1}(2:end);%reference number of the stimulus
Types = EEGtraining_info(:,2);%Scaned{4}(2:end);%category of the stim
Excerpt = cell2mat(EEGtraining_info(:,3));%Scaned{2}(2:end);%excerpt 1 or 2
Voice = cell2mat(EEGtraining_info(:,4));%Scaned{3}(2:end);%voice m or f
