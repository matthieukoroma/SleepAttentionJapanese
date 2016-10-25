function [Codes,Categories,NamesJ,NamesF,Freq,acrosslistpairs,attentionpairs] = load_stim_info(pathF) 

% fid = fopen('C:\\Data\\SleepAttentionJapanese\\Expe_Folder\\Expe_Material\\CorrespondancesTextest.csv');
% Scaned = textscan(fid,'%s%s%s%s%s');

cd(pathF)
load('trad_pairs.mat');

Codes = cell2mat(trad_pairs(:,1));%Scaned{1}(2:end);%reference number of the stimulus
Categories = cell2mat(trad_pairs(:,4));%Scaned{4}(2:end);%category of the word
NamesJ = trad_pairs(:,3);%Scaned{2}(2:end);%names in japanese
NamesF = trad_pairs(:,2);%Scaned{3}(2:end);%names in french
Freq = cell2mat(trad_pairs(:,5));%Scaned{5}(2:end);%freq of the french name
acrosslistpairs = cell2mat(trad_pairs(:,7));%Scaned{7}(2:end);
attentionpairs = cell2mat(trad_pairs(:,6));%Scaned{6}(2:end);

%guillaume
% fid = fopen('C:\\Data\\SleepAttentionJapanese\\Expe_Folder\\Expe_Material\\CorrespondancesTextest.csv');
% Scaned = textscan(fid,'%s%s%s%s%s');
% Codes = Scaned{1}(2:end);%reference number of the stimulus
% Categories = Scaned{4}(2:end);%category of the word
% NamesJ = Scaned{2}(2:end);%names in japanese
% NamesF = Scaned{3}(2:end);%names in french
% Freq = Scaned{6}(2:end);%freq of the french name
% acrosslistpairs = Scaned{8}(2:end);
% attentionpairs = Scaned{7}(2:end);