acrosspairs=[];
for i=1:4:12
for j=1:16
acrosspairs(list_sav(j,i)) = pair;
acrosspairs(list_sav(j,i+2)) = pair;
pair=pair+1;
end
end
for i=2:4:12
for j=1:16
acrosspairs(list_sav(j,i)) = pair;
acrosspairs(list_sav(j,i+2)) = pair;
pair=pair+1;
end
end
acrosspairs;
jo=mat2cell(acrosspairs,1,repmat(1,1,288));
trad_pairs(:,7)=jo';


for i=1:1
for j=1:16
acrosspairs(list_test(j,i)) = pair;
acrosspairs(list_test(j,i+2)) = pair;
pair=pair+1;
end
end
for i=2:2
for j=1:16
acrosspairs(list_test(j,i)) = pair;
acrosspairs(list_test(j,i+2)) = pair;
pair=pair+1;
end
end
acrosspairs;
jo=mat2cell(acrosspairs,1,repmat(1,1,288));
trad_pairs(:,7)=jo';