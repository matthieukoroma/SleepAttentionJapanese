%latin_sq

function lat_sq = generate_latin_square(list_sav)

MM = perms(1:3);%list
MM2 = perms(1:2);%testlist
MM3 = perms(1:2);%attended_ingored
%%
g=1;
lat_sqr=([]);

for k=1:length(MM3)
    for j =1:length(MM2)
        for i=1:length(MM)
            lat_sqr(g).wake_attended = list_sav(:,(MM(i,1)-1)*4+(MM2(j,1)-1)*2+MM3(k,1));
            lat_sqr(g).wake_ignored = list_sav(:,(MM(i,1)-1)*4+(MM2(j,1)-1)*2+MM3(k,2));
            lat_sqr(g).N2_attended = list_sav(:,(MM(i,1)-1)*4+(MM2(j,1)-1)*2+MM3(k,1));
            lat_sqr(g).N2_ignored = list_sav(:,(MM(i,1)-1)*4+(MM2(j,1)-1)*2+MM3(k,2));
            lat_sqr(g).N3_attended = list_sav(:,(MM(i,1)-1)*4+(MM2(j,1)-1)*2+MM3(k,1));
            lat_sqr(g).N3_ignored = list_sav(:,(MM(i,1)-1)*4+(MM2(j,1)-1)*2+MM3(k,2));
            g=g+1;
        end
    end
end

