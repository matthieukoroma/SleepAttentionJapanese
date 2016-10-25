%latin_sq

function latin_sq = generate_latin_square(lists)

MM = perms(1:3);%list
MM2 = perms(1:2);%testlist
MM3 = perms(1:2);%attended_ingored
%%
g=1;
lat_sqr_order=zeros(24,2);
sides = {'L','R'};
for i=1:length(MM)
    for j=1:length(MM2)
        for k =1:length(MM3)
            lat_sqr_order(g,1) = sides{MM3(k,1)};
            lat_sqr_order(g,2) = sides{MM3(k,2)};
            g=g+1;
        end
    end
end

