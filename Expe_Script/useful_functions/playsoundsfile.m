%%
sr=44100;
for i=1:288
    for voice = ['f','m']
        load([voice,num2str(i)])
        voice
        i
        sound(Stim.Mat(:,1:2),sr)
      
    end
end