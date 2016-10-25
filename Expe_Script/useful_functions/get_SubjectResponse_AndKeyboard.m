function [mouseKey timingResp keyboardKey timingKey]= get_SubjectResponse_AndKeyboard

Resp=[];
% [globalX, globalY, buttons]=Screen('GetMouseHelper', 3);
[keyDown, secs, keyboardCode] = KbCheck; timingKey=secs;
keyboardKey=find(keyboardCode);
if isempty(keyboardKey)
    keyboardKey=0;
else
    keyboardKey=keyboardKey(1);
end
% if length(find(buttons))==1
%     Resp=[find(buttons) GetSecs];
% end
% 
% % Find actual Resp key and timing
if ~isempty(Resp)
    mouseKey=Resp(1,1);
    timingResp = Resp(1,2);
else
    mouseKey=0;
    timingResp = NaN;
end



