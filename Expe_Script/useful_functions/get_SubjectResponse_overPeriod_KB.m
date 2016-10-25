function [mouseKey timingResp keyboardKey]= get_SubjectResponse_overPeriod(timeStart,timeMax,breakFlag)

if nargin<3
    breakFlag=0;
end

keyboardKey=0;
% Check buttons are not yet pressed
% [globalX, globalY, buttons]=Screen('GetMouseHelper', 3);
% while any(buttons) % if already down, wait for release
%     disp('Subject should release button')
%     [globalX, globalY, buttons]=Screen('GetMouseHelper', 3);
% end
buttons=[];

% Let Subject timeMax seconds to reply
Resp=[];
fprintf('\n');

progressi = sprintf('you have ... %3.1f s left',round((timeStart+timeMax-GetSecs)*10)/10);
fprintf(progressi);
while GetSecs<timeStart+timeMax % wait for press
    pause(0.1)
    progress = [repmat('\b',1,length(progressi)-1),sprintf('you have ... %3.1f s left',round((timeStart+timeMax-GetSecs)*10)/10)];
    fprintf(progress);
%     [globalX, globalY, buttons]=Screen('GetMouseHelper', 3);
    if keyboardKey==0
        [keyboardKey, secs, keyCode] = KbCheck;
    end
    if length(find(keyCode))==1
        Resp=[Resp ; [find(keyCode) GetSecs]];
        break;
    end
end
pause(0.1)
while breakFlag==0 && GetSecs<timeStart+timeMax
    pause(0.1)
    progress = [repmat('\b',1,length(progressi)-1),sprintf('you have ... %3.1f s left',round((timeStart+timeMax-GetSecs)*10)/10)];
    fprintf(progress);
    if keyboardKey==0
        [keyboardKey, secs, keyCode] = KbCheck;
    end
end
fprintf('\n');
% Find actual Resp key and timing
if ~isempty(Resp)
    mouseKey=Resp(1,1);
    timingResp = Resp(1,2);
else
    mouseKey=0;
    timingResp = NaN;
end


