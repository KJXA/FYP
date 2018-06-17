function windingcount = phasefix_windingcount(array)
winding=0;
for i=1:(length(array)-1)
    loopCounter=0;
    while abs(array(i)-array(i+1)) > 1.9*pi
        if array(i)-array(i+1) > 0
            array(i+1:end)=array(i+1:end)+2*pi;
            winding = winding +1;
            if loopCounter>=30
                break;
            end
            loopCounter = loopCounter +1;
        elseif array(i)-array(i+1) < 0
            array(i+1:end)=array(i+1:end)-2*pi;
            winding = winding -1;
            if loopCounter>=30
                break;
            end
            loopCounter = loopCounter +1;
        end
    end
end
windingcount=winding;