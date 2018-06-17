function fixed_array = phasefix(array)
for i=1:(length(array)-1)
    while abs(array(i)-array(i+1)) > 1.9*pi
        if array(i)-array(i+1) > 0
            array(i+1:end)=array(i+1:end)+2*pi;
        elseif array(i)-array(i+1) < 0
            array(i+1:end)=array(i+1:end)-2*pi;
        end
    end
end
fixed_array=array;