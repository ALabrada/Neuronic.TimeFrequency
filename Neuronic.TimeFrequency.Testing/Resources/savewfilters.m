function savewfilters(wname)
[wtype,fname,family,bounds] =  ...
    wavemngr('fields',wname,'type','file','fn','bounds');
switch wtype
    case 1
        [Lo_R,Hi_R] = wfilters(wname,'r');
        F = [Lo_R;Hi_R];
        save(['wfilters_', wname, '.mat'], 'F');
    case 2
        [~,~,Lo_R2,Hi_R1,~,~,Lo_R1,Hi_R2] = wfilters(wname);
        F = [Lo_R1;Hi_R2;Lo_R2;Hi_R1];
        save(['wfilters_', wname, '.mat'], 'F');
end
end