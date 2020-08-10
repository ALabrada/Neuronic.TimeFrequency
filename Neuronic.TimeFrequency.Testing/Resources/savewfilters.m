function savewfilters(wname)
[wtype,fname,family,bounds] =  ...
    wavemngr('fields',wname,'type','file','fn','bounds');
switch wtype
    case 1
        [Lo_R,Hi_R] = wfilters(wname,'r');
    case 2
        [~,~,~,~,~,~,Lo_R,Hi_R] = wfilters(wname);
end
csvwrite(['wfilters_', wname, '.csv'], [Lo_R;Hi_R]);
end