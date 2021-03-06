function evalwav(wname)
[wtype,fname,family,bounds] =  ...
    wavemngr('fields',wname,'type','file','fn','bounds');
switch wtype
    case 1
        [~,psi,x] = wavefun(wname, 10);
    case 2
        [~,psi,~,~,x] = wavefun(wname, 10);
    case 4
        [psi,x] = wavefun(wname, 10);
end
save(['wavefun_', wname, '.mat'], 'x', 'psi');
end