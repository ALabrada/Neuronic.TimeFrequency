function savedwt(S,wname)
[a,d] = dwt(S,wname);
save(['dwt_', wname, '.mat'], 'S', 'a', 'd');
end