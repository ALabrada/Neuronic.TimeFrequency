function saveidwt(S,wname)
[a,d] = dwt(S,wname);
S2 = idwt(a,d,wname);
save(['idwt_', wname, '.mat'], 'S', 'S2');
end