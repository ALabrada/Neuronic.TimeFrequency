function savedwt(S,wname)
[a,d] = dwt(S,wname);
dlmwrite(['dwt_', wname, '.csv'], S);
dlmwrite(['dwt_', wname, '.csv'], [a;d], '-append');
end