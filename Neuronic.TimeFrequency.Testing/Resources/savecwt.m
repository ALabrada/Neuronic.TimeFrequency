function savecwt(S,wname)
scales = 1:16;
c = cwt(S,scales,wname);
csvwrite(['cwt_', wname, '.csv'], [S; c]);
end