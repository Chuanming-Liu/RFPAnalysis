clear
clc

model=mkreadnd('ak135.nd');
dangle=mksmarttakeoff('P',model,0,2);
mkrayfan('P',model,0,dangle,'b');