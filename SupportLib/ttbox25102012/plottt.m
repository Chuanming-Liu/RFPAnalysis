clear
clc

model=mkreadnd('prem.nd');
phase='P';
h=0;
dangle=0.5;
tt=mkttcurves(model,h,dangle);
handles=mkplotttcurves(tt,'b');