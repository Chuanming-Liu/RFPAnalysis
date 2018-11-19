clear
clc

clr=mkreadclr('iasp91.clr','silent');
model=mkclr2model(clr,10,'spherical');
mkwritend(model,'iasp91.nd');
