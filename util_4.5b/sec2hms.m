function hms = sec2hms(tox);

sec=mod(tox,60);
tox=(tox-sec)/60;
mn=mod(tox,60);
hr=(tox-mn)/60;

hms = [int2str(hr),'h: ',int2str(mn),'m : ',int2str(sec), 's'];