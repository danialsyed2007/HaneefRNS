%Read RNS data

[dex,dat]=ReadECoGData('D:\Haneefproject\PAtientA\Baylor College of Medicine, Houston_WL_11969371 EXTERNAL #PHI\Baylor College of Medicine, Houston_WL_11969371 Data EXTERNAL #PHI\133349710605950000.dat', 'test.lay');
dat=dat{1};
dat=dat';
[xdn,~,~,spike]=dbtDenoise(dat, 250, 0.1);


[b,a]=butter(4,[4 8]/125,'bandpass')
filter(b,a,xdn)
josh=filter(b,a,xdn)