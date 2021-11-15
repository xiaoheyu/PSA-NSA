filename indata 'H:\SA\data.csv';
filename sparseij 'H:\SA\Rook.csv';
OPTIONS LINESIZE=72;
TITLE 'PSA-NSA mixtures: SAR-SMA';

data step1; infile indata dlm=','; input EV1 PopDen2 tr_PopDen2 tr_Respira tr_Medin new_id; 
if new_id='.' then delete;
y=tr_Respira;
sid=new_id;
run;
proc standard data=step1 out=step1(replace=yes) mean=0; var y; run;
proc means noprint; var y; output out=outn n=n; run;
data outn(replace=yes); set outn(keep=n); run;
proc univariate data=step1 normal; var y; run;

data swmo; infile sparseij dlm=','; input NEW_ID NID WEIGHT;
if new_id='.' then delete;
sid=new_id;
weight=1;
drop NEW_ID;
run;
data swm; set swmo(keep=sid); set swmo(keep=nid); set swmo(keep=weight); run;

proc spatialreg data=step1 Wmat=swm corrb; model y = /type=SAR; spatialid sid; run;
ods trace on;
proc spatialreg data=step1 Wmat=swm corrb; model y = /type=SARMA; spatialid sid; output out=outsar pred=yhat resid=yr; 
ods output parameterestimates=outparms;
run;
ods trace off;
proc transpose data=outparms out=outsa prefix=rho; var estimate; run;
data outsa(replace=yes); set outsa; rho=rho2; theta=rho3; drop rho1-rho4 _name_; run;
proc print data=outsa; run;

data outsar(replace=yes); set outsar(keep=y); set outsar(keep=yhat); run;
data sparseij; set swm(keep=sid nid); run;
PROC IML;
use outsar; read all into yin;
yo=yin[,1]; mu=yin[,2];
y=yo;
use outn; read all into n;
use outsa; read all into saparms;
rho=saparms[1];
theta=saparms[2];
use sparseij; read all into ij;

indx =n*(ij[,1]-1) + ij[,2];
iden=i(n);
c=shape(0,n,n);
c[indx]=1;

IDEN=I(n);
ONE=J(n,1);
sumc=one`*c*one;
sumc2= one`*c*c*one;
print sumc sumc2;

y=(inv(i(n) - theta*diag(1/(one`*c))*c))*((i(n) - rho*diag(1/(one`*c))*c)*yo - mu);

mc=n/(one`*c*one)*(y-(one`*y/n)*one)`*c*(y-(one`*y/n)*one)/((y-(one`*y/n)*one)`*(y-(one`*y/n)*one));
emc=-1/(n-1);
semc=sqrt(2/(one`*c*one));
zmc=(mc-emc)/semc;
pr=1-probnorm(abs(zmc));
print mc emc semc zmc pr;

gro=((n-1)/(2*one`*c*one))*(one`*diag(c*( (y*one`-one*y`)#(y*one`-one*y`))) *one)/((y-(one`*y/n)*one)`*(y-(one`*y/n)*one));
segr=sqrt(2/(one`*c*one) + 2*(one`*c*c*one)/(one`*c*one)**2);
z=(gro-1)/segr;
pr=1-probnorm(abs(z));

gr=((n-1)/(one`*c*one))*(diag((y-(one`*y/n)*one)*(y-(one`*y/n)*one)`)*one)`*c*one/((y-(one`*y/n)*one)`*(y-(one`*y/n)*one)) - ((n-1)/n)*mc;
print gro segr z pr;

y=(i(n) - rho*diag(1/(one`*c))*c)*yo - mu;
mc=n/(one`*c*one)*(y-(one`*y/n)*one)`*c*(y-(one`*y/n)*one)/((y-(one`*y/n)*one)`*(y-(one`*y/n)*one));
emc=-1/(n-1);
semc=sqrt(2/(one`*c*one));
zmc=(mc-emc)/semc;
pr=1-probnorm(abs(zmc));
print mc emc semc zmc pr;

gro=((n-1)/(2*one`*c*one))*(one`*diag(c*( (y*one`-one*y`)#(y*one`-one*y`))) *one)/((y-(one`*y/n)*one)`*(y-(one`*y/n)*one));
segr=sqrt(2/(one`*c*one) + 2*(one`*c*c*one)/(one`*c*one)**2);
z=(gro-1)/segr;
pr=1-probnorm(abs(z));

gr=((n-1)/(one`*c*one))*(diag((y-(one`*y/n)*one)*(y-(one`*y/n)*one)`)*one)`*c*one/((y-(one`*y/n)*one)`*(y-(one`*y/n)*one)) - ((n-1)/n)*mc;
print gro segr z pr;

vars={y};
create step3 var vars;
append;

quit;


/*
data step3(replace=yes); set step3; set step1(keep=sid); run;
proc spatialreg data=step3 Wmat=swm corrb; model y = /type=sMA; spatialid sid; run;
proc spatialreg data=step3 Wmat=swm corrb; model y = /type=sAr; spatialid sid; run;


proc spatialreg data=step1 Wmat=swm outest=parms(replace=yes) corrb; *bounds _lambda<0; model y = /type=SARMA; *restrict _rho=0.592885; spatialid sid; run;

*/
