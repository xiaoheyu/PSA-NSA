filename indata 'H:\SA\data.csv';
filename sparseij 'H:\SA\Rook.csv';
OPTIONS LINESIZE=72;
TITLE 'PSA-NSA mixtures: SAR-SMA';

data step1; infile indata dlm=','; input EV1 PopDen2 tr_PopDen2 tr_Respira tr_Medin new_id; 
if new_id='.' then delete;
y=tr_Medin;
sid=new_id;
run;
proc standard data=step1 out=step1(replace=yes) mean=0; var y; run;
proc means noprint; var y; output out=outn n=n; run;
data outn(replace=yes); set outn(keep=n); run;
proc univariate data=step1 normal; var y; run;

data swm; infile sparseij dlm=','; input NEW_ID NID WEIGHT;
i=new_id; j=nid;
if i='.' then delete;
drop NEW_ID NID WEIGHT;
run;
PROC IML;
use swm; read all into ij;
use outn; read all into n;
print n;
c=shape(0,n,n);
indx =n*(ij[,1]-1) + ij[,2];
c[indx]=1;
*print c;

IDEN=I(n);
ONE=J(n,1);
sumc=one`*c*one;
print n sumc;

mcm=shape(0,n,n);
mco=shape(0,n,1);
mcm=(iden-one*one`/n)*c*(iden-one*one`/n);
call eigen(evals,evecs,mcm);
mco=n*evals/(one`*c*one);
mcmax=mco[1];
mcmin=mco[n];
print mcmin mcmax;
id=shape(0,n,1);
do i=1 to n;
id[i]=i;
if mco[i]>0 then mco[i]=mco[i]/mcmax;
if mco[i]<0 then mco[i]=mco[i]/abs(mcmin);
end;
print id mco;

create step2 from evecs;
append from evecs;
quit;

data step1(replace=yes); set step1; set step2; run;
proc glmselect maxmacro=1 outdesign=outmesf seed=1234567;
    model y=col1-col99 col212-col402/selection=stepwise(select=sl sle=0.1 sls=0.1001 choose=cv) cvMethod=random(13) stats=press;
	 output out=tempesf p=yhat resid=yr;
run;
proc reg; model y=yhat; run;
proc univariate normal; var yr; run;
proc reg; model y=COL2 COL3 COL4 COL5 COL6 COL7 COL8 COL9 COL10 COL12 COL13 COL14 COL15 COL16 COL18 COL19 COL20 COL25 COL26 COL29 COL32 
COL39 COL42 COL45 COL46 COL47 COL49 COL50 COL53 COL54 COL56 COL61 COL64 COL65 COL67 COL72 COL73 COL77 COL79 COL81 COL98; run;
proc reg; model y=COL213 COL219 COL223 COL227 COL232 COL236 COL238 COL248 COL249 COL250 COL252 COL254 COL260 COL271 COL275 COL294 COL306 
COL309 COL311 COL324 COL332 COL335 COL337 COL338 COL341 COL342 COL343 COL351 COL359 COL361 COL380 COL383 COL389; run;

data x; set outmesf; drop y; run;
proc standard data=tempesf out=tempesf mean=0; var yhat; run;
data esf; set tempesf(keep=yhat); run;
data esfr; set tempesf(keep=yr); run;
PROC IML;
use swm; read all into ij;
use outn; read all into n;
use esf; read all into esf;
use esfr; read all into yr;
use x; read all into x;
print n;
c=shape(0,n,n);
indx =n*(ij[,1]-1) + ij[,2];
c[indx]=1;
ONE=J(n,1);

y=yr;
mc=n/(one`*c*one)*(y-(one`*y/n)*one)`*c*(y-(one`*y/n)*one)/((y-(one`*y/n)*one)`*(y-(one`*y/n)*one));
gr=((n-1)/(2*one`*c*one))*(one`*diag(c*( (y*one`-one*y`)#(y*one`-one*y`))) *one)/((y-(one`*y/n)*one)`*(y-(one`*y/n)*one));
print mc gr;
p=ncol(x)-1;
emc=-(n/(one`*c*one))*trace((inv(x`*x))*x`*c*x)/(n-p-1);
smc=sqrt(2/(one`*c*one));
zmc=(mc-emc)/smc;
p=1-cdf('normal',zmc);
print emc smc zmc p;

y=esf;
mc=n/(one`*c*one)*(y-(one`*y/n)*one)`*c*(y-(one`*y/n)*one)/((y-(one`*y/n)*one)`*(y-(one`*y/n)*one));
gr=((n-1)/(2*one`*c*one))*(one`*diag(c*( (y*one`-one*y`)#(y*one`-one*y`))) *one)/((y-(one`*y/n)*one)`*(y-(one`*y/n)*one));
print mc gr;

quit;

*/
