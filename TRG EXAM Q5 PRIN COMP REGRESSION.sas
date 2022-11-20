********************************
TRG880 PRIN COMP REGRESSION;
*******************************;

proc import datafile="C:/Users/Matthew/Documents/1 UNI STUFF/TRG880/TRG_EXAM ACTUAL/QUESTION 5/q5.csv"
dbms=csv out=pcr_init replace;
run;

proc print data=pcr_init (obs=10);
run;

proc reg data=pcr_init;
model y = x1 x2 x3 x4 x5 x6 x7 x8/ vif tol collin;
title 'OLS Regression for VIF';
run;

proc sgplot data=pcr_init;
title "Scatter Plot of Y and X1";
	scatter y=y x=x1;
run;

proc sgplot data=pcr_init;
title "Scatter Plot of Y and X2";
	scatter y=y x=X2;
run;

proc sgplot data=pcr_init;
title "Scatter Plot of Y and X3";
	scatter y=y x=X3;
run;
proc sgplot data=pcr_init;
title "Scatter Plot of Y and X4";
	scatter y=y x=X4;
run;
proc sgplot data=pcr_init;
title "Scatter Plot of Y and X5";
	scatter y=y x=X5;
run;
proc sgplot data=pcr_init;
title "Scatter Plot of Y and X6";
	scatter y=y x=X6;
run;
proc sgplot data=pcr_init;
title "Scatter Plot of Y and X7";
	scatter y=y x=X7;
run;
proc sgplot data=pcr_init;
title "Scatter Plot of Y and X8";
	scatter y=y x=X8;
run;

proc iml;
*****INVESTIGATING MC;
use pcr_init;
read all into yx;

n=nrow(yx);
p=ncol(x_orig);
y_orig=yx[,ncol(yx)];
x_orig=yx[,1:(ncol(yx)-1)];

x=x_orig;
y=y_orig;
print (y_orig[1:5,]) (x_orig[1:5,]);




start LR(bh, r2, vif, se, t, pval, x , y);
    n=nrow(x);
	p=ncol(x);
    xm=j(n,1,1)||x;
	bh=inv(xm`*xm)*xm`*y;
	yh=xm*bh;
	sse=(y-yh)`*(y-yh);
	sigma2=(1/(n-p-1))*sse;
	se_beta=inv(xm`*xm)*sigma2;
	se=sqrt(vecdiag(se_beta));
	t=bh/se;
	df=j(p,1,1);
	pval=2*(1-probt(abs(t), n-p));
	tss=sum((y-mean(y))##2);
	r2=1-(sse/tss);
	vif=j(p+1, 1, .);
	vif[1]=0;
	do jj=1 to p;
	   id=remove(1:p, jj);
	   yi=x[,jj];
	   xi=j(n,1,1)||x[,id];
	   bhi=inv(xi`*xi)*xi`*yi;
	   yhi=xi*bhi;
	   ssei=(yi-yhi)`*(yi-yhi);
	   tssi=sum((yi-mean(yi))##2);
	   r2i=1-(ssei/tssi);
	   store=jj+1;
	   vif[store]=1/(1-r2i);
	end;
finish;

*Simple Linear Regression;
x=x_orig;
corx=corr(x);
n=nrow(x);
p=ncol(x);
y=y_orig;
nmI="Intercept";
nm1="X1":"X8";
nmq=nmI||nm1;
run LR(bhat, r2, vif, se, t ,pval, x , y);
print bhat[r=nmq c="BHAT"] se t pval vif r2;

S=cov(x_orig);
R=cov2corr(S);
print S[L="Covariance" r=nm1 c=nm1];
print R[L="Correlation" r=nm1 c=nm1];
quit;


******BEGINING PCR REGRESSION;

proc iml;
use pcr_init;
read all into yx;
n=nrow(yx);
p=ncol(x_orig);
y_orig=yx[,ncol(yx)];
x_orig=yx[,1:(ncol(yx)-1)];


xc=(x_orig-mean(x_orig))/(std(x_orig));
yc=(y_orig-mean(y_orig))/(std(y_orig));
S=cov(xc);
R=cov2corr(S);
y=y_orig;


print "Covariance:" S;
print "Correlation:" R;
call eigen (D,V,R);
v_nam={"x1" "x2" "x3" "x4" "x5" "x6" "x7"};
print D[r=v_nam] V[r=v_nam c=("Prin1":"Prin7")];
nfac=ncol(xc)-0;
print "Number of Prin Comps:" nfac;
eig_v=V[,1:(ncol(xc)-0)];
print eig_v;
Z=xc*eig_v;
print "Prin Comps:" (z[1:5,]);
corr_z=corr(z);
print corr_z;
bh_z=inv(z`*z)*z`*yc;


bh_x=eig_v*inv(Z`*Z)*Z`*yc;
beta0=mean(y_orig)-mean(x_orig)*beta;
bh=beta0//bh_x;
nmx="Bx1":"Bx8";
nmz="Bz1":"Bz8";
print bh_z[r=nmz] bh_x[r=nmx];
quit;
*CHECKING RESULTS;
proc standard data=pcr_init out=pcr_init_std m=0 s=1;
var x1-x5 y;
run ;
***PCA REGRESSION;
title "PRINCIPAL COMPONENT REGRESSION";
proc reg data=pcr_init_std plots=none outest=PE; 
   model y = x1 x2 x3 x4 x5 x6 x7 x8  / PCOmit=2;       
quit;
 
proc print data=PE(where=(_Type_="IPC")) noobs;
   var Intercept--x8;
run;




******STEPWISE PCR;
*ADD ONE AT A TIME;

proc iml;
use pcr_init;
read all into yx;
n=nrow(yx);
p=ncol(x_orig);
y_orig=yx[,ncol(yx)];
x_orig=yx[,1:(ncol(yx)-1)];


xc=(x_orig-mean(x_orig))/(std(x_orig));
yc=(y_orig-mean(y_orig))/(std(y_orig));
S=cov(xc);
R=cov2corr(S);
y=y_orig;


print "Covariance:" S;
print "Correlation:" R;
call eigen (D,V,R);
v_nam={"x1" "x2" "x3" "x4" "x5" "x6" "x7" "x8"};
print D[r=v_nam] V[r=v_nam c=("Prin1":"Prin8")];
*THIS IS THE FULL EIGENVECTOR;

*THIS FUNCTION WILL CALCULATE THE PROPORTION OF VARIANCE IN THE PREDICTORS EXPLAINED BY THE PRINCOMPS ;
start pred_var(part_ev,eigvals);
eigsum=eigvals[+,];
en=nrow(part_ev);
eig_prop={};
do i=1 to en;
	e_cons=part_ev[i,];
	int=e_cons/eigsum;
	eig_prop=eig_prop//int;
end;
	cum_ep=cusum(eig_prop);
	return(cum_ep);
finish;
eig_prop=pred_var(D[1:3],D);
print eig_prop;

betaz_res=J(ncol(V),nrow(V),.);
betax_res=J(ncol(V),nrow(V),.);
predvar_res=J(nrow(V),1,.);
respvar_res=J(nrow(V),1,.);

stp={};
do jj=1 to ncol(V);
	stp=stp||J(1,1,jj);
	nfac=jj;
	print "Number of Prin Comps:" nfac;
	eig_v=V[,1:nfac];
	val_d=D[1:nfac,];
	predvar=pred_var(val_d,D);
	predvar_res[jj,]=predvar[jj,];
	print eig_v;
	Z=xc*eig_v;
	print "Prin Comps:" (z[1:5,]);
	corr_z=corr(z);
	print corr_z;
	bh_z=inv(z`*z)*z`*yc;

	yh_z=z*bh_z;
	sse=(yc-yh_z)`*(yc-yh_z);
	tss=sum((yc-mean(yc))##2);
	r2=1-(sse/tss);


	respvar_res[jj,]=r2;

	bh_x=eig_v*inv(Z`*Z)*Z`*yc;
	beta0=mean(y_orig)-mean(x_orig)*bh_x;
	bh=beta0//bh_x;
	nmx="Bx1":"Bx8";
	nmz="Bz1":"Bz8";
	nmr="Step1":"Step8";
	print bh_z[r=nmz] bh_x[r=nmx];
	dim=nrow(bh_z);
	betaz_res[1:dim,jj]=bh_z;
	betax_res[,jj]=bh_x;
end;

print predvar_res[r=nmr] respvar_res;
var_res1=stp`||predvar_res||respvar_res;
cnn={"Step" "Predictor_Var" "Response_Var"};

create var_res1dat from var_res1[c=cnn];
append from var_res1;
close var_res1dat;

print betaz_res[c=nmr r=nmz];
print betax_res[c=nmr r=nmx];
qn="Step"//nmz`//nmx`;
mix_parm=stp//betaz_res//betax_res;
mix_parm=mix_parm`;
print mix_parm[c=qn];

create mix_dat from mix_parm[c=qn];
append from mix_parm;
close mix_dat;


quit;

proc print data=mix_dat;
run;

proc sgplot data=var_res1dat;
title "Proportion of Variance Explain";
	series y=Predictor_Var x=step/ lineattrs=(thickness=2);
	series y=Response_Var x=step/ lineattrs=(thickness=2);
run;

proc sgplot data=mix_dat;
title "Tracking Changes in Parameter Values";
	series y=bz1 x=step/ lineattrs=(thickness=2);
	series y=bz2 x=step/ lineattrs=(thickness=2);
	series y=bz3 x=step/ lineattrs=(thickness=2);
	series y=bz4 x=step/ lineattrs=(thickness=2);
	series y=bz5 x=step/ lineattrs=(thickness=2);
	series y=bz6 x=step/ lineattrs=(thickness=2);
	series y=bz7 x=step/ lineattrs=(thickness=2);
	series y=bz8 x=step/ lineattrs=(thickness=2);
	series y=bx1 x=step/ lineattrs=(pattern=dash thickness=4);
	series y=bx2 x=step/ lineattrs=(pattern=dash thickness=4);
	series y=bx3 x=step/ lineattrs=(pattern=dash thickness=4);
	series y=bx4 x=step/ lineattrs=(pattern=dash thickness=4);
	series y=bx5 x=step/ lineattrs=(pattern=dash thickness=4);
	series y=bx6 x=step/ lineattrs=(pattern=dash thickness=4);
	series y=bx7 x=step/ lineattrs=(pattern=dash thickness=4);
	series y=bx8 x=step/ lineattrs=(pattern=dash thickness=4);
run;



****FEATURE SELECTION;
proc iml;

use pcr_init;
read all into yx;

n=nrow(yx);
p=ncol(x_orig);
y_orig=yx[,ncol(yx)];
x_orig=yx[,1:(ncol(yx)-1)];

x=x_orig;
y=y_orig;




start LR(bh, r2, vif, se, t, pval, x , y);
    n=nrow(x);
	p=ncol(x);
    xm=j(n,1,1)||x;
	bh=inv(xm`*xm)*xm`*y;
	yh=xm*bh;
	sse=(y-yh)`*(y-yh);
	sigma2=(1/(n-p-1))*sse;
	se_beta=inv(xm`*xm)*sigma2;
	se=sqrt(vecdiag(se_beta));
	t=bh/se;
	df=j(p,1,1);
	pval=2*(1-probt(abs(t), n-p));
	tss=sum((y-mean(y))##2);
	r2=1-(sse/tss);
	vif=j(p+1, 1, .);
	vif[1]=0;
	do jj=1 to p;
	   id=remove(1:p, jj);
	   yi=x[,jj];
	   xi=j(n,1,1)||x[,id];
	   bhi=inv(xi`*xi)*xi`*yi;
	   yhi=xi*bhi;
	   ssei=(yi-yhi)`*(yi-yhi);
	   tssi=sum((yi-mean(yi))##2);
	   r2i=1-(ssei/tssi);
	   store=jj+1;
	   vif[store]=1/(1-r2i);
	end;
finish;

*Simple Linear Regression;
x=x_orig;
corx=corr(x);
n=nrow(x);
p=ncol(x);
y=y_orig;
nmI="Intercept";
nm1="X1":"X8";
nmq=nmI||nm1;


parms_res=J(ncol(x)+1,ncol(x)-1,.);
vif_res=J(ncol(x)+1,ncol(x)-1,.);
rsq_res=J(ncol(x),1,.);
do ii=2 to ncol(x);
	x_int=x_orig[,1:ii];
	x=x_int;
	run LR(bhat, r2, vif, se, t ,pval, x , y);
	dim=nrow(bhat);
	nmr="Step2":"Step8";
	nmb="B0":"B8";
	parms_res[1:dim,(ii-1)]=bhat;
	vif_res[1:dim,(ii-1)]=vif;
	rsq_res[ii]=r2;
end;

print parms_res[c=nmr r=nmb];
print vif_res[c=nmr r=nmb];
nmrr="Step1"||nmr;
cn={"RSQ"};
print vif_res[c=nmr r=nmb] rsq_res[r=nmrr c=cn];

quit;

******************************************
POST FEATURE SELECTION;
*****************************************;


proc iml;
use pcr_init;
read all into yx;
n=nrow(yx);
p=ncol(x_orig);
y_orig=yx[,ncol(yx)];
x_orig=yx[,1:2];


xc=(x_orig-mean(x_orig))/(std(x_orig));
yc=(y_orig-mean(y_orig))/(std(y_orig));
S=cov(xc);
R=cov2corr(S);
y=y_orig;


print "Covariance:" S;
print "Correlation:" R;
call eigen (D,V,R);
v_nam={"x1" "x2" "x3" "x4" "x5" "x6" "x7" "x8"};
print D[r=v_nam] V[r=v_nam c=("Prin1":"Prin8")];
*THIS IS THE FULL EIGENVECTOR;

*THIS FUNCTION WILL CALCULATE THE PROPORTION OF VARIANCE IN THE PREDICTORS EXPLAINED BY THE PRINCOMPS ;
start pred_var(part_ev,eigvals);
eigsum=eigvals[+,];
en=nrow(part_ev);
eig_prop={};
do i=1 to en;
	e_cons=part_ev[i,];
	int=e_cons/eigsum;
	eig_prop=eig_prop//int;
end;
	cum_ep=cusum(eig_prop);
	return(cum_ep);
finish;
eig_prop=pred_var(D[1:3],D);
print eig_prop;

betaz_res=J(ncol(V),nrow(V),.);
betax_res=J(ncol(V),nrow(V),.);
predvar_res=J(nrow(V),1,.);
respvar_res=J(nrow(V),1,.);

stp={};
do jj=1 to ncol(V);
	stp=stp||J(1,1,jj);
	nfac=jj;
	print "Number of Prin Comps:" nfac;
	eig_v=V[,1:nfac];
	val_d=D[1:nfac,];
	predvar=pred_var(val_d,D);
	predvar_res[jj,]=predvar[jj,];
	print eig_v;
	Z=xc*eig_v;
	print "Prin Comps:" (z[1:5,]);
	corr_z=corr(z);
	print corr_z;
	bh_z=inv(z`*z)*z`*yc;

	yh_z=z*bh_z;
	sse=(yc-yh_z)`*(yc-yh_z);
	tss=sum((yc-mean(yc))##2);
	r2=1-(sse/tss);


	respvar_res[jj,]=r2;

	bh_x=eig_v*inv(Z`*Z)*Z`*yc;
	beta0=mean(y_orig)-mean(x_orig)*bh_x;
	bh=beta0//bh_x;
	nmx="Bx1":"Bx2";
	nmz="Bz1":"Bz2";
	nmr="Step1":"Step8";
	print bh_z[r=nmz] bh_x[r=nmx];
	dim=nrow(bh_z);
	betaz_res[1:dim,jj]=bh_z;
	betax_res[,jj]=bh_x;
end;

print predvar_res[r=nmr] respvar_res;
var_res1=stp`||predvar_res||respvar_res;
cnn={"Step" "Predictor_Var" "Response_Var"};

create var_res1dat from var_res1[c=cnn];
append from var_res1;
close var_res1dat;

print betaz_res[c=nmr r=nmz];
print betax_res[c=nmr r=nmx];
qn="Step"//nmz`//nmx`;
mix_parm=stp//betaz_res//betax_res;
mix_parm=mix_parm`;
print mix_parm[c=qn];

create mix_dat from mix_parm[c=qn];
append from mix_parm;
close mix_dat;


quit;

proc print data=mix_dat;
run;

proc sgplot data=var_res1dat;
title "Proportion of Variance Explain";
	series y=Predictor_Var x=step/ lineattrs=(thickness=2);
	series y=Response_Var x=step/ lineattrs=(thickness=2);
run;
