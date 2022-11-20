*******************************
TRG 880 EXAM u18119543
Non-Linear Regression and bootstrap;

quit;
proc import datafile="C:/Users/Matthew/Documents/1 UNI STUFF/TRG880/TRG_EXAM ACTUAL/QUESTION 1/q1.csv"
dbms=csv out=nlin_dat replace;
run;

proc print data=nlin_dat (obs=10);
run;

proc sgplot data=nlin_dat;
title "Initial Graphical Representation of Time Series Data";
	scatter x=t y=pt/ filledoutlinedmarkers 
   markerfillattrs=(color=blue) 
   markeroutlineattrs=(color=blue thickness=5)
   markerattrs=(symbol=circlefilled size=5);
run;
******USING NLIN TO CHECK;
proc nlin data=nlin_dat  ; 
	parms K=15 P0=10 r=10 ;
	model pt = (K*P0)/(P0+(K-P0)*exp(-r*t));
	run ;

****NLIN;
	quit;
proc iml;

/* _FUN and _DER are text strings that define model and deriv */
/* _parm contains parm names                                  */
/* _beta contains initial values for parameters               */
/* _k is the number of parameters                             */

*THIS IS WHERE WE DEFINE OUR FUNCTION;
start nlinit;
   _dep = "Pt";                      /* dependent variable */
   _fun = "(K*P0)/(P0+(K-P0)*exp(-r*t))"; /* nonlinear regression model */
   _der  = {"(P0**2 * (exp(r*t)-1)*exp(r*t))/(K+P0*exp(r*t)-P0)**2",
			"(K**2 * exp(r*t))/((exp(r*t)-1)*P0+K)**2",
			"-((K*P0)*(P0-K)*t*exp(r*t))/(P0*exp(r*t)-P0+K)**2"};
   _parm = {"K", "P0", "r"};               /* names of parameters */
   _beta = {45, 15, 2};          *INITIAL PARAMETER VALUES. NOTE! IF THESE ARE NOT 'CLOSE ENOUGH' THEN THE PROCESS WILL FAIL;
   _k= nrow(_parm);                   /* number of parameters */
finish nlinit;


start nlgen;
   call change(_fun, '*', '#', 0);  /* substitute '#' for '*' */
   call change(_der, '*', '#', 0);

   /* Write the NLFIT module at run time */
   call queue('start nlfit;');
   do i=1 to _k;
      call queue(_parm[i], "=_beta[", char(i,2), "];");
   end;
   call queue("_y = ", _dep, ";",
              "_p = ", _fun, ";",
              "_r = _y - _p;",
              "_sse = ssq(_r);",
              "finish;" );

   /* Write the NLDERIV function at run time */
   call queue('start nlderiv; free _NULL_; _x = ');
   do i=1 to _k;
      call queue("(", _der[i], ")||");
   end;
   call queue("_NULL_; finish;");

   call queue("resume;");    /* Pause to compile the functions */
   pause *;
finish nlgen;

/*NRA HAPPENS HERE*/
start nlest;
   run nlfit;                /* f, r, and sse for initial beta */

   /* Gauss-Newton iterations to estimate parameters */
   do _iter=1 to 30 until(_eps<1e-8);
      run nlderiv;             /* subroutine for derivatives */
      _lastsse = _sse;
      _xpxi=sweep(_x`*_x);
      _delta = _xpxi*_x`*_r;            /* correction vector */
      _old = _beta;              /* save previous parameters */
      _beta = _beta + _delta;        /* apply the correction */
      run nlfit;                         /* compute residual */
      _eps = abs((_lastsse-_sse)) / (_sse+1e-6);
      /* Hartley subiterations */
      do _subit=1 to 10 while(_sse>_lastsse);
         _delta = _delta/2;   /* halve the correction vector */
         _beta = _old+_delta; /* apply the halved correction */
         run nlfit;                        /* find sse et al */
      end;
      /* if no improvement after 10 halvings, exit iter loop */
      if _subit>10 then _eps=0;
   end;

   /* display table of results  */
   if _iter < 30 then do; /* convergence */
      _dfe = nrow(_y) - _k;
      _mse = _sse/_dfe;
      _std = sqrt(vecdiag(_xpxi)#_mse);
      _t = _beta/_std;
      _prob = 1 - cdf("F", _t#_t, 1, _dfe);
      print _beta[label="Estimate"] _std[label="Std Error"]
            _t[label="t Ratio"]     _prob[format=pvalue6.];
      print _iter[label="Iterations"] _lastsse[label="SSE"];
   end;
   else print "Convergence failed";
finish nlest;


use nlin_dat;
read all var "pt" into Pt;
read all var "t"  into t;
read all var "year" into year;



run nlinit; /* define strings that define the regression model */
run nlgen;  /* write modules that evaluate the model           */
run nlest;  /* compute param estimates, std errs, and p-values */


**********************COMPUTING ESTIMATES;
Parm_beta=_beta;
print Parm_beta;
K=Parm_beta[1];
P0=Parm_beta[2];
r=Parm_beta[3];

Pt_hat=(K*P0)/(P0+(K-P0)*exp(-r*t));
print year t Pt Pt_hat;
cnn={"year" "t" "Pt" "Pt_hat"};
nlin_est=year||t||Pt||Pt_hat;

start corr_yx(y,x);
	p=ncol(x);
	joint=y||x;
	corr_check=corr(joint);
	correl=corr_check[1,2:(p+1)];
	return (correl);
finish;
cor_p=corr_yx(Pt,Pt_hat);
R2=cor_p**2;
print cor_p R2;

create nlin_out from nlin_est[c=cnn];
append from nlin_est;
close nlin_out;
quit;


*****BOOTSTRAP;
proc iml;

use nlin_dat;
read all var "pt" into Pt;
read all var "t"  into t;
read all var "year" into year;


/* _FUN and _DER are text strings that define model and deriv */
/* _parm contains parm names                                  */
/* _beta contains initial values for parameters               */
/* _k is the number of parameters                             */

*THIS IS WHERE WE DEFINE OUR FUNCTION;
start nlinit;
   _dep = "Pt";                      /* dependent variable */
   _fun = "(K*P0)/(P0+(K-P0)*exp(-r*t))"; /* nonlinear regression model */
   _der  = {"(P0**2 * (exp(r*t)-1)*exp(r*t))/(K+P0*exp(r*t)-P0)**2",
			"(K**2 * exp(r*t))/((exp(r*t)-1)*P0+K)**2",
			"-((K*P0)*(P0-K)*t*exp(r*t))/(P0*exp(r*t)-P0+K)**2"};
   _parm = {"K", "P0", "r"};               /* names of parameters */
   _beta = {90, 20, 0.1};          *INITIAL PARAMETER VALUES. NOTE! IF THESE ARE NOT 'CLOSE ENOUGH' THEN THE PROCESS WILL FAIL;
   _k= nrow(_parm);                   /* number of parameters */
finish nlinit;


start nlgen;
   call change(_fun, '*', '#', 0);  /* substitute '#' for '*' */
   call change(_der, '*', '#', 0);

   /* Write the NLFIT module at run time */
   call queue('start nlfit;');
   do i=1 to _k;
      call queue(_parm[i], "=_beta[", char(i,2), "];");
   end;
   call queue("_y = ", _dep, ";",
              "_p = ", _fun, ";",
              "_r = _y - _p;",
              "_sse = ssq(_r);",
              "finish;" );

   /* Write the NLDERIV function at run time */
   call queue('start nlderiv; free _NULL_; _x = ');
   do i=1 to _k;
      call queue("(", _der[i], ")||");
   end;
   call queue("_NULL_; finish;");

   call queue("resume;");    /* Pause to compile the functions */
   pause *;
finish nlgen;

/*NRA HAPPENS HERE*/
start nlest;
   run nlfit;                /* f, r, and sse for initial beta */

   /* Gauss-Newton iterations to estimate parameters */
   do _iter=1 to 30 until(_eps<1e-8);
      run nlderiv;             /* subroutine for derivatives */
      _lastsse = _sse;
      _xpxi=sweep(_x`*_x);
      _delta = _xpxi*_x`*_r;            /* correction vector */
      _old = _beta;              /* save previous parameters */
      _beta = _beta + _delta;        /* apply the correction */
      run nlfit;                         /* compute residual */
      _eps = abs((_lastsse-_sse)) / (_sse+1e-6);
      /* Hartley subiterations */
      do _subit=1 to 10 while(_sse>_lastsse);
         _delta = _delta/2;   /* halve the correction vector */
         _beta = _old+_delta; /* apply the halved correction */
         run nlfit;                        /* find sse et al */
      end;
      /* if no improvement after 10 halvings, exit iter loop */
      if _subit>10 then _eps=0;
   end;

   /* display table of results  */
   if _iter < 30 then do; /* convergence */
      _dfe = nrow(_y) - _k;
      _mse = _sse/_dfe;
      _std = sqrt(vecdiag(_xpxi)#_mse);
      _t = _beta/_std;
      _prob = 1 - cdf("F", _t#_t, 1, _dfe);
/*      print _beta[label="Estimate"] _std[label="Std Error"]*/
/*            _t[label="t Ratio"]     _prob[format=pvalue6.];*/
/*      print _iter[label="Iterations"] _lastsse[label="SSE"];*/
   end;
   else print "Convergence failed";
finish nlest;

start corr_yx(y,x);
	p=ncol(x);
	joint=y||x;
	corr_check=corr(joint);
	correl=corr_check[1,2:(p+1)];
	return (correl);
finish;



pt_org=PT;
t_org=t;
n=nrow(t); ss_n=n;
ss_m=J(ss_n,1,0);
do i=1 to ss_n;
	ss_m[i,]=i;
end;

/*print ss_m;*/
sam_size=round(0.8*n);
print sam_size;

NSAM=500; *Number of Bootstraps;
cor_res=J(NSAM,1,.);
R2_res=J(NSAM,1,.);
Rat_res=J(NSAM,1,.);
do jj=1 to NSAM;

	PT=pt_org;
	t=t_org;
	ss=sample((ss_m),sam_size, "NoRep");
	ss=ss`;
	call sort(ss);
	t_sam=t[ss,];
	PT_sam=pt[ss,];

	t=t_sam;
	pt=pt_sam;

	run nlinit; 
	run nlgen;  
	run nlest;  

	Parm_beta=_beta;
	K=Parm_beta[1];
	P0=Parm_beta[2];
	r=Parm_beta[3];

	Pt_hat=(K*P0)/(P0+(K-P0)*exp(-r*t));
	cor_p=corr_yx(Pt,Pt_hat);
	R2=cor_p**2;
	Rat=K/P0;

	cor_res[jj]=cor_p;
	R2_res[jj]=R2;
	rat_res[jj]=Rat;

end;

*cor_res CONFIDENCE;
call sort(cor_res);
	a=1-0.10;
	p_vec={0.05,0.95};
    call qntl(CI_COR,cor_res,p_vec);
	print "Confidence Interval:" a ;
	print p_vec CI_COR;

*R2_res CONFIDENCE;
call sort(R2_res);
	a=1-0.10;
	p_vec={0.05,0.95};
    call qntl(CI_R2,R2_res,p_vec);
	print "Confidence Interval:" a ;
	print p_vec CI_R2;

*RATIO CONFIDENCE;
call sort(rat_res);
	a=1-0.10;
	p_vec={0.05,0.95};
    call qntl(CI_ratio,rat_res,p_vec);
	print "Confidence Interval:" a ;
	print p_vec CI_ratio;

quit;

proc sgplot data=nlin_out;
title "Graphical Representation of Initial and Estimated Population Growth";
	scatter x=t y=pt/ filledoutlinedmarkers 
   markerfillattrs=(color=blue) 
   markeroutlineattrs=(color=blue thickness=5)
   markerattrs=(symbol=circlefilled size=5);
   	scatter x=t y=Pt_hat/ filledoutlinedmarkers 
   markerfillattrs=(color=red) 
   markeroutlineattrs=(color=red thickness=5)
   markerattrs=(symbol=circlefilled size=5);
run;



*POINTWISE PREDICTION;
proc iml;
use nlin_dat;
read all var "pt" into Pt;
read all var "t"  into t;
read all var "year" into year;
pt_org=pt;
pt_max=max(pt_org);
bnd=J(10,1,.);
do jj=1 to 10;
	bnd[jj]=pt_max-(10*ranuni(5));
end;

pt_morg=pt_org//bnd;


start nlinit;
   _dep = "Pt";                      /* dependent variable */
   _fun = "(K*P0)/(P0+(K-P0)*exp(-r*t))"; /* nonlinear regression model */
   _der  = {"(P0**2 * (exp(r*t)-1)*exp(r*t))/(K+P0*exp(r*t)-P0)**2",
			"(K**2 * exp(r*t))/((exp(r*t)-1)*P0+K)**2",
			"-((K*P0)*(P0-K)*t*exp(r*t))/(P0*exp(r*t)-P0+K)**2"};
   _parm = {"K", "P0", "r"};               /* names of parameters */
   _beta = {90, 20, 0.1};          *INITIAL PARAMETER VALUES. NOTE! IF THESE ARE NOT 'CLOSE ENOUGH' THEN THE PROCESS WILL FAIL;
   _k= nrow(_parm);                   /* number of parameters */
finish nlinit;


start nlgen;
   call change(_fun, '*', '#', 0);  /* substitute '#' for '*' */
   call change(_der, '*', '#', 0);

   call queue('start nlfit;');
   do i=1 to _k;
      call queue(_parm[i], "=_beta[", char(i,2), "];");
   end;
   call queue("_y = ", _dep, ";",
              "_p = ", _fun, ";",
              "_r = _y - _p;",
              "_sse = ssq(_r);",
              "finish;" );

   call queue('start nlderiv; free _NULL_; _x = ');
   do i=1 to _k;
      call queue("(", _der[i], ")||");
   end;
   call queue("_NULL_; finish;");

   call queue("resume;");    /* Pause to compile the functions */
   pause *;
finish nlgen;

/*NRA HAPPENS HERE*/
start nlest;
   run nlfit;                

   /* Gauss-Newton iterations to estimate parameters */
   do _iter=1 to 30 until(_eps<1e-8);
      run nlderiv;             
      _lastsse = _sse;
      _xpxi=sweep(_x`*_x);
      _delta = _xpxi*_x`*_r;            
      _old = _beta;             
      _beta = _beta + _delta;       
      run nlfit;                       
      _eps = abs((_lastsse-_sse)) / (_sse+1e-6);
      do _subit=1 to 10 while(_sse>_lastsse);
         _delta = _delta/2;   
         _beta = _old+_delta; 
         run nlfit;                        
      end;
      if _subit>10 then _eps=0;
   end;


   if _iter < 30 then do; 
      _dfe = nrow(_y) - _k;
      _mse = _sse/_dfe;
	  _cvar = _xpxi;
      _std = sqrt(vecdiag(_xpxi)#_mse);
      _t = _beta/_std;
      _prob = 1 - cdf("F", _t#_t, 1, _dfe);
   end;
   else print "Convergence failed";
finish nlest;


run nlinit; /* define strings that define the regression model */
run nlgen;  /* write modules that evaluate the model           */
run nlest;  /* compute param estimates, std errs, and p-values */


**********************COMPUTING ESTIMATES;
Parm_beta=_beta;
Parm_std=_std;
Parm_cov=_cvar;
print Parm_cov; *VARIANCE//COVARIANCE MATRIX;
print Parm_beta Parm_std;
K=Parm_beta[1];
P0=Parm_beta[2];
r=Parm_beta[3];
t_max=max(t);
y_max=max(year);
y_min=min(year);
t_newl=do(0,(t_max+10),1); t_new=t_newl`;
y_newl=do(y_min,(y_max+10),1); y_new=y_newl`;

print t_new y_new;


parm_var=sum(vecdiag(parm_cov));
parm_std=sqrt(parm_var);
print parm_var parm_std;
Pt_hat=(K*P0)/(P0+(K-P0)*exp(-r*t_new));
Poi_var=J(nrow(Pt_hat),1,.);
do ii=1 to nrow(Pt_hat);
	t_int=t_new[ii];
	poi_var[ii]=t_int*parm_var*t_int;
end;
poi_std=sqrt(poi_var);
print Poi_var poi_std;
Pt_below=Pt_hat-(1.96)#poi_std;
Pt_above=Pt_hat+(1.96)#poi_std;
print Pt_below Pt_hat Pt_above ;
pwise_est=Pt_below||Pt_hat||Pt_above||t_new||y_new||pt_morg;
cnn={"Lower_CI" "Pt_hat" "Upper_CI" "t" "year" "Original_Data"};
create pwise_out from pwise_est[c=cnn];
append from pwise_est;
close pwise_out;

quit;

proc sgplot data=pwise_out;
title "Graphical Representation of point-wise estimates";
	series x=year y=Lower_CI/ lineattrs=(color=red thickness=2);
	series x=year y=Pt_hat/ lineattrs=(color=black thickness=4);
	series x=year y=Upper_CI/ lineattrs=(color=red thickness=2);
	scatter x=year y=Original_Data/ filledoutlinedmarkers 
   markerfillattrs=(color=blue) 
   markeroutlineattrs=(color=blue thickness=4)
   markerattrs=(symbol=circlefilled size=4);
run;
