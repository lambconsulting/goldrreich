

%let main =U:\Consulting\KEL\Fox\Goldreich;
%let indir =&main.\input;
%let outdir =main\output;
libname gold "&indir." ;
libname out "&outdir.";
/*filename mylst "&outdir.\my_listing.lst";*/
/*filename mylog "&outdir.\my_listing.log";*/

/* Original individual isolate data prior to 20200526 phone conversation and drop cefazolin to avoid numeric character conflict */
/*data gold; set gold.gold; run;*/
data gold_indiv; set gold.gold (drop=cefazolin); run;

/* Bring in total isolate data per 20200526 phone conversation */
data total; set gold.total; run;

/* Stack individual bacteria and total isolate */
data gold; set gold_indiv total; run;





/*cefazolin*/
%let cont_vars = bacitracin	chloramphenicol ciprofloxacin erythromycin
gentamicin marbofloxacin neomycin ofloxacin oxytetracycline
polymyxin_b	tobramycin;
%let cont_vars_all = &cont_vars. cefazolin;
%let bac = staph strep pseudo past;


%macro antibiods;
%let ds_vars = &cont_vars. ;
%let ds_count=%sysfunc(countw(&ds_vars.));
%do i = 1 %to &ds_count.; * &ds_count.;
data xxx (drop=%scan(&ds_vars,&i.)); length antibio_nm $50.;
set gold (keep=yr id bacteria %scan(&ds_vars,&i.)) ;
antibio_nm = "%scan(&ds_vars,&i.)";
antibio = %scan(&ds_vars,&i.);
run;
data antibio; set antibio xxx; run;
%end;
%mend;
data antibio; set _null_; run;
%antibiods;

ods pdf file="C:\Junk\panel_survs.pdf";
proc sgpanel data= antibio;
by id;
panelby id / /*layout=lattice novarname */ uniscale= row ;
vline yr / response=antibio group=antibio_nm;
/*rowaxis fitpolicy = splitalways;*/
/*colaxis fitpolicy =staggerthin;*/
/*axis interval=even;*/

run;

ods pdf close;

proc freq data=gold;
tables yr bacteria;
run;

/* Regression */

%macro mixed(iv,bacteria,cvars);
/*%let ds_vars = strep staph past pseudo ;*/
%let contvars = &cvars.;
/*%let ds_count=%sysfunc(countw(&ds_vars.));*/
%let count=%sysfunc(countw(&contvars.));
/*%do j = 1 %to &ds_count.; * &ds_count.;*/
data ds; set gold;
where bacteria ="&bacteria."; run;
%do i = 1 %to &count.; /* &count. */
/*ods trace on;*/
proc mixed data=ds;
model %scan(&contvars,&i.) = &iv. / solution;
ods output solutionf=sol tests3=t3;
run;
data sol; length ds dv iv $50.; set sol; where effect not in ("Intercept");
ds = "&bacteria."; dv = "%scan(&contvars,&i.)"; iv = "&iv."; run;
/*ods trace off;*/
data reg (keep=ds dv iv Estimate StdErr tvalue fvalue ProbF); merge sol t3; run;
data reg2 (rename=(tvalue = reg_t_test fvalue=type3_f_test)); set reg; run;
data reg_all; set reg_all reg2; run;
data _null_;
   set reg;
   call symput('eqn',"Parameter estimate="||estimate);
run;
%put &eqn.;
proc sgplot data=ds;
ODS GRAPHICS / RESET IMAGENAME = "%scan(&contvars,&i.)_&bacteria." IMAGEFMT =PNG ;
ODS LISTING IMAGE_DPI= 300 GPATH = 'C:\Junk' ; 
reg y =  %scan(&contvars,&i.) x = &iv. / clm;
   inset "&eqn." "DV = %scan(&contvars,&i.)" "Bacteria = &bacteria.";
run;
	%end;
proc datasets nolist nodetails; delete ds sol reg reg2 t3; run;

%mend;
data reg_all; set _null_; run;



%mixed(yr,staph,&cont_vars);
%mixed(yr,strep,&cont_vars);
%mixed(yr,pseudo,&cont_vars);
%mixed(yr,past,&cont_vars);
%mixed(yr,total,&cont_vars_all);




%macro penorm (ds,iv,bacteria);

/*%let ds_vars = strep staph past pseudo ;*/
%let contvars = &cont_vars_all.;
/*%let ds_count=%sysfunc(countw(&ds_vars.));*/
%let count=%sysfunc(countw(&contvars.));
/*%do j = 1 %to &ds_count.; * &ds_count.;*/
data ds; set &ds.;
where bacteria ="&bacteria."; run;

%do i = 1 %to &count.;
/*ods trace on;*/
proc reg data=ds;
  model %scan(&contvars,&i.) = &iv. ;
  output out=%scan(&contvars,&i.)res rstudent=%scan(&contvars,&i.)r  ;
  ods output ParameterEstimates=%scan(&contvars,&i.)pe FitStatistics = fs;
run;
quit;
%if %sysfunc(exist(%scan(&contvars,&i.)pe)) %then %do;

data r_val (keep=r2 cvalue2); set fs (keep=label2 cvalue2);
where label2 = "R-Square";
drop label2;
r2 = input(cvalue2,10.4);
run;


data pe0;
length var bacteria $50;
set %scan(&contvars,&i.)pe;
var = "%scan(&contvars,&i.)";
bacteria ="&bacteria.";
where variable not in ("Intercept");
run;

data pe_r; merge pe0 r_val; run;


data pe;
set pe pe_r;
run;


* Examine residuals for normality ;
proc univariate data=%scan(&contvars,&i.)res plots plotsize=30 normal;
  var %scan(&contvars,&i.)r;
  ods output TestsForNormality=%scan(&contvars,&i.)norm;
run;
data %scan(&contvars,&i.)norm (drop=varname);
length var bacteria $50;
set %scan(&contvars,&i.)norm;
var = "%scan(&contvars,&i.)";
bacteria ="&bacteria.";
run;
data norm;
set norm %scan(&contvars,&i.)norm;
run;

%end;
	%end;
proc datasets nolist nodetails; delete ds ; run;

%mend;
data norm; set _null_; data pe; set _null_; data pe_norm; set _null_; data npar; set _null_; 
data pe_r; set _null_; run;


%penorm(gold,yr,staph);
%penorm(gold,yr,strep);
%penorm(gold,yr,pseudo);
%penorm(gold,yr,past);
%penorm(gold,yr,total);