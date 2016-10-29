$Title BESTA Model _main structure

$STitle Example model definitions
$eolcom //
*option profile =3;
Option Limrow=0;
option Limcol=0;

$include 1.multi_Set.gms
$include 3.muti_Data.gms
$include 2.multi_Assumption.gms //this contains scalar maxdis (maximum distance) and picked candidate
$include 4.multi_Parameter.gms



************ 7. OBJECTIVE FUNCTION AND CONSTRAINTS ************

Variables
z(k) Total cost and GHG emissions to support a  ethanol plant
A(i,crop,b) Acres of switchgrass produced annually
XC(i,crop,b) Tons of Switchgrass produced annually
AH(i,m,crop,b) Acres of land harvested through Nov to Feb
XH(i,m,crop,b) Tons of switchgrass harvested annually (Nov ~ Feb)
NXS(i,m,crop,b,t) Tons of switchgrass newly stored in each month (Nov ~ Feb)
XS(i,m,crop,b,t) Tons of switchgrass cumulatively stored in each month (Nov ~ Oct)
XTN(i,m,crop,b) Tons of newly harvested switchgrass that is transported in month "m"
XTO(i,m,crop,b,t) Tons of stored switchgrass that is transported in month "m"
Q(m) Monthly output of ethanol
LOA  Total acreage harvested
LOXH Total tonnage harvested
Numb(b,m) Number of balers used in harvest
Numl(b,m) Number of loaders used in harvest
Numw(b,m) Number of mowers used in harvest
Numrk(b,m) Number of rakes used in harvest
Numt(b,m) Number of tractors used in harvest
truckload(i,m) truckload used for delievery
;

Positive Variables A, AH, XC, XH, NXS,XS, XTN, XTO, Q
;
Binary variable NumbBio number of biorefineries built;


Equations
OBJECTIVE1  Suppy chain cost
OBJECTIVE2  Greenhouse gas emission
Acreage(i,crop) Acreage converted from cropland "crop" and available in Hexagon "i"
Cultivate(i,crop,b) Tonnage of switchgrass available in Hexagon "i"
HarLoss(i,crop,b) No tonnage loss occurs during harvest activity (production tonnage = harvest tonnage)
AAH(i,m,crop,b) Acreage harvested*Yield = Tonnage harvested
Acreage3 Total acerage of switchgrass harvested
Tonnage1 Total tonnage of switchgrass harvested
Tonnage2(m) The ratio of monthly harvested tonnage to the total harvest tonnage = the ratio of available monthly working hours of harvest machinery over the total working hours of each harvest machinery
Fieldayb(b,m) Machine time of baler is under the available field hours in each month
Fieldayl(b,m) Machine time of loader is under the available field hours in each month
Fieldayw(b,m) Machine time of mower is under the available field hours in each month
Fieldayk(b,m) Machine time of raker is under the available field hours in each month
Fieldayt(b,m) Machine time of tractor is under the available field hours in each month
HarTran1(i,m,crop,b) Newly harvested switchgrass should be no less than newly harvested switchgrass that is transported
IHST1(i,m,crop,b) Tonnage of newly stored switchgrass equals the total tonnage harvested minus tonnage sent to the biorefinery directly
IHST2(i,m,crop,b) Tonnage harvested is greater than the tonnage put into storage in month "m"
IHST3(i,m,crop,b,t) Tonnage of newly stored switchgrass equals to the cumulative stored tonnage in November
HST1(i,m,crop,b,t) Tonnage of cumulative storage in month "November"
HST2(i,m,crop,b,t) Tonnage of cumulative storage in month "m" of harvest season
HST3(i,m,crop,b,t) Tonnage of cumulative storage in month "m" of off-harvest season
StoTran1(i,m,crop,b,t) Tonnage of storage in month "m" is greater than the tonnage sent from the storage in month "m"
Conversion(m) Quantity of ethanol produced
Demand(m) Quantity of ethanol produced is no less than the demand of ethanol in each month
delievery(i,m) define truckload used in delievery
;


OBJECTIVE1.. z('cost') =e= sum((i,crop,b),BEPswi(i,crop,b)*XC(i,crop,b))
             + sum((i,m,crop,b,t),Gamma(i,b,t)* NXS(i,m,crop,b,t))
             + sum((i,b),Theta1(i,b)*(sum((m,crop),XTN(i,m,crop,b))+ sum((m,crop,t),XTO(i,m,crop,b,t)))/(1-DMLtrans))
;

OBJECTIVE2.. z('CO2emission') =e= sum(crop,sum((i,m,b),AH(i,m,crop,b))*(CO2(crop)+N2O(crop)))
             + sum((i,m,crop,b),XH(i,m,crop,b))*StorE+ sum(b,sum((i,m,crop),AH(i,m,crop,b))*(HarE(b)+EME(b)))
             + sum((i,m,crop,b),AH(i,m,crop,b)*(FertiE+HerbiE+SeedE))
             + sum((i,m),truckload(i,m)*TransE2(i,m))
;

*** Production
Acreage(i,crop)..                 sum(b,A(i,crop,b)) =l= avaAcrecrops(i,crop) ;
Cultivate(i,crop,b)..             XC(i,crop,b) =l= A(i,crop,b)*YieldS(i,b);


*** Harvest
HarLoss(i,crop,b)..               XC(i,crop,b) =e= sum(m,XH(i,m,crop,b))/(1-DML_har) ;
AAH(i,m,crop,b)..                 XH(i,m,crop,b) =e= AH(i,m,crop,b)*YieldS(i,b);
Acreage3..                        sum((i,m,crop,b),AH(i,m,crop,b)) =e= LOA ;
Tonnage1..                        sum((i,m,crop,b),XH(i,m,crop,b)) =e= LOXH ;
Tonnage2(m)$m2(m)..               sum((i,crop,b),XH(i,m,crop,b)) =e= CapUnit/Lambda*rateava(m) ;
Fieldayb(b,m)..                   sum((i,crop),MTb(i,b)*AH(i,m,crop,b)) =l= Numb(b,m)*avaday(m) ;
Fieldayl(b,m)..                   sum((i,crop),MTl(i,b)*AH(i,m,crop,b)) =l= Numl(b,m)*avaday(m) ;
Fieldayw(b,m)..                   sum((i,crop),MTw(i,b)*AH(i,m,crop,b)) =l= Numw(b,m)*avaday(m) ;
Fieldayk(b,m)..                   sum((i,crop),MTk(i,b)*AH(i,m,crop,b)) =l= Numrk(b,m)*avaday(m) ;
Fieldayt(b,m)..                   sum((i,crop),MTt(i,b)*AH(i,m,crop,b)) =l= Numt(b,m)*avaday(m) ;


*** Transportation and Storage during Harvest Season
HarTran1 (i,m,crop,b)..           XH(i,m,crop,b) =g= XTN(i,m,crop,b) ;
IHST1(i,m,crop,b)$m4(m)..         sum((t),NXS(i,m,crop,b,t)) =e= XH(i,m,crop,b)- XTN(i,m,crop,b)/(1-DMLtrans) ;
IHST2(i,m,crop,b)..               XH(i,m,crop,b) =g= sum(t,NXS(i,m,crop,b,t)) ;
IHST3(i,m,crop,b,t)$m1(m)..       NXS(i,m,crop,b,t) =e= XS(i,m,crop,b,t) ;
HST1(i,m,crop,b,t)$m1(m)..        XS(i,m+1,crop,b,t)=e= (1-slm(b,t,m))*XS(i,m,crop,b,t)- XTO(i,m+1,crop,b,t)/(1-DMLtrans)
                                                         + NXS(i,m+1,crop,b,t) ;
HST2(i,m,crop,b,t)$m2(m)..        XS(i,m+1,crop,b,t) =e= (1-slm(b,t,m))*XS(i,m,crop,b,t)- XTO(i,m+1,crop,b,t)/(1-DMLtrans)
                                                         + NXS(i,m+1,crop,b,t) ;

*** Transportation and Storage during Off-Harvest Season
HST3(i,m,crop,b,t)$m3(m)..        XS(i,m+1,crop,b,t) =e= (1-slm(b,t,m))*XS(i,m,crop,b,t)- XTO(i,m+1,crop,b,t)/(1-DMLtrans) ;
StoTran1(i,m,crop,b,t)..          XS(i,m,crop,b,t) =g= XTO(i,m+1,crop,b,t)/(1-DMLtrans) ;

*** Biorefinery Demand
Conversion(m)..                   Q(m) =e= Lambda*(sum((i,crop,b),XTN(i,m,crop,b))+ sum((i,crop,b,t),XTO(i,m,crop,b,t))) ;
Demand(m)..                       Q(m) =e= Dd*Cap(m) ;

delievery(i,m)..                  truckload(i,m) =e=  sum(b,(sum((crop,t),XTO(i,m,crop,b,t))+ sum((crop),XTN(i,m,crop,b)))/loadwt2(i,m,b));


option      MIP=CPLEX ;
Option      Limrow=0;
option      limcol=0;
Option IterLim = 1000000;
option reslim = 1000000;

Model example / all /;

*---------------------------------------------------------------------------
$STitle eps-constraint method

Set k1(k) the first element of k, km1(k) all but the first elements of k;
k1(k)$(ord(k)=1) = yes; km1(k)=yes; km1(k1) = no;
Set kk(k)     active objective function in constraint allobj

Parameter
   rhs(k)     right hand side of the constrained obj functions in eps-constraint
   maxobj(k)  maximum value from the payoff table
   minobj(k)  minimum value from the payoff table

scalar
iter   total number of iterations
infeas total number of infeasibilities
elapsed_time elapsed time for payoff and e-sonstraint
start start time
finish finish time

Variables
   a_objval   auxiliary variable for the objective function
   obj        auxiliary variable during the construction of the payoff table
Positive Variables
   sl(k)      slack or surplus variables for the eps-constraints
Equations
   con_obj(k) constrained objective functions
   augm_obj   augmented objective function to avoid weakly efficient solutions
   allobj     all the objective functions in one expression;

con_obj(km1)..   z(km1) - dir(km1)*sl(km1) =e= rhs(km1);

* We optimize the first objective function and put the others as constraints
* the second term is for avoiding weakly efficient points
augm_obj..
  sum(k1,dir(k1)*z(k1))+1e-3*sum(km1,sl(km1)/(maxobj(km1)-minobj(km1))) =e= a_objval;

allobj..  sum(kk, dir(kk)*z(kk)) =e= obj;

Model mod_payoff    / example, allobj / ;
Model mod_epsmethod / example, con_obj, augm_obj / ;

Parameters
  payoff(k,k)  payoff tables entries
  jpayoff(j,k,k)
  payoffindi(k,k) includes the indirect GHG emissions
  jpayoffindi(j,k,k) includes the indirect GHG emissions
  hayratio(k)
  jhayratio(j,k);

Alias(k,kp);

set dynamic(i,j);

loop(j,
Theta1(i,b)= Theta(i,j,b) ;
dynamic(i,j)=(RegionDistance(i,j)>maxdis);
transE2(i,m)=TransEE(i,j,m)/1000;

XTN.fx(i,m,crop,b1)=0;   //b1=0 means square bale; b2=0 means round bale
XTO.fx(i,m,crop,b1,t)=0;
XTO.fx(i,m,crop,b,t2)=0; //t1=0 means untarp bale; t2=0 means tarp bale
XH.fx(i,m3,crop,b)=0;    //no harvest during off-harvest season
XTO.fx(i,m1,crop,b,t)=0; //no inventory left. stored switchgrass deliveried in the first month is 0

XH.fx(i,m,crop,b)$dynamic(i,j)=0;
AH.fx(i,m,crop,b)$dynamic(i,j)=0;
XTN.fx(i,m,crop,b)$dynamic(i,j)=0;
XTO.fx(i,m,crop,b,t)$dynamic(i,j)=0;
NXS.fx(i,m,crop,b,t)$dynamic(i,j)=0;
XS.fx(i,m,crop,b,t)$dynamic(i,j)=0;
mod_payoff.holdfixed=1;

  loop(kp,
    kk(kp)=yes;
    solve mod_payoff using mip maximizing obj;
    payoff(kp,k) = z.l(k);
    payoffindi(kp,'cost')=z.l('cost');
    payoffindi(kp,'co2emission')=z.l('co2emission')
                               +sum(b,sum(m,Numb.L(b,m))*BalerE(b))+sum((b,m),Numl.L(b,m))*LoadE
                               +sum((b,m),Numw.L(b,m))*MowerE + sum((b,m),Numrk.L(b,m))*RakeE
                               +sum((b,m),Numt.L(b,m))*TracE ;
    hayratio(kp)= sum((i,m,b),AH.l(i,m,'hay',b))/sum((i,m,crop,b),AH.l(i,m,crop,b));
    kk(k++1) = kk(k);
    );
  jpayoff(j,kp,k)= payoff(kp,k);
  jpayoffindi(j,kp,k)= payoffindi(kp,k);
  jhayratio(j,kp)=hayratio(kp)
  if (mod_payoff.modelstat<>1 and mod_payoff.modelstat<>8,
abort 'no optimal solution for mod_payoff');

option clear=XH;
option clear=AH;
option clear=XTN;
option clear=XTO;
option clear=NXS;
option clear=XS;  // release the constrant of variables for the next loop
);

Set g grid points /g0*g%gridpoints%/
    g1(g) grid points subset /g0,g%gridpoints%/
    g2(g) grid points subset2
    grid(k,g) grid
Parameter
    gridrhs(k,g) rhs of eps-constraint at grid point
    jgridrhs(j,k,g) j dimension added
    maxg(k) maximum point in grid for objective
    firstOffMax, lastZero some counters
    numg(g) ordinal value of g starting with 0
    kpayoff(k,k)
    landuseratio(g)
    jlanduseratio(j,g) the percentage of hayland in total acres of land used;
numg(g) = ord(g)-1;
grid(km1,g) = yes;
g2(g)=yes;
g2(g1)=no;


parameter solutions(g,k), solutions2(j,g,k), AHL(i,m,crop,b),consl(g,k),jconsl(j,g,k);
loop(j,
 kpayoff(kp,k)=jpayoff(j,kp,k);
 minobj(k)=smin(kp,kpayoff(kp,k));
 maxobj(k)=smax(kp,kpayoff(kp,k));
 maxg(km1) = smax(grid(km1,g), numg(g));
 jgridrhs(j,km1,g)$(dir(km1)=-1) = (maxobj(km1) - numg(g)/maxg(km1)*(maxobj(km1)- minobj(km1)));
 jgridrhs(j,km1,g)$(dir(km1)=1) = (minobj(km1) + numg(g)/maxg(km1)*(maxobj(km1)- minobj(km1)));
 jgridrhs(j,km1,g)$(numg(g)=maxg(km1) and dir(km1)=-1) = minobj(km1)+0.5;
 jgridrhs(j,km1,g)$(numg(g)=maxg(km1) and dir(km1)=1) = maxobj(km1)-0.5;

);


loop(j,
 gridrhs(km1,g)=jgridrhs(j,km1,g);
 Theta1(i,b)= Theta(i,j,b) ;
 transE2(i,m)=TransEE(i,j,m)/1000;
 dynamic(i,j)=(RegionDistance(i,j)>maxdis);

XTN.fx(i,m,crop,b1)=0;
XTO.fx(i,m,crop,b1,t)=0;
XTO.fx(i,m,crop,b,t2)=0;
XH.fx(i,m3,crop,b)=0;
XTO.fx(i,m1,crop,b,t)=0;

XH.fx(i,m,crop,b)$dynamic(i,j)=0;
AH.fx(i,m,crop,b)$dynamic(i,j)=0;
XTN.fx(i,m,crop,b)$dynamic(i,j)=0;
XTO.fx(i,m,crop,b,t)$dynamic(i,j)=0;
NXS.fx(i,m,crop,b,t)$dynamic(i,j)=0;
XS.fx(i,m,crop,b,t)$dynamic(i,j)=0;
mod_epsmethod.holdfixed=1;

    loop(g2,
    rhs(km1)=gridrhs(km1,g2);
    solve mod_epsmethod maximizing a_objval using mip;
    consl(g2,k)=z.l(k);
    solutions(g2,'cost')=z.l('cost');
    solutions(g2,'CO2emission')=z.l('CO2emission')
                               +sum(b,sum(m,Numb.L(b,m))*BalerE(b))+sum((b,m),Numl.L(b,m))*LoadE
                               +sum((b,m),Numw.L(b,m))*MowerE + sum((b,m),Numrk.L(b,m))*RakeE
                               +sum((b,m),Numt.L(b,m))*TracE;
    landuseratio(g2)=sum((i,m,b),AH.l(i,m,'hay',b))/sum((i,m,crop,b),AH.l(i,m,crop,b));
    );
    solutions2(j,g2,k)=solutions(g2,k);
    jconsl(j,g2,k)=consl(g2,k);
    jlanduseratio(j,g2)=landuseratio(g2);
    AHL(i,m,crop,b)=AH.l(i,m,crop,b);
option clear=XH;
option clear=AH;
option clear=XTN;
option clear=XTO;
option clear=NXS;
option clear=XS;
);

solutions2(j,'g0',k)=jpayoffindi(j,'cost',k);
solutions2(j,'g%gridpoints%',k)=jpayoffindi(j,'co2emission',k);
jlanduseratio(j,'g0')=jhayratio(j,'cost');
jlanduseratio(j,'g%gridpoints%')=jhayratio(j,'co2emission');

set upper /1*10000/;
parameter solutions3(upper,k);

loop ((j,g),
solutions3(upper,k)$(ord(upper)eq ord(j)*10+ord(g))=solutions2(j,g,k));

set jg(upper);
jg(upper)$(solutions3(upper,'co2emission')>0)=yes;


alias(jg,jgd);
parameter count(upper,upper),count2(upper);

loop(jg,
 loop(jgd,
   count(jg,jgd)$(solutions3(jg,'cost')-solutions3(jgd,'cost')>=0
   and solutions3(jg,'CO2emission')-solutions3(jgd,'CO2emission')>=0
   and(solutions3(jg,'cost')-solutions3(jgd,'cost')+solutions3(jg,'CO2emission')-solutions3(jgd,'CO2emission'))>0)=1;
 );
);

count2(jg)=sum(jgd,count(jg,jgd));

set jg2(upper);
jg2(upper)$(count2(upper)>=1)=yes;

set pareto(upper);
pareto(jg)=yes; pareto(jg2)=no;

parameter optimalgrids(upper,k);
optimalgrids(pareto,k)=solutions3(pareto,k);

execute_unload 'AllTN6_12' payoff, jpayoff,solutions, solutions2, solutions3,count, count2,optimalgrids, jconsl;
