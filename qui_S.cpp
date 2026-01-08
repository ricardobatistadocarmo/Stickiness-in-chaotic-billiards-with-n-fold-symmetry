#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXRAND (4294967296)  
#define FNORM (2.3283064365e-10)  
#define RANDOM ((ira[ip++] = ira[ip1++] + ira[ip2++]) ^ ira[ip3++])  
#define FRANDOM (FNORM * RANDOM)  
#define pm1 ((FRANDOM > 0.5) ? 1 : -1)  
#define pm0 ((FRANDOM > 0.5) ? 1 : 0) 
 
unsigned myrand, ira[256];  
unsigned char ip, ip1, ip2, ip3;  
 
unsigned rand4init(void);  
void Init_Random(void);  
 
unsigned rand4init (void) {

  unsigned long long y;

  y = myrand * 16807LL;
  myrand = (y & 0x7fffffff) + (y >> 31);

  if (myrand & 0x80000000) myrand = (myrand & 0x7fffffff) + 1;
  
  return myrand;

  }  
 
void Init_Random(void) {

  int i;
     
  ip = 128;
  ip1 = ip - 24;
  ip2 = ip - 55;
  ip3 = ip - 61;

  for (i=ip3; i<ip; i++) ira[i] = rand4init();

  }

double R,r,PI,sqrt3;
int i,j,n,sectors,lad, pointsR,pontosline,pointsr,pointstotal,ladaux,rest,divi,sectorcol,sectoraux,cont,ok,sectorprev,ladmax,k;
double aloop[10000], bloop[10000],a,b,xv[100000],yv[100000],xc[100000],yc[100000],gam1,gam2,alf,t4,t3,t2,t1,dx;
double  dtR,dtr,xrot,yrot;
double dify,difx;
double vx, vy, ytraj,xtraj, atraj,btraj, xtest, ytest, xcollision, ycollision;
double x1,x2,xfinal,yfinal,ar[100000],br[100000],xsol[5],ysol[5],ysolr[5],dist[5],dmin;
double y11, y22,deriy,angy, vxfinal,vyfinal,xmin,xmax,ymin,ymax,sum,per,pertotal,vtan,angmax;
double ncvel, varvel, dvel,ncper,varper,dper,medrel,fillingbox,cont_med,xi,yi,distreg,sumreg,disttotal,sumtotal, sumchaos,medchaos,medreg,qsi[100],sumqsi,dregfinal[100],dchaoticfinal[100],sumregfinal,sumchaosfinal, rho, mediafinalqsi,sumdeviation,deviationqsi,deviation,mediafinalreg,mediafinalcao,deviationreg,deviationcao;
int box[1001][1001],ifinal,jfinal,ICs,IC,IC_cont,timetraj,conttrajfinal,contfinal,contstick,count, contreg, contchaos ;

double *xdisc, *ydisc, *coefang, *coelinear, *ang, *sum_mr, *xvet, *yvet ;
int *icont, *jcont, *vetsector;

int pborder;

int main(){

  pborder = 50000001;
  
  xdisc = (double *) malloc(pborder * sizeof(double));
  ydisc = (double *) malloc(pborder * sizeof(double));
  coefang = (double *) malloc(pborder * sizeof(double));
  coelinear = (double *) malloc(pborder * sizeof(double));
  ang = (double *) malloc(pborder * sizeof(double));
  sum_mr = (double *) malloc(pborder * sizeof(double));
  icont = (int *) malloc(pborder * sizeof(int));
  jcont = (int *) malloc(pborder * sizeof(int));
  vetsector = (int *) malloc(pborder * sizeof(int));
  xvet = (double *) malloc(pborder * sizeof(double));
  yvet = (double *) malloc(pborder * sizeof(double));

  myrand=1314;  //Seed for random number generator. int number.
  Init_Random();

  timetraj = 20000000; //Maximum number of collisions
  ICs = 5; //Number of initial conditions

  PI=3.1415926535897932384626433832795028841971693993751058209749445923;
  sqrt3=sqrt(3);	

  R = 0.5;  //Radius of one of the circunference that make up the boundary
  r = 0.25; // Radius of the second circunference
  n = 2; // Billiard symmetry

  sectors = 3*n;  //lad from 1 to 3*n

  xc[1] = 1-R;  //before xc5
  yc[1] = 0;  //after yc5

  xc[sectors] = 1 - r;  //before xc3. Sectors is the last lad in counter-clockwise
  yc[sectors] = 0;  //antes yc3

  xc[3] = xc[sectors]*cos(2*PI/n) - yc[sectors]*sin(2*PI/n); 
  yc[3] = xc[sectors]*sin(2*PI/n) + yc[sectors]*cos(2*PI/n); 

  printf("%d  \n",sectors);	

  FILE* output;
  output = fopen("c2_stick_20_1000x1000_20M_.dat","w");	

  FILE* output3;
  output3 = fopen("c2_rho_part5.dat","w");

  FILE *input1;
  input1 = fopen("a.dat","r"); //File containing the values of the angular coefficients (slopes) of the lines that connect the circunferences

  FILE *input2;
  input2 = fopen("b.dat","r"); //File containing the values of the linear coefficients (intercepts) of the lines that connect the circunferences

// The choice of values for a and b will depend on the chosen symmetry n

  j=1;
  while (j<1001){
  
    fscanf(input1,"%lf",&aloop[j]);
    j++;

    }

  i=1;
  while (i<1001){
  
    fscanf(input2,"%lf",&bloop[i]);
    i++;

    }
	
  a = aloop[n];		
  b = bloop[n];
	
  //Three basal points
  gam1 = ( (3/(r*r))  + (1/(r*r))   )*a - ( -(sqrt3/(r*r))  + (sqrt3/(r*r))    );
  gam2 = ( (3/(r*r))  + (1/(r*r))   )*(b-yc[3]) + ( -(sqrt3/(r*r))  + (sqrt3/(r*r)))*xc[3];
  alf = 4/(r*r*r*r);

  xv[1] = 1;
  yv[1] = 0
  
  xv[2] = -(a*b - xc[1]*pow(R/R,2) )/( a*a + pow(R/R,2)   );
  yv[2] =  sqrt( R*R - pow((xv[2]-xc[1])*(R/R),2)  )  ;

  xv[3] =  - (gam1*gam2 - 4*xc[3]*alf)/(gam1*gam1 + 4*alf);
  yv[3] =    (( - ((sqrt3*(xv[3]-xc[3]))/(r*r))  + (sqrt3/(r*r))*(xv[3]-xc[3]) + 2*sqrt( (3/(r*r)) + (1/(r*r)) - (4/(r*r*r*r))*pow(xv[3]-xc[3],2)    ))/( (3/(r*r)) + (1/(r*r))  ))   + yc[3]  ;

   //Points between sectors' border
  i = 4;
  while(i<(sectors+1)){
  
    xv[i] = xv[i-3]*cos(2*PI/n) - yv[i-3]*sin(2*PI/n); 
    yv[i] = xv[i-3]*sin(2*PI/n) + yv[i-3]*cos(2*PI/n);
    
    i++;
    
    }
	
  xv[sectors+1] = xv[1];
  yv[sectors+ 1] = yv[1];	
	
  //Discretization of sectros' border 1, 2, and 3
  t4 = fabs(acos((xv[4]-xc[3])/r));

  t3 = fabs(acos((xv[3]-xc[3])/r));

  t2 = fabs(acos((xv[2]-xc[1])/R));

  t1 =0;

  printf(" tempos %.12f %.12f %.12f %.12f \n",t4, t3, t2,t1);

  pointsR = 400 ; //100 for almost symmetries
  dtR = (t1-t2)/pointsR;

  i=1;
  xdisc[i] = xv[1];
  ydisc[i] = yv[1];
  
  while(i<(pointsR+1)){
  
    xdisc[i+1] = xc[1]  + R*cos(t1 - dtR*i);
    ydisc[i+1] = sqrt( R*R - (xdisc[i+1] - xc[1] )*(xdisc[i+1] - xc[1] )   ) +yc[1] ;
    
    i++;
    
    }	

  //Discretization of sector 2
  pontosline = 400;
  dx = (xv[2] - xv[3])/pontosline;

  i =  1;
  xdisc[pointsR + i] = xv[2];
  ydisc[pointsR + i] = yv[2];

  while(i<(pontosline+1)){
  
    xdisc[i+pointsR+1] = xv[2] - dx*i;
    ydisc[i+pointsR+1] = a*xdisc[i+pointsR+1] + b;
    
    i++;
    	
    }

  //Discretization of sector 3
  pointsr = 400 ;
  dtr = (t3-t4)/pointsr;

  i=1;
  xdisc[i+ pointsR+ pontosline] = xv[3];
  ydisc[i+ pointsR+ pontosline] = yv[3];

  while(i<(pointsr+1)){
  
    xdisc[i+pointsR+pontosline+1] = xc[3]  + r*cos(t3 - dtr*i);
    ydisc[i+pointsR+pontosline+1] = sqrt( r*r - (xdisc[i+pointsR+pontosline+1] - xc[3] )*(xdisc[i+pointsR+pontosline+1] - xc[3] )   ) +yc[3] ;
  
    i++;
  
    }		

  pointstotal = pointsR+pontosline+pointsr;	
	
  //Rotation
  j = 1;
  while(j<n){
  
    i=1;
    while(i<(pointsR + pointsr + pontosline +1) ){
    
      xrot = xdisc[i]*cos(2*PI*j/n) - ydisc[i]*sin(2*PI*j/n);
      yrot = xdisc[i]*sin(2*PI*j/n) + ydisc[i]*cos(2*PI*j/n);
      
      xdisc[(pointsR + pointsr + pontosline)*j + i] = xrot;
      ydisc[(pointsR + pointsr + pontosline)*j + i] = yrot;
      
      i++;
      
      }
      
    j++;
	
    }	

  xdisc[n*pointstotal +1] = xdisc[1];
  ydisc[n*pointstotal +1] = ydisc[1];

  i = 1;	
  while(i<(  n*pointstotal + 1 )  ){
  
    coefang[i] = 0;
    coelinear[i] = 0;
    ang[i] = 0;
    
    i++;
    
    }
    
  i = 1;
  angmax=0;
  while(i<(  n*pointstotal + 1 )  ){

    dify = (ydisc[i+1] - ydisc[i]);
    difx = (xdisc[i+1] - xdisc[i] );
    coefang[i] = dify/difx;
       
    coelinear[i] = ydisc[i] - coefang[i]*xdisc[i];
    ang[i] = atan(coefang[i]); 

    if((ang[i]>angmax) && (i<(pointstotal*n))){
  
      angmax = ang[i];
      ladmax = i;
	
      }
    
    if(i==(n*pointstotal)){
  
      coefang[i] = (ydisc[1] - ydisc[i])/(xdisc[1] - xdisc[i] );
      coelinear[i] = ydisc[i] - coefang[i]*xdisc[i];
      ang[i] = atan(coefang[i]);
    
      }
    
    i++;
    
    }		

  //Radius 'R' circunferences centers
  i = 1;
  while(i<n){
  
    xc[1+3*i] = xc[3*i-2]*cos(2*PI/n) - yc[3*i-2]*sin(2*PI/n); 
    yc[1+3*i] = xc[3*i-2]*sin(2*PI/n) + yc[3*i-2]*cos(2*PI/n);
    
    i++;
    
    }	
	
	
  //Radius 'r' circunferences centers
  i = 1;
  while(i<(n-1)){
  
    xc[3+3*i] = xc[3*i]*cos(2*PI/n) - yc[3*i]*sin(2*PI/n); 
    yc[3+3*i] = xc[3*i]*sin(2*PI/n) + yc[3*i]*cos(2*PI/n);
    
    i++;
    
    }		
	
  i = 0;
  while(i<(timetraj+1)){
  
    sum_mr[i] = 0;
    i++;
    
    }	
	
  IC = 1;	

  sumqsi = 0;
  sumregfinal = 0;
  sumchaosfinal = 0;
	
  while(IC<(ICs+1)){

    sumreg = 0;
    sumtotal = 0;
    contreg=0;

    fillingbox = 0;

    i=1;
    while(i<1001){
    
      j=1;
      while(j<1001){
      
        box[i][j]=0;
        j++;
	
	}
	
      i++;
      
      }
  
    xi = FRANDOM*0.99998;
    yi = 0;
			
    vx = 2*FRANDOM -1;
    vy = sqrt(1-vx*vx);

    atraj =  (vy/vx);
    btraj = -atraj*xi;

    sectoraux = 0;
    cont = 1;

    while(cont<(timetraj+1)){

      ok = 0;

      i = 1;
      while(i<(n*pointstotal +1)){
      
        if(i==ladaux) i++;

	xtest = (btraj - coelinear[i] )/(coefang[i]- atraj);
	ytest = coefang[i]*xtest + coelinear[i];	

	if(xdisc[i+1]>xdisc[i]){
	
	  xmax = xdisc[i+1];
	  xmin = xdisc[i];
	
	  }

	if(xdisc[i+1]<xdisc[i]){
	
	  xmax = xdisc[i];
	  xmin = xdisc[i+1];

	  }

	if(ydisc[i+1]>ydisc[i]){
	
  	  ymax = ydisc[i+1];
	  ymin = ydisc[i];

	  }

	if(ydisc[i+1]<ydisc[i]){
	
	  ymax = ydisc[i];
	  ymin = ydisc[i+1];
	  
	  }

	if (( (xmin < xtest ) && (xmax > xtest )  )  && ( (ymin < ytest ) && (ymax > ytest )   ) ){
	
	  xcollision = xtest;
	  ycollision = ytest;
	  ladaux = i;
	  
	  rest = ladaux % 400; // The number 400 is the same one used in the discretization of each sector.
	  divi = ladaux / 400; // The number 400 is the same one used in the discretization of each sector.
	  
	  if(rest==0) sectoraux = divi;
	  if(rest!=0)msectoraux = divi + 1;
	  
	  j = 0;  //Interaction with radius 'R' circunferences sectors. Sectors 1, 4, 7, 10, ...
	  while(j<n){
	  
	    sectorcol = 1 + 3*j;
	    
	    if(sectorcol==sectoraux){
	      
	      if(sectorprev!=sectoraux){
	        
	        xsol[1] = (sqrt( -atraj*atraj*xc[sectorcol]*xc[sectorcol] + atraj*atraj*R*R - 2*atraj*btraj*xc[sectorcol] + 2*atraj*xc[sectorcol]*yc[sectorcol] - btraj*btraj + 2*btraj*yc[sectorcol] - yc[sectorcol]*yc[sectorcol] +R*R   )  -atraj*btraj  + atraj*yc[sectorcol]  + xc[sectorcol]  )/( atraj*atraj +1  );
	        xsol[2] = xsol[1];
	          
	          xsol[3] = (-sqrt( -atraj*atraj*xc[sectorcol]*xc[sectorcol] + atraj*atraj*R*R - 2*atraj*btraj*xc[sectorcol] + 2*atraj*xc[sectorcol]*yc[sectorcol] - btraj*btraj + 2*btraj*yc[sectorcol] - yc[sectorcol]*yc[sectorcol] +R*R   )  -atraj*btraj  + atraj*yc[sectorcol]  + xc[sectorcol]  )/( atraj*atraj +1  );
	        xsol[4] = xsol[3];
	          
	        ysol[1] = sqrt( R*R - (xsol[1] -xc[sectorcol] )*(xsol[1] -xc[sectorcol] )   )  + yc[sectorcol];
	          
	        ysol[2] = -sqrt( R*R - (xsol[2] -xc[sectorcol] )*(xsol[2] -xc[sectorcol] )   )  + yc[sectorcol];
	          
	        ysol[3] = sqrt( R*R - (xsol[3] -xc[sectorcol] )*(xsol[3] -xc[sectorcol] )   )  + yc[sectorcol];
	          
	        ysol[4] = -sqrt( R*R - (xsol[4] -xc[sectorcol] )*(xsol[4] -xc[sectorcol] )   )  + yc[sectorcol];
          
		dist[1] = sqrt(  (xsol[1]-xmin)*(xsol[1]-xmin) +  (ysol[1]-ymin)*(ysol[1]-ymin)   )	+  sqrt(  (xsol[1]-xmax)*(xsol[1]-xmax) +  (ysol[1]-ymax)*(ysol[1]-ymax)   );
		  
		dist[2] = sqrt(  (xsol[2]-xmin)*(xsol[2]-xmin) +  (ysol[2]-ymin)*(ysol[2]-ymin)   )	+  sqrt(  (xsol[2]-xmax)*(xsol[2]-xmax) +  (ysol[2]-ymax)*(ysol[2]-ymax)   );
		  
		dist[3] = sqrt(  (xsol[3]-xmin)*(xsol[3]-xmin) +  (ysol[3]-ymin)*(ysol[3]-ymin)   )	+  sqrt(  (xsol[3]-xmax)*(xsol[3]-xmax) +  (ysol[3]-ymax)*(ysol[3]-ymax)   );
		  
		dist[4] = sqrt(  (xsol[4]-xmin)*(xsol[4]-xmin) +  (ysol[4]-ymin)*(ysol[4]-ymin)   )	+  sqrt(  (xsol[4]-xmax)*(xsol[4]-xmax) +  (ysol[4]-ymax)*(ysol[4]-ymax)   );
		  
		ysolr[1] = atraj*xsol[1] + btraj;
		ysolr[2] = atraj*xsol[2] + btraj;
		ysolr[3] = atraj*xsol[3] + btraj;
		ysolr[4] = atraj*xsol[4] + btraj;
		  
		dmin = 9999;
		  
		if(fabs((ysolr[1]-ysol[1]))<0.00000001 ){
		  
		  if(dist[1]<dmin){
		    
		    dmin = dist[1];	   	  
		    xfinal = xsol[1];
		    yfinal = ysol[1];
		  
		    deriy = - (xfinal-xc[sectorcol])/( sqrt(R*R  - ( xc[sectorcol] - xfinal )*( xc[sectorcol] - xfinal )  )  ); 
	  	    angy = atan(deriy);
	  	      
	  	    }
	  	      
	  	  }
	  	    
	  	if(fabs((ysolr[2]-ysol[2]))<0.00000001 ){
	  	  
	  	  if(dist[2]<dmin){
	  	    
	  	    dmin = dist[2];
	  	    xfinal = xsol[2];
	  	    yfinal = ysol[2];
	  	    deriy =  (xfinal-xc[sectorcol])/( sqrt(R*R  - ( xc[sectorcol] - xfinal )*( xc[sectorcol] - xfinal )  )  );
	  	    angy = atan(deriy);
	  	      
	  	    }
	  	      
	  	  }
	  	    
	  	if(fabs((ysolr[3]-ysol[3]))<0.00000001 ){
	  	  
	  	  if(dist[3]<dmin){
	  	    
	  	    dmin = dist[3];
	  	    xfinal = xsol[3];
	  	    yfinal = ysol[3];
	  	    deriy = - (xfinal-xc[sectorcol])/( sqrt(R*R  - ( xc[sectorcol] - xfinal )*( xc[sectorcol] - xfinal )  )  );
	  	    angy = atan(deriy);
	  	      
	  	    }
	  	    
	  	  }
	  	  
	        if(fabs((ysolr[4]-ysol[4]))<0.00000001 ){
	  	
	  	  if(dist[4]<dmin){
	  	  
	  	    dmin = dist[4];
	  	    xfinal = xsol[4];
	  	    yfinal = ysol[4];
	  	    deriy =  (xfinal-xc[sectorcol])/( sqrt(R*R  - ( xc[sectorcol] - xfinal )*( xc[sectorcol] - xfinal )  )  );
	  	    angy = atan(deriy);
	  	    
	  	    }
	  	  
	  	  }
	  	
	        }
	      
	      if(sectorprev==sectoraux){
	    
	        xsol[1] = (sqrt( -atraj*atraj*xc[sectorcol]*xc[sectorcol] + atraj*atraj*R*R - 2*atraj*btraj*xc[sectorcol] + 2*atraj*xc[sectorcol]*yc[sectorcol] - btraj*btraj + 2*btraj*yc[sectorcol] - yc[sectorcol]*yc[sectorcol] +R*R   )  -atraj*btraj  + atraj*yc[sectorcol]  + xc[sectorcol]  )/( atraj*atraj +1  );
	        xsol[2] = xsol[1];
	      
	        xsol[3] = (-sqrt( -atraj*atraj*xc[sectorcol]*xc[sectorcol] + atraj*atraj*R*R - 2*atraj*btraj*xc[sectorcol] + 2*atraj*xc[sectorcol]*yc[sectorcol] - btraj*btraj + 2*btraj*yc[sectorcol] - yc[sectorcol]*yc[sectorcol] +R*R   )  -atraj*btraj  + atraj*yc[sectorcol]  + xc[sectorcol]  )/( atraj*atraj +1  );
	        xsol[4] = xsol[3];
	      
	        ysol[1] = sqrt( R*R - (xsol[1] -xc[sectorcol] )*(xsol[1] -xc[sectorcol] )   )  + yc[sectorcol];
	      
	        ysol[2] = -sqrt( R*R - (xsol[2] -xc[sectorcol] )*(xsol[2] -xc[sectorcol] )   )  + yc[sectorcol];
	      
	        ysol[3] = sqrt( R*R - (xsol[3] -xc[sectorcol] )*(xsol[3] -xc[sectorcol] )   )  + yc[sectorcol];
	      
	        ysol[4] = -sqrt( R*R - (xsol[4] -xc[sectorcol] )*(xsol[4] -xc[sectorcol] )   )  + yc[sectorcol];
	      
	        ysolr[1] = atraj*xsol[1] + btraj;
	        ysolr[2] = atraj*xsol[2] + btraj;
	        ysolr[3] = atraj*xsol[3] + btraj;
	        ysolr[4] = atraj*xsol[4] + btraj;
	      
	        dist[1] = sqrt(  (xsol[1]-xmin)*(xsol[1]-xmin) +  (ysol[1]-ymin)*(ysol[1]-ymin)   )	+  sqrt(  (xsol[1]-xmax)*(xsol[1]-xmax) +  (ysol[1]-ymax)*(ysol[1]-ymax)   );
	      
	        dist[2] = sqrt(  (xsol[2]-xmin)*(xsol[2]-xmin) +  (ysol[2]-ymin)*(ysol[2]-ymin)   )	+  sqrt(  (xsol[2]-xmax)*(xsol[2]-xmax) +  (ysol[2]-ymax)*(ysol[2]-ymax)   );
	      
	        dist[3] = sqrt(  (xsol[3]-xmin)*(xsol[3]-xmin) +  (ysol[3]-ymin)*(ysol[3]-ymin)   )	+  sqrt(  (xsol[3]-xmax)*(xsol[3]-xmax) +  (ysol[3]-ymax)*(ysol[3]-ymax)   );
	      
	        dist[4] = sqrt(  (xsol[4]-xmin)*(xsol[4]-xmin) +  (ysol[4]-ymin)*(ysol[4]-ymin)   )	+  sqrt(  (xsol[4]-xmax)*(xsol[4]-xmax) +  (ysol[4]-ymax)*(ysol[4]-ymax)   );
	      
	        if((fabs((xfinal-xsol[1]))<0.00000001 ) && (fabs((yfinal-ysol[1]))<0.00000001 )) dist[1] = 9999;
	        if((fabs((xfinal-xsol[2]))<0.00000001 ) && (fabs((yfinal-ysol[2]))<0.00000001 )) dist[2] = 9999;
	        if((fabs((xfinal-xsol[3]))<0.00000001 ) && (fabs((yfinal-ysol[3]))<0.00000001 )) dist[3] = 9999;
	        if((fabs((xfinal-xsol[4]))<0.00000001 ) && (fabs((yfinal-ysol[4]))<0.00000001 )) dist[4] = 9999;
	      
	        dmin = 999;
	      
	        if(fabs((ysolr[1]-ysol[1]))<0.00000001 ){
	      
	          if(dist[1]<dmin) {
	        
	            dmin = dist[1];
	            xfinal = xsol[1];
	            yfinal = ysol[1];
	          
	            deriy = - (xfinal-xc[sectorcol])/( sqrt(R*R  - ( xc[sectorcol] - xfinal )*( xc[sectorcol] - xfinal )  )  );
	            angy = atan(deriy);
	          
	            }
	          
	          }
	        
	        if(fabs((ysolr[2]-ysol[2]))<0.00000001 ){
	      
	          if(dist[2]<dmin){
	        
	            dmin = dist[2];
	            xfinal = xsol[2];
	            yfinal = ysol[2];
	            deriy =  (xfinal-xc[sectorcol])/( sqrt(R*R  - ( xc[sectorcol] - xfinal )*( xc[sectorcol] - xfinal )  )  );
	            angy = atan(deriy);
	          
	            }
	          
	          }
	        
	        if(fabs((ysolr[3]-ysol[3]))<0.00000001 ){	
		  
		  if(dist[3]<dmin){
		  
		    dmin = dist[3];
		    xfinal = xsol[3];
		    yfinal = ysol[3];
		    deriy = - (xfinal-xc[sectorcol])/( sqrt(R*R  - ( xc[sectorcol] - xfinal )*( xc[sectorcol] - xfinal )  )  );
		    angy = atan(deriy);
		    
		    }
		  
		  }
		
	        if(fabs((ysolr[4]-ysol[4]))<0.00000001 ){
	      
	          if(dist[4]<dmin){
	        
	            dmin = dist[4];
	            xfinal = xsol[4];
	            yfinal = ysol[4];
	            deriy =  (xfinal-xc[sectorcol])/( sqrt(R*R  - ( xc[sectorcol] - xfinal )*( xc[sectorcol] - xfinal )  )  );
	            angy = atan(deriy);
		  
		    }
		  
		  }
		
	        }
	      
	      vxfinal=vx*cos(2.0*angy)+vy*sin(2.0*angy);
	      vyfinal=vx*sin(2.0*angy)-vy*cos(2.0*angy);
	    
	      i = n*pointstotal +1;
	    
	      }
	    
	    j++;
	  
	    }
	
	  j = 0;  //Interaction with line sectors. Sectors 2, 5, 8, 11, ...
	  while(j<n){
	  
	    sectorcol = 2 + 3*j;
	    
	    if(sectorcol==sectoraux){
	    
	      ar[sectorcol] = coefang[ladaux];
	      br[sectorcol] = coelinear[ladaux];
	      
	      yfinal = ycollision;
	      xfinal = xcollision;
	      
	      deriy = ar[sectorcol];
	      angy = atan(deriy);
	      
	      vxfinal=vx*cos(2.0*angy)+vy*sin(2.0*angy);
	      vyfinal=vx*sin(2.0*angy)-vy*cos(2.0*angy); 
             
             i = n*pointstotal +1;
             
             }
             
           j++;
           
           }
           
         j = 0;  //Interaction with radius 'r' circunferences sectors. Sectors 3, 6, 9, 12, ...
         while(j<n){
         
           sectorcol = 3 + 3*j;
           
           if(sectorcol==sectoraux){
           
             if(sectorprev!=sectoraux){
             
               xsol[1] = (sqrt( -atraj*atraj*xc[sectorcol]*xc[sectorcol] + atraj*atraj*r*r - 2*atraj*btraj*xc[sectorcol] + 2*atraj*xc[sectorcol]*yc[sectorcol] - btraj*btraj + 2*btraj*yc[sectorcol] - yc[sectorcol]*yc[sectorcol] +r*r   )  -atraj*btraj  + atraj*yc[sectorcol]  + xc[sectorcol]  )/( atraj*atraj +1  );
               xsol[2] = xsol[1];
               
               xsol[3] = (-sqrt( -atraj*atraj*xc[sectorcol]*xc[sectorcol] + atraj*atraj*r*r - 2*atraj*btraj*xc[sectorcol] + 2*atraj*xc[sectorcol]*yc[sectorcol] - btraj*btraj + 2*btraj*yc[sectorcol] - yc[sectorcol]*yc[sectorcol] +r*r   )  -atraj*btraj  + atraj*yc[sectorcol]  + xc[sectorcol]  )/( atraj*atraj +1  );
               xsol[4] = xsol[3];
               
               ysol[1] = sqrt( r*r - (xsol[1] -xc[sectorcol] )*(xsol[1] -xc[sectorcol] )   )  + yc[sectorcol];
               
               ysol[2] = -sqrt( r*r - (xsol[2] -xc[sectorcol] )*(xsol[2] -xc[sectorcol] )   )  + yc[sectorcol];
               
               ysol[3] = sqrt( r*r - (xsol[3] -xc[sectorcol] )*(xsol[3] -xc[sectorcol] )   )  + yc[sectorcol];
               
               ysol[4] = -sqrt( r*r - (xsol[4] -xc[sectorcol] )*(xsol[4] -xc[sectorcol] )   )  + yc[sectorcol];
               
               ysolr[1] = atraj*xsol[1] + btraj;
               ysolr[2] = atraj*xsol[2] + btraj;
               ysolr[3] = atraj*xsol[3] + btraj;
               ysolr[4] = atraj*xsol[4] + btraj;
               
               dist[1] = sqrt(  (xsol[1]-xmin)*(xsol[1]-xmin) +  (ysol[1]-ymin)*(ysol[1]-ymin)   )	+  sqrt(  (xsol[1]-xmax)*(xsol[1]-xmax) +  (ysol[1]-ymax)*(ysol[1]-ymax)   );
               
               dist[2] = sqrt(  (xsol[2]-xmin)*(xsol[2]-xmin) +  (ysol[2]-ymin)*(ysol[2]-ymin)   )	+  sqrt(  (xsol[2]-xmax)*(xsol[2]-xmax) +  (ysol[2]-ymax)*(ysol[2]-ymax)   );
               
               dist[3] = sqrt(  (xsol[3]-xmin)*(xsol[3]-xmin) +  (ysol[3]-ymin)*(ysol[3]-ymin)   )	+  sqrt(  (xsol[3]-xmax)*(xsol[3]-xmax) +  (ysol[3]-ymax)*(ysol[3]-ymax)   );
               
               dist[4] = sqrt(  (xsol[4]-xmin)*(xsol[4]-xmin) +  (ysol[4]-ymin)*(ysol[4]-ymin)   )	+  sqrt(  (xsol[4]-xmax)*(xsol[4]-xmax) +  (ysol[4]-ymax)*(ysol[4]-ymax)   );
               
               dmin = 9999;
               
               if(fabs((ysolr[1]-ysol[1]))<0.00000001 ){
               
                 if(dist[1]<dmin){
                 
                   dmin = dist[1];
                   xfinal = xsol[1];
                   yfinal = ysol[1];
                   deriy = - (xfinal-xc[sectorcol])/( sqrt(r*r  - ( xc[sectorcol] - xfinal )*( xc[sectorcol] - xfinal )  )  );
                   angy = atan(deriy);
                   
                   }
                   
                 }
                 
               if(fabs((ysolr[2]-ysol[2]))<0.00000001 ){
               
                 if(dist[2]<dmin){
                 
                   dmin = dist[2];
                   xfinal = xsol[2];
                   yfinal = ysol[2];
                   deriy =  (xfinal-xc[sectorcol])/( sqrt(r*r  - ( xc[sectorcol] - xfinal )*( xc[sectorcol] - xfinal )  )  );
                   angy = atan(deriy);
                   
                   }
                   
                 }
                 
               if(fabs((ysolr[3]-ysol[3]))<0.00000001 ){
               
                 if(dist[3]<dmin){
                 
                   dmin = dist[3];
                   xfinal = xsol[3];
                   yfinal = ysol[3];
                   deriy = - (xfinal-xc[sectorcol])/( sqrt(r*r  - ( xc[sectorcol] - xfinal )*( xc[sectorcol] - xfinal )  )  );
                   angy = atan(deriy);
                   
                   }
                   
                 }
                 
               if(fabs((ysolr[4]-ysol[4]))<0.00000001 ){
               
                 if(dist[4]<dmin){
                 
                   dmin = dist[4];
                   xfinal = xsol[4];
                   yfinal = ysol[4];
                   deriy =  (xfinal-xc[sectorcol])/( sqrt(r*r  - ( xc[sectorcol] - xfinal )*( xc[sectorcol] - xfinal )  )  );
                   angy = atan(deriy);
                   
                   }
                   
                 }
                 
               }
               
             if(sectorprev==sectoraux){
             
               xsol[1] = (sqrt( -atraj*atraj*xc[sectorcol]*xc[sectorcol] + atraj*atraj*r*r - 2*atraj*btraj*xc[sectorcol] + 2*atraj*xc[sectorcol]*yc[sectorcol] - btraj*btraj + 2*btraj*yc[sectorcol] - yc[sectorcol]*yc[sectorcol] +r*r   )  -atraj*btraj  + atraj*yc[sectorcol]  + xc[sectorcol]  )/( atraj*atraj +1  );
               xsol[2] = xsol[1];
               
               xsol[3] = (-sqrt( -atraj*atraj*xc[sectorcol]*xc[sectorcol] + atraj*atraj*r*r - 2*atraj*btraj*xc[sectorcol] + 2*atraj*xc[sectorcol]*yc[sectorcol] - btraj*btraj + 2*btraj*yc[sectorcol] - yc[sectorcol]*yc[sectorcol] +r*r   )  -atraj*btraj  + atraj*yc[sectorcol]  + xc[sectorcol]  )/( atraj*atraj +1  );
               xsol[4] = xsol[3];
               
               ysol[1] = sqrt( r*r - (xsol[1] -xc[sectorcol] )*(xsol[1] -xc[sectorcol] )   )  + yc[sectorcol];
               
               ysol[2] = -sqrt( r*r - (xsol[2] -xc[sectorcol] )*(xsol[2] -xc[sectorcol] )   )  + yc[sectorcol];
               
               ysol[3] = sqrt( r*r - (xsol[3] -xc[sectorcol] )*(xsol[3] -xc[sectorcol] )   )  + yc[sectorcol];
               
               ysol[4] = -sqrt( r*r - (xsol[4] -xc[sectorcol] )*(xsol[4] -xc[sectorcol] )   )  + yc[sectorcol];
               
               ysolr[1] = atraj*xsol[1] + btraj;
               ysolr[2] = atraj*xsol[2] + btraj;
               ysolr[3] = atraj*xsol[3] + btraj;
               ysolr[4] = atraj*xsol[4] + btraj;
               
               dist[1] = sqrt(  (xsol[1]-xmin)*(xsol[1]-xmin) +  (ysol[1]-ymin)*(ysol[1]-ymin)   )	+  sqrt(  (xsol[1]-xmax)*(xsol[1]-xmax) +  (ysol[1]-ymax)*(ysol[1]-ymax)   );
               
               dist[2] = sqrt(  (xsol[2]-xmin)*(xsol[2]-xmin) +  (ysol[2]-ymin)*(ysol[2]-ymin)   )	+  sqrt(  (xsol[2]-xmax)*(xsol[2]-xmax) +  (ysol[2]-ymax)*(ysol[2]-ymax)   );
               
               dist[3] = sqrt(  (xsol[3]-xmin)*(xsol[3]-xmin) +  (ysol[3]-ymin)*(ysol[3]-ymin)   )	+  sqrt(  (xsol[3]-xmax)*(xsol[3]-xmax) +  (ysol[3]-ymax)*(ysol[3]-ymax)   );
               
               dist[4] = sqrt(  (xsol[4]-xmin)*(xsol[4]-xmin) +  (ysol[4]-ymin)*(ysol[4]-ymin)   )	+  sqrt(  (xsol[4]-xmax)*(xsol[4]-xmax) +  (ysol[4]-ymax)*(ysol[4]-ymax)   );
               
               if((fabs((xfinal-xsol[1]))<0.00000001 ) && (fabs((yfinal-ysol[1]))<0.00000001 )) dist[1] = 9999;
               if((fabs((xfinal-xsol[2]))<0.00000001 ) && (fabs((yfinal-ysol[2]))<0.00000001 )) dist[2] = 9999;
               if((fabs((xfinal-xsol[3]))<0.00000001 ) && (fabs((yfinal-ysol[3]))<0.00000001 )) dist[3] = 9999;
               if((fabs((xfinal-xsol[4]))<0.00000001 ) && (fabs((yfinal-ysol[4]))<0.00000001 )) dist[4] = 9999;
               
               dmin = 999;
               
               if(fabs((ysolr[1]-ysol[1]))<0.00000001 ){
               
                 if(dist[1]<dmin){
                 
                   dmin = dist[1];
                   xfinal = xsol[1];
                   yfinal = ysol[1];
                   
                   deriy = - (xfinal-xc[sectorcol])/( sqrt(r*r  - ( xc[sectorcol] - xfinal )*( xc[sectorcol] - xfinal )  )  );
                   angy = atan(deriy);
                   
                   }
                   
                 }
                 
               if(fabs((ysolr[2]-ysol[2]))<0.00000001 ){
               
                 if(dist[2]<dmin){
                 
                   dmin = dist[2];
                   xfinal = xsol[2];
                   yfinal = ysol[2];
                   deriy =  (xfinal-xc[sectorcol])/( sqrt(r*r  - ( xc[sectorcol] - xfinal )*( xc[sectorcol] - xfinal )  )  );
                   angy = atan(deriy);
                   
                   }
                   
                 }
                 
               if(fabs((ysolr[3]-ysol[3]))<0.00000001 ){
               
                 if(dist[3]<dmin){
                 
                   dmin = dist[3];
                   xfinal = xsol[3];
                   yfinal = ysol[3];
                   deriy = - (xfinal-xc[sectorcol])/( sqrt(r*r  - ( xc[sectorcol] - xfinal )*( xc[sectorcol] - xfinal )  )  );
                   angy = atan(deriy);
                   
                   }
                   
                 }
                 
               if(fabs((ysolr[4]-ysol[4]))<0.00000001 ){
               
                 if(dist[4]<dmin){
                 
                   dmin = dist[4];
                   xfinal = xsol[4];
                   yfinal = ysol[4];
                   deriy =  (xfinal-xc[sectorcol])/( sqrt(r*r  - ( xc[sectorcol] - xfinal )*( xc[sectorcol] - xfinal )  )  );
                   angy = atan(deriy);
                   
                   }
                   
                 }
                 
               }
               
             vxfinal=vx*cos(2.0*angy)+vy*sin(2.0*angy);
             vyfinal=vx*sin(2.0*angy)-vy*cos(2.0*angy); 
             
             i = n*pointstotal +1;
             
             }
             
           j++;
           
           }
           
         }
         
       i++;
       
       }
       
     sectorprev = sectoraux;
     
     atraj = (vyfinal/vxfinal);
     btraj = yfinal - atraj*xfinal;
     
     sum = 0;
     i = 1 ;
     while(i<(  n*pointstotal + 1 )  ){
     
       sum = sum + sqrt(  ( xdisc[i+1] - xdisc[i]  )*( xdisc[i+1] - xdisc[i]  ) + ( ydisc[i+1] - ydisc[i] )*( ydisc[i+1] - ydisc[i] )  );
	i++;
	
	}
	
      pertotal = sum;
      
      //Collision perimeter
      sum = 0;
      i = 1;
      while(i<ladaux){
      
        sum = sum + sqrt(  ( xdisc[i+1] - xdisc[i]  )*( xdisc[i+1] - xdisc[i]  ) + ( ydisc[i+1] - ydisc[i] )*( ydisc[i+1] - ydisc[i] )  );
        i++;
        
        }
        
      if(ladaux>1) per = sum +  sqrt(  ( xdisc[ladaux] - xfinal  )*( xdisc[ladaux] - xfinal   ) + ( ydisc[ladaux] - yfinal  )*( ydisc[ladaux] -  yfinal  )  );
      
      if(ladaux==1) per = sum +  sqrt(  ( xdisc[ladaux] - xfinal  )*( xdisc[ladaux] - xfinal   ) + ( ydisc[ladaux] - yfinal  )*( ydisc[ladaux] - yfinal   )  );
      
      per = per/pertotal;
      
      //Tangetial velocity
   
     if(ladaux<(ladmax+1)) vtan = vx*cos(ang[ladaux]+PI) + vy*sin(ang[ladaux]+PI);
     if(ladaux>(ladmax)) vtan = vx*cos(ang[ladaux]) + vy*sin(ang[ladaux]);
     
     if(cont>1){
     
       //Relative Measure
       ncvel = 1000; //Number of cells in the discretization of the tangential velocity
       
       varvel = -1;
       dvel = 2/ncvel;
       
       i = 1;
       while(i<(ncvel+1)){ //Selection of the cell to store a given velocity
       
         if((vtan>(varvel+dvel*(i-1))) && (vtan<(varvel+dvel*i))){
         
           ifinal = i;
           icont[cont] = i;
           i = ncvel;
           
           }
           
         i++;
         
         }
         
       ncper = 1000; //Number of cells in the discretization of the perimeter
       
       varper = 0;
       
       dper = 1/ncper; // Total perimeter normalized to 1;
       
       j = 1;
       while(j<(ncper+1)){ //Selection of the cell to store a given perimeter
       
         if((per>(varper+dper*(j-1))) && (per<(varper+dper*j))){
       
           jfinal = j;
           jcont[cont] = j;
           j = ncper;
         
           }
           
         j++;
         
         }
         
       }
       
     if((cont%1000000)==0) printf("%d %d %d %.12f %.12f %.12f %.12f\n",IC,cont,ladaux,xfinal,yfinal,vxfinal,vyfinal,per);
     
     xvet[cont] = xfinal;
     yvet[cont] = yfinal;
     
     vx = vxfinal;
     vy = vyfinal;
     
     cont++;
     
     }  //End of a IC
     
   cont = 1;
   contstick = 1;
   count = 1;
   
   while(cont<(timetraj+1)){  //Identification of occupied boxes when the trajectory presents a minimum of 20 collisions between the straight sectors of the billiard
   
     if(((vetsector[cont])%3)==2){
     
       if(cont==(contstick+1)){
       
         count++;
         if(count>20){
         
           contfinal=count; //Time selected at a trapping moment
           conttrajfinal = cont;
           
           }
           
         }
         
       else{
       
         if(count>20){
         
           i = conttrajfinal - contfinal +1;
           while(i<(conttrajfinal+1)){
           
             if(box[jcont[i]][icont[i]] ==0){
             
               fillingbox = fillingbox + 1;
               box[jcont[i]][icont[i]] =  box[jcont[i]][icont[i]] + 1;
               
               }
               
             if(i>(conttrajfinal - contfinal +1)){
             
               distreg = sqrt(  (xvet[i]-xvet[i-1])*(xvet[i]-xvet[i-1]) +  (yvet[i]-yvet[i-1])*(yvet[i]-yvet[i-1])    );
               sumreg = sumreg + distreg;
               contreg = contreg + 1;
               
               }
               
             i++;
             
             }
             
           }
           
         count = 1;
         
         }
         
       contstick = cont;
       
       }
       
     if(cont>1){
     
       disttotal = sqrt(  (xvet[cont]-xvet[cont-1])*(xvet[cont]-xvet[cont-1]) +  (yvet[cont]-yvet[cont-1])*(yvet[cont]-yvet[cont-1])    );
       sumtotal = sumtotal + disttotal;
       
       }
       
     cont++;
     
     }
     
   qsi[IC] = fillingbox/1000000; //Occupied boxes at trapping moments for each initial condition
   
   sumqsi = sumqsi + qsi[IC];
   
   contchaos = (timetraj-1) - contreg;
   
   sumchaos = sumtotal - sumreg;
   
   printf("cont %d \n",contchaos);
   printf("somacaos %f \n",sumchaos);
   
   medchaos = sumchaos/contchaos;
   medreg = sumreg/contreg;
   
   dchaoticfinal[IC] = medchaos;
   sumchaosfinal = dchaoticfinal[IC] + sumchaosfinal;
   
   dregfinal[IC] = medreg;
   sumregfinal = sumregfinal + dregfinal[IC];
   
   printf("medchaos %f \n",medchaos);
   printf("medreg %f \n",medreg);
   printf("rho %f \n",rho);
   
   fprintf(output3, "%d %.12f %.12f %.12f %.12f %.12f %.12f %.12f \n",IC,fillingbox,qsi[IC],sumqsi,medchaos,sumchaosfinal,medreg,sumregfinal);
   
   IC++;
   
   }  //End of all ICs

  mediafinalqsi = sumqsi/ICs; //Average number of occupied boxes for different initial conditions that exhibited trapping
  mediafinalreg =  sumregfinal/ICs;
  mediafinalcao =  sumchaosfinal/ICs;

  IC = 1;
  sumdeviation = 0;
  while(IC<=ICs){
  
    deviation = qsi[IC] - mediafinalqsi;
    
    sumdeviation = sumdeviation + deviation*deviation;
    IC++;
    
    }

  deviationqsi = sqrt(sumdeviation/(ICs-1));
  
  IC = 1;
  sumdeviation = 0;
  while(IC<=ICs){
  
    deviation = dchaoticfinal[IC] - mediafinalcao;
    
    sumdeviation = sumdeviation + deviation*deviation;
    IC++;
    
    }

  deviationcao = sqrt(sumdeviation/(ICs-1));

  IC = 1;
  sumdeviation = 0;
  while(IC<=ICs){
  
    deviation = dregfinal[IC] - mediafinalreg;
    
    sumdeviation = sumdeviation + deviation*deviation;
    IC++;
    
    }

  deviationreg = sqrt(sumdeviation/(ICs-1));

  rho = mediafinalqsi/( mediafinalqsi + (1-mediafinalqsi)*(mediafinalreg/mediafinalcao)    );

  fprintf(output3, "%.12f %.12f %.12f %.12f %.12f %.12f %.12f\n",rho,mediafinalqsi,deviationqsi,mediafinalcao,deviationcao,mediafinalreg,deviationreg);
		
  fclose(output);
  fclose(output3);
  fclose(input1);
  fclose(input2);
  
  }
