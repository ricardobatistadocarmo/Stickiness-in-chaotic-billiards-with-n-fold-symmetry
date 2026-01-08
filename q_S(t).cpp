/*******************************
 *
 *  Author: Ricardo B. do Carmo
 *  Affiliation: Instituto Federal de Alagoas
 *  Email: ricardo.carmo@ifal.edu.br
 *
 *  Description:
 *  This code performs numerical calculations to obtain phase-space
 *  points related to trapping events along the same trajectory.
 *
 *  The program was developed for scientific and academic purposes, and is intended
 *  to support research in nonlinear dynamics, chaos, and billiard-like systems.
 *
 *  If you use or modify this code for academic work, please cite or acknowledge
 *  the author appropriately.
 *
 *  Created: [June 2024]
 *  Last modified: [July 2025]
 *
 *******************************/

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
int i,j,n,sectors,lad, pointsR,pointsline,pointsr,pontostotal,ladaux,resto,divi,sectorcol,sectoraux,cont,ok,sectorprev,ladmax,k;
double aloop[10000], bloop[10000],a,b,xv[100000],yv[100000],xc[100000],yc[100000],gam1,gam2,alf,t4,t3,t2,t1,dx;
double  dtR,dtr,xrot,yrot;
double dify,difx;
double vx, vy, ytraj,xtraj, atraj,btraj, xtest, ytest, xcollision, ycollision;
double x1,x2,xfinal,yfinal,ar[100000],br[100000],xsol[5],ysol[5],ysolr[5],dist[5],dmin;
double y11, y22,deriy,angy, vxfinal,vyfinal,xmin,xmax,ymin,ymax,sum,per,pertotal,vtan,angmax;
double ncvel, varvel, dvel,ncper,varper,dper,medrel,box,cont_med,xi,yi;
int box_matrix[3001][3001],ifinal,jfinal,ICs,IC,IC_cont,temptraj,conttrajfinal,contfinal,contstick,count ;

double *xdisc, *ydisc, *coefang, *coelinear, *ang, *sum_rm ;
int *icont, *jcont, *vetsetor;

int pborder;

int main(){

  pborder = 50000001;
  
  xdisc = (double *) malloc(pborder * sizeof(double));
  ydisc = (double *) malloc(pborder * sizeof(double));
  coefang = (double *) malloc(pborder * sizeof(double));
  coelinear = (double *) malloc(pborder * sizeof(double));
  ang = (double *) malloc(pborder * sizeof(double));
  sum_rm = (double *) malloc(pborder * sizeof(double));
  icont = (int *) malloc(pborder * sizeof(int));
  jcont = (int *) malloc(pborder * sizeof(int));
  vetsetor = (int *) malloc(pborder * sizeof(int));

  myrand=1314;  //Seed for random numbers. Any integer number.
  Init_Random();

  temptraj = 50000000; //Maximum number of collisions
  ICs = 1;

  PI=3.1415926535897932384626433832795028841971693993751058209749445923;
  sqrt3=sqrt(3);	

  R = 0.5;  //Radius of one of the circunference that make up the boundary
  r = 0.25; // Radius of the second circunference
  n = 3; // Billiard symmetry
  sectors = 3*n;  //lad from 1 to 3*n

  xc[1] = 1-R;  //before xc5
  yc[1] = 0;  //before yc5

  xc[sectors] = 1 - r;  //Before xc3. Sectors is the last lad in counter-clockwise
  yc[sectors] = 0;  //before yc3

  xc[3] = xc[sectors]*cos(2*PI/n) - yc[sectors]*sin(2*PI/n); 
  yc[3] = xc[sectors]*sin(2*PI/n) + yc[sectors]*cos(2*PI/n); 

  printf("%d  \n",sectors);	

  FILE* output;
  output = fopen("c3_1CI_stick_20_3000x3000_50mi.dat","w");	

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
  	
  //Three basal points Cálculo dos três primeiros pontos:
  gam1 = ( (3/(r*r))  + (1/(r*r))   )*a - ( -(sqrt3/(r*r))  + (sqrt3/(r*r))    );
  gam2 = ( (3/(r*r))  + (1/(r*r))   )*(b-yc[3]) + ( -(sqrt3/(r*r))  + (sqrt3/(r*r)))*xc[3];
  alf = 4/(r*r*r*r);

  xv[1] = 1;
  yv[1] = 0; 

  xv[2] = -(a*b - xc[1]*pow(R/R,2) )/( a*a + pow(R/R,2)   );
  yv[2] =  sqrt( R*R - pow((xv[2]-xc[1])*(R/R),2)  )  ;

  xv[3] =  - (gam1*gam2 - 4*xc[3]*alf)/(gam1*gam1 + 4*alf);
  yv[3] =    (( - ((sqrt3*(xv[3]-xc[3]))/(r*r))  + (sqrt3/(r*r))*(xv[3]-xc[3]) + 2*sqrt( (3/(r*r)) + (1/(r*r)) - (4/(r*r*r*r))*pow(xv[3]-xc[3],2)    ))/( (3/(r*r)) + (1/(r*r))  ))   + yc[3]  ;

  //Points between sectors' borders
  i = 4;
  while(i<(sectors+1)){
  
    xv[i] = xv[i-3]*cos(2*PI/n) - yv[i-3]*sin(2*PI/n); 
    yv[i] = xv[i-3]*sin(2*PI/n) + yv[i-3]*cos(2*PI/n);
    
    i++;
    
    }
	
  xv[sectors+1] = xv[1];
  yv[sectors+ 1] = yv[1];	
	
  //Discretization of sectors' borders 1, 2, and 3
  t4 = fabs(acos((xv[4]-xc[3])/r));

  t3 = fabs(acos((xv[3]-xc[3])/r));

  t2 = fabs(acos((xv[2]-xc[1])/R));

  t1 =0;

  printf(" tempos %.12f %.12f %.12f %.12f \n",t4, t3, t2,t1);

  pointsR = 400 ; //We use 100 for almost all symmetries
  dtR = (t1-t2)/pointsR;

  i=1;
  xdisc[i] = xv[1];
  ydisc[i] = yv[1];
  while(i<(pointsR+1)){
  
    xdisc[i+1] = xc[1]  + R*cos(t1 - dtR*i);
    ydisc[i+1] = sqrt( R*R - (xdisc[i+1] - xc[1] )*(xdisc[i+1] - xc[1] )   ) +yc[1] ;
    
    i++;
    
    }	

  //Discretization od sector 2 Discretização do setor 2
  pointsline = 400;
  dx = (xv[2] - xv[3])/pointsline;

  i = 1;
  xdisc[pointsR + i] = xv[2];
  ydisc[pointsR + i] = yv[2];
  while(i<(pointsline+1)){
	
    xdisc[i+pointsR+1] = xv[2] - dx*i;
    ydisc[i+pointsR+1] = a*xdisc[i+pointsR+1] + b;
    
    i++;	
    
    }

  //Discretization of sector 3
  pointsr = 400;
  dtr = (t3-t4)/pointsr;

  i=1;
  xdisc[i+ pointsR+ pointsline] = xv[3];
  ydisc[i+ pointsR+ pointsline] = yv[3];
  while(i<(pointsr+1)){
  
    xdisc[i+pointsR+pointsline+1] = xc[3]  + r*cos(t3 - dtr*i);
    ydisc[i+pointsR+pointsline+1] = sqrt( r*r - (xdisc[i+pointsR+pointsline+1] - xc[3] )*(xdisc[i+pointsR+pointsline+1] - xc[3] )   ) +yc[3] ;
    
    i++;
    
    }	

  pontostotal = pointsR+pointsline+pointsr;	
	
  //Rotation
  j = 1;
  while(j<n){
  
    i=1;
    while(i<(pointsR + pointsr + pointsline +1) ){
    
      xrot = xdisc[i]*cos(2*PI*j/n) - ydisc[i]*sin(2*PI*j/n);
      yrot = xdisc[i]*sin(2*PI*j/n) + ydisc[i]*cos(2*PI*j/n);
      
      xdisc[(pointsR + pointsr + pointsline)*j + i] = xrot;
      ydisc[(pointsR + pointsr + pointsline)*j + i] = yrot;
      
      i++;
      
      }
      
    j++;
    
    }
    
  xdisc[n*pontostotal +1] = xdisc[1];
  ydisc[n*pontostotal +1] = ydisc[1];	

  i = 1;	
  while(i<(  n*pontostotal + 1 )  ){
  
    coefang[i] = 0;
    coelinear[i] = 0;
    ang[i] = 0;
    
    i++;
    
    }
    
  i = 1;
  angmax=0;
  while(i<(  n*pontostotal + 1 )  ){
  
    dify = (ydisc[i+1] - ydisc[i]);
    difx = (xdisc[i+1] - xdisc[i] );
    coefang[i] = dify/difx;
    
    coelinear[i] = ydisc[i] - coefang[i]*xdisc[i];
    ang[i] = atan(coefang[i]); 

    if((ang[i]>angmax) && (i<(pontostotal*n))){
    
      angmax = ang[i];
      ladmax = i;
      
      }
      
    if(i==(n*pontostotal)){
    
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
    
    i++
    
    }
    
  //Radius 'r' circunferences centers
  i = 1;
  while(i<(n-1)){
  
    xc[3+3*i] = xc[3*i]*cos(2*PI/n) - yc[3*i]*sin(2*PI/n); 
    yc[3+3*i] = xc[3*i]*sin(2*PI/n) + yc[3*i]*cos(2*PI/n);
    
    i++;
    
    }
    	
  i = 0;
  while(i<(temptraj+1)){
  
    sum_rm[i] = 0;
    i++;
    
    }
    
  IC = 1;	
  while(IC<(ICs+1)){
  
    box = 0;
    i=1;
    j=1;
    
    while(i<3001){
    
      j=1;
      while(j<3001){
      
        box_matrix[i][j]=0;
        j++;
        
	}
	
      i++;
      
      }
    
    //Initial Conditions
    xi = FRANDOM*0.9999;
    yi = 0;
			
    vx = 2*FRANDOM -1;
    vy = sqrt(1-vx*vx); //Must be positive

    atraj =  (vy/vx);
    btraj = -atraj*xi;

    sectoraux = 0;
    cont = 1;
    while(cont<(temptraj+1)){
    
      ok = 0;
      i = 1;
      while(i<(n*pontostotal +1)){
      
        if(i==ladaux) i++
	
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
	  
	  resto = ladaux % 400; // The number 400 is the same one used in the discretization of each sector.
	  divi = ladaux / 400; // The number 100 is the same one used in the discretization of each sector.
	  
	  if(resto==0) sectoraux = divi;
	  if(resto!=0) sectoraux = divi + 1;
	  
	  j = 0;  //Interaction between radius 'R' circunferences sectors. Sectors 1, 4, 7, 10, ...
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
              
            vxfinal=vx*cos(2.0*angy)+vy*sin(2.0*angy);
            vyfinal=vx*sin(2.0*angy)-vy*cos(2.0*angy);
            
            i = n*pontostotal +1;
            
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
            
            i = n*pontostotal +1;
            
            }
            
          j++;
          
          }
          
        j = 0; //Interaction with radius 'r' circunferences sectors. Sectors 3, 6, 9, 12, ...
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
            
            i = n*pontostotal +1;
            
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
    while(i<(  n*pontostotal + 1 )  ){
    
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
    
    //Tangential velocity
    if(ladaux<(ladmax+1)) vtan = vx*cos(ang[ladaux]+PI) + vy*sin(ang[ladaux]+PI);
    if(ladaux>(ladmax)) vtan = vx*cos(ang[ladaux]) + vy*sin(ang[ladaux]);
    
    if(cont>1){
    
      //Relative Measure
      ncvel = 3000; //Number of cells in the discretization of the tangential velocity
      
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
        
      ncper = 3000; //Number of cells in the discretization of the perimeter
      
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
        
      vetsetor[cont] = sectoraux;
      
      }
    
    if((cont%100000)==0) printf("%d %d %d %.12f %.12f %.12f %.12f\n",IC,cont,ladaux,xfinal,yfinal,vxfinal,vyfinal,per);
  
    vx = vxfinal;
    vy = vyfinal;
  
    cont = cont + 1;
  
    } //End of a IC
    
  IC++;
    
  }  //End of all ICs

  cont = 1;
  contstick = 1;
  count = 1;
  while(cont<(temptraj+1)){ //Identification of occupied boxes when the trajectory presents a minimum of 20 collisions between the straight sectors of the billiard
  
    if(((vetsetor[cont])%3)==2){
    
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
          
            if(box_matrix[jcont[i]][icont[i]] ==0){
            
              box++;
              box_matrix[jcont[i]][icont[i]] = box_matrix[jcont[i]][icont[i]] + 1;
              
              fprintf(output, "%d %d %d %d %f\n",conttrajfinal,contfinal,jcont[i],icont[i],box);
              
              }
              
            i++;
            
            }
            
          }
          
        count = 1;
        
        }
        
      contstick = cont;
      
    }
    
  cont++;	

  }
		
  fclose(output);	
  fclose(input1);
  fclose(input2);

  }
