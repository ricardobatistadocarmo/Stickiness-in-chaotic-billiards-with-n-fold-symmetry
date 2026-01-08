/*******************************
 *
 *  Author: Ricardo B. do Carmo
 *  Affiliation: Instituto Federal de Alagoas
 *  Email: ricardo.carmo@ifal.edu.br
 *
 *  Description:
 *  This code numerically calculates the perimeter of symmetric billiards.
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

double R,r,PI,sqrt3;
int i,j,n,sectors,lad, pointsR,pointsline,pointsr,pointstotal,ladaux,rest,divi,sectorcol,sectoraux,cont,ok,sectorprev,ladmax;
double aloop[10000], bloop[10000],a,b,xv[100000],yv[100000],xc[100000],yc[100000],gam1,gam2,alf,t4,t3,t2,t1,dx,per_line,per_circ,per_line_nor;
double  dtR,dtr,xrot,yrot;
double dify,difx;
double pertotal,angmax;

double *xdisc, *ydisc, *coefang, *coelinear, *ang;
int pborder;

int main(){

  pborder = 50000001;
  
  xdisc = (double *) malloc(pborder * sizeof(double));
  ydisc = (double *) malloc(pborder * sizeof(double));
  coefang = (double *) malloc(pborder * sizeof(double));
  coelinear = (double *) malloc(pborder * sizeof(double));
  ang = (double *) malloc(pborder * sizeof(double));

  FILE* output;
  output = fopen("c2_a_c1000.dat","w");	

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

  PI=3.1415926535897932384626433832795028841971693993751058209749445923;
  sqrt3=sqrt(3);	

  R = 0.9;  //Radius of one of the circunference that make up the boundary
  r = 0.45; // Radius of the second circunference
  n = 2; // Billiard symmetry

  while(n<1001){
  
    sectors = 3*n;  //lad from 1 to 3*n
    
    xc[1] = 1-R;  //before xc5
    yc[1] = 0;  //before yc5
    
    xc[sectors] = 1 - r;  //Before xc3. Sectors is the last lad in counter-clockwise
    yc[sectors] = 0;  //antes yc3
    
    xc[3] = xc[sectors]*cos(2*PI/n) - yc[sectors]*sin(2*PI/n);
    yc[3] = xc[sectors]*sin(2*PI/n) + yc[sectors]*cos(2*PI/n);
    
    printf("%d  \n",sectors);
    
    a = aloop[n];
    b = bloop[n];
    
    printf("%d %.12f %.12f \n",n,a,b);
    
    //Three basal points:
    
    gam1 = ( (3/(r*r))  + (1/(r*r))   )*a - ( -(sqrt3/(r*r))  + (sqrt3/(r*r))    );
    gam2 = ( (3/(r*r))  + (1/(r*r))   )*(b-yc[3]) + ( -(sqrt3/(r*r))  + (sqrt3/(r*r)))*xc[3];
    alf = 4/(r*r*r*r);
    
    xv[1] = 1;
    yv[1] = 0;
    
    xv[2] = -(a*b - xc[1]*pow(R/R,2) )/( a*a + pow(R/R,2)   );
    yv[2] =  sqrt( R*R - pow((xv[2]-xc[1])*(R/R),2)  )  ;
    
    xv[3] =  - (gam1*gam2 - 4*xc[3]*alf)/(gam1*gam1 + 4*alf);
    yv[3] =    (( - ((sqrt3*(xv[3]-xc[3]))/(r*r))  + (sqrt3/(r*r))*(xv[3]-xc[3]) + 2*sqrt( (3/(r*r)) + (1/(r*r)) - (4/(r*r*r*r))*pow(xv[3]-xc[3],2)    ))/( (3/(r*r)) + (1/(r*r))  ))   + yc[3]  ;

    i = 4; //Points between sectors' border
    while(i<(sectors+1)){
    
      xv[i] = xv[i-3]*cos(2*PI/n) - yv[i-3]*sin(2*PI/n); 
      yv[i] = xv[i-3]*sin(2*PI/n) + yv[i-3]*cos(2*PI/n);
      
      i++;
      
      }
      
    xv[sectors+1] = xv[1];
    yv[sectors+ 1] = yv[1];
    
    //Discretization of sectors' border 1, 2, and, 3
    t4 = fabs(acos((xv[4]-xc[3])/r));
    
    t3 = fabs(acos((xv[3]-xc[3])/r));
    
    t2 = fabs(acos((xv[2]-xc[1])/R));
    
    t1 = 0;
    
    pointsR = 100;
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
    pointsline = 100;
    dx = (xv[2] - xv[3])/pointsline;
    i =  1;
    
    xdisc[pointsR + i] = xv[2];
    ydisc[pointsR + i] = yv[2];
    
    while(i<(pointsline+1)){
    
      xdisc[i+pointsR+1] = xv[2] - dx*i;
      ydisc[i+pointsR+1] = a*xdisc[i+pointsR+1] + b;
      
      i++;
      
      }
      
    //Discretization of sector 3
    pointsr = 100;
    dtr = (t3-t4)/pointsr;
    
    i=1;
    xdisc[i+ pointsR+ pointsline] = xv[3];
    ydisc[i+ pointsR+ pointsline] = yv[3];
    
    while(i<(pointsr+1)){
    
      xdisc[i+pointsR+pointsline+1] = xc[3]  + r*cos(t3 - dtr*i);
      ydisc[i+pointsR+pointsline+1] = sqrt( r*r - (xdisc[i+pointsR+pointsline+1] - xc[3] )*(xdisc[i+pointsR+pointsline+1] - xc[3] )   ) +yc[3];
      
      i++;
      
      }
      
    pointstotal = pointsR+pointsline+pointsr;	
	
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

    xdisc[n*pointstotal +1] = xdisc[1];
    ydisc[n*pointstotal +1] = ydisc[1];

    i = 1;	
    while(i<(  n*pointstotal + 1 )  ){
  
      coefang[i] = 0;coelinear[i] = 0;
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
    
    //Total perimeter
    soma = 0;
    i = 1;
    while(i<(  n*pointstotal + 1 )  ){
  
      soma = soma + sqrt(  ( xdisc[i+1] - xdisc[i]  )*( xdisc[i+1] - xdisc[i]  ) + ( ydisc[i+1] - ydisc[i] )*( ydisc[i+1] - ydisc[i] )  );
      i++;
    
      }
    
    pertotal = soma; 

    per_line =  (n*sqrt(  ( xv[3] - xv[2]  )*( xv[3] - xv[2]  ) + ( yv[3] - yv[2] )*( yv[3] - yv[2] )  )); //Perimeter of the straight sectors

    per_circ = pertotal - per_line; //Perimeter of the circular sectors

    per_line_nor = per_line/pertotal;

    fprintf(output,"%d %.12f %.12f %.12f %.12f\n",n,per_line,per_circ,pertotal,per_line_nor);
    printf("%d %.12f %.12f %.12f %.12f \n",n,per_line,per_circ,pertotal,per_line_nor);  

    n++;	

    }
			
  fclose(output);	
  fclose(input1);
  fclose(input2);

  }
