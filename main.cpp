#include <iostream>
#include <stdio.h>      
#include <stdlib.h>     
#include <time.h>       
#include <math.h>
#include <cmath>
#include <fstream>
#include <vector>
#include <gsl/gsl_linalg.h>


//parameters
#define PI 3.14159265
#define n 15                               //number of cell
#define r 10                               //radius of each cell (micro-meter/sec)
#define sizeX 120                          //length of tissue
#define sizeY 100                          //width of tissue (porrosity α = 1 - ((n*PI*r^2)/(sizeX*sizeY)))
#define probabilityOfGettingCapture 0.1    //probability of uptake of nano particles (in paper denoted by 'ρ')
#define radiOfNanoParticle 100             //radius of nano particle in nanometer (in paper denoted by 'a', unit nano-meter)
#define bulkDelt 0.1                       //large time step in bulk (sec)
#define clearance 0.25                     //clearance around the cell to avoid overlap (rejniak model)
#define numb 400000000                     //maximum number of discrete points
#define tht 2                              //angle in cell to create discrete points
#define LBP 240                            //number of lower boundary discrete points
#define UBP 240                            //number of Upper boundary discrete points
#define IPF 200                            //number of input discrete points 



using namespace std;

double* x =                (double*)malloc(n*sizeof(double));        // x coordinate of center of cell
double* y =                (double*)malloc(n*sizeof(double));        // y coordinate of center of cell
double* aPointx =          (double*)malloc(numb*sizeof(double));     // x coordinate of the points
double* aPointy =          (double*)malloc(numb*sizeof(double));     // y coordinate of the point
double* alpha =            (double*)malloc(numb*sizeof(double));  
double* beta =             (double*)malloc(numb*sizeof(double)); 
double* ak =               (double*)malloc(numb*sizeof(double)); 
double* bk =               (double*)malloc(numb*sizeof(double));  
double* rowValue =         (double*)malloc(numb*sizeof(double));          //matrix
double* constantOfMatrix = (double*)malloc(numb*sizeof(double));          //constant
double* forceValue =       (double*)malloc(numb*sizeof(double)); 
double* sak =              (double*)malloc(numb*sizeof(double)); 
double* sbk =              (double*)malloc(numb*sizeof(double));  
double* salpha =           (double*)malloc(numb*sizeof(double));  
double* sbeta =            (double*)malloc(numb*sizeof(double));
double* cbeta =            (double*)malloc(numb*sizeof(double));


//Global variables
int count=0, count1=0, count2=0, count3=0, constantOfMatrixCounter=0;         //count is discrete point counter
double t=0, ypx1=0, ypy1=0, ypx=0, ypy=0, ypx2=0, ypy2=0, timechecker=10;     //t is time and ypx are positions at different time

//Initializing the function
void pointCounter();
void cell();
void discretePoints();
void linearequtionMatrix();
void linearEquationSolver();
void mainSimulation();
void removeAll();


int main ()
{
	
//remove previous data
removeAll();
     	
//count total point in the tissue	
pointCounter();


// initialize random seed:
srand (time(NULL));

//random cell position
cell();


//discretizing the whole system
 discretePoints();


//calculating the value of each element the matrix to solve linear set of equation 
linearequtionMatrix();


//solving linear set of equation using GSL package
linearEquationSolver();

//main simulation
mainSimulation();

  return 0;
}





//functions

//function for remove previous data
void removeAll()
{
	remove("cell_center_positions.txt");
	remove( "particle_capture_cell_identity.txt " );
	remove("particle1_travelpath.xlsx");
	remove("particle2_travelpath.xlsx");
	remove("particle3_travelpath.xlsx");
	remove("particle4_travelpath.xlsx");
	remove("particle5_travelpath.xlsx");
}


//function for count total point in the tissue
void pointCounter()
{
    count  = (360/tht)*n+LBP+UBP+IPF;
    count1 = (360/tht)*n;
    count2 = (360/tht)*n+LBP;
    count3 = (360/tht)*n+LBP+UBP;	
}


//function for random cell position
void cell()
{
	 /* generate random center point of the cell: */
	 ofstream myfile;
     myfile.open ("cell_center_positions.txt");
     ofstream out("cell_center_positions.txt", ios::app);
     srand (time(NULL));
     
     for(int i=0; i<50; i++){
         x[0] = (double)rand()/(RAND_MAX + 0)+0+(rand()%sizeX);
         y[0] = (double)rand()/(RAND_MAX + 0)+0+(rand()%sizeY);
            if(x[0]>=(r+clearance) && x[0]<=(sizeX-(r+clearance)) && y[0]>=(r+clearance) && y[0]<=(sizeY-(r+clearance)))     //is used to make sure no half cell
	       break;
     }
     
     out<< "centre is:  ";
	 out<< x[0] <<"    ";
     out<< y[0] << "\n\n";
      
     for(int i=0; i<n; i++){
	     int l=0, p=0;
	     x[i] = (double)rand()/(RAND_MAX + 0)+0+(rand()%sizeX);
         y[i] = (double)rand()/(RAND_MAX + 0)+0+(rand()%sizeY);
	     for(int j=0; j<=i; j++){
		     l= ((x[i]-x[j])*(x[i]-x[j])) + ((y[i]-y[j])*(y[i]-y[j])); // center to center distacnce = sum of radii
		     p=j;
		     if(l <= (4*(r + clearance)* (r + clearance)))
			    break;
		  if(x[i]<=(r+clearance) || x[i]>=(sizeX-(r+clearance)) || y[i]<=(r+clearance) || y[i]>=(sizeY-(r+clearance)))         //is used to make sure no half cell
			  break;
	      }
	      
	      if(i==p){
		  
		     out<<"\n\ncentre is:  ";
	         out<< x[i] <<"    ";
             out<< y[i] << "\n\n";
         
	     }
	     else{
		      i= i-1;
	     }	  
     }      
     myfile.close();     	
}


//function for discretizing the whole system
void discretePoints()
{     	
	 int pointx=0, pointy=0;
     for(double theta=tht; theta<=360; theta=theta+tht){
	     aPointx[pointx]= double(x[0]) + double((r*cos ( theta* PI / 180.0 )));
	     aPointy[pointy]= double(y[0]) + double((r*sin ( theta* PI / 180.0 )));
	     pointx=pointx+1;
	     pointy=pointy+1;
     }

     for(int i=1; i<n; i++){
		 for(double theta=tht; theta<=360; theta=theta+tht){
	         aPointx[pointx]= double(x[i]) + double((r*cos ( theta* PI / 180.0 )));
	         aPointy[pointy]= double(y[i]) + double((r*sin ( theta* PI / 180.0 )));
	         pointx=pointx+1;
	         pointy=pointy+1;

         }
	  }

     double sizep, sizeq,lbp,ipf;
     sizep=sizeX-1;
     sizeq=sizeY-2;
     lbp=LBP;
     ipf=IPF;
     double increment1 = sizep/lbp;
     double increment2 = sizeq/ipf;

     double forceLowerBoundaryX=0.5;
     for(int forceBL=count1; forceBL<count2; forceBL++){
	     aPointx[forceBL]=forceLowerBoundaryX;
	     aPointy[forceBL]=0;
	     forceLowerBoundaryX=forceLowerBoundaryX+increment1;
     }  
 

     double forceUpperBoundaryX=0.5;
     for(int forceBU=count2; forceBU<count3; forceBU++){
	     aPointx[forceBU]=forceUpperBoundaryX;
	     aPointy[forceBU]=sizeY;
	     forceUpperBoundaryX=forceUpperBoundaryX+increment1;
     }  


     double forceIny=1;
     for(int forceIn=count3; forceIn<count; forceIn++){
	     aPointx[forceIn]=0;
	     aPointy[forceIn]=forceIny;
	     forceIny=forceIny+increment2;
     }
}


//function for calculating the value of each element the matrix to solve linear set of equation 
void linearequtionMatrix()
{
     //force calculation     
     double rk=0, rk1=0, aPointx1=0, aPointy1=0; 
     int rowValueCounter=0, constantTerm = 0, constantTerm2 = 0, constantTerm3 = 0;
     for(int distance=0; distance<count; distance++){
	     aPointx1=aPointx[distance];
	     aPointy1=aPointy[distance];

         for(int distance1=0; distance1<count; distance1++){
	         ak[constantTerm] = aPointx1-aPointx[distance1];
	         bk[constantTerm] = aPointy1-aPointy[distance1];
	         rk1 = ak[constantTerm]*ak[constantTerm]+bk[constantTerm]*bk[constantTerm];
	         rk = sqrt(rk1);
	         alpha[constantTerm] = (1/(4*PI*(2.5)))*(-(0.5*(log((rk*rk)+(0.25))))+((0.25)/((rk*rk)+0.25)));
	         beta[constantTerm] = (1/(4*PI*(2.5)))*(1/((rk*rk)+(0.25)));
	         constantTerm = constantTerm + 1;
	 
         }
	 
	     for(int distance2=0; distance2<count; distance2++){
		     rowValue[rowValueCounter] = alpha[constantTerm2] + beta[constantTerm2]*ak[constantTerm2]*ak[constantTerm2];
		     rowValueCounter = rowValueCounter + 1;
		     rowValue[rowValueCounter] = beta[constantTerm2]*ak[constantTerm2]*bk[constantTerm2];
		     rowValueCounter = rowValueCounter + 1;
		     constantTerm2 = constantTerm2 + 1;
	     }
	     if(distance<count3){
		  constantOfMatrix[constantOfMatrixCounter] = 0;
		  constantOfMatrixCounter = constantOfMatrixCounter + 1;
	     }
	     for(int distance3=0; distance3<count; distance3++){
		     rowValue[rowValueCounter] = beta[constantTerm3]*ak[constantTerm3]*bk[constantTerm3];
		     rowValueCounter = rowValueCounter + 1;
		     rowValue[rowValueCounter] = alpha[constantTerm3] + beta[constantTerm3]*bk[constantTerm3]*bk[constantTerm3];
		     rowValueCounter = rowValueCounter + 1;
		     constantTerm3 = constantTerm3 + 1;
	     }
	     if(distance<count3){
		    constantOfMatrix[constantOfMatrixCounter] = 0;
		    constantOfMatrixCounter = constantOfMatrixCounter + 1;
	     }
     }

     for (int forceIn=count3; forceIn<count; forceIn++){
		  constantOfMatrix[constantOfMatrixCounter] = 1;
		  constantOfMatrixCounter = constantOfMatrixCounter + 1;
		  constantOfMatrix[constantOfMatrixCounter] = 0;
		  constantOfMatrixCounter = constantOfMatrixCounter + 1;
     }
}


//function for solving linear set of equation using GSL package
void linearEquationSolver()
{
     gsl_matrix_view m 
    = gsl_matrix_view_array (rowValue, constantOfMatrixCounter, constantOfMatrixCounter);
    
  gsl_vector_view b
    = gsl_vector_view_array (constantOfMatrix, constantOfMatrixCounter);

  gsl_vector *fx = gsl_vector_alloc (constantOfMatrixCounter);
  
  int s;

  gsl_permutation * p = gsl_permutation_alloc (constantOfMatrixCounter);

  gsl_linalg_LU_decomp (&m.matrix, p, &s);

  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, fx);

  unsigned int constantOfMatrixCounter1;
  constantOfMatrixCounter1 = constantOfMatrixCounter;
  vector<double> force;
  for (size_t t = 0; t != constantOfMatrixCounter1; t++) {
	  force.push_back(gsl_vector_get (fx, t)); 
  }

  gsl_permutation_free (p);
  gsl_vector_free (fx); 
 
  int forceCounter = 0;
  for(vector<double>::const_iterator t = force.begin(); t !=force.end(); t++) {
	  forceValue[forceCounter] = *t;
	  forceCounter = forceCounter + 1;
		
  }	
}


//function for main simulation
void mainSimulation()
{
	 ofstream myfile;
     myfile.open ("particle_capture_cell_identity.txt ");
     ofstream out3("particle_capture_cell_identity.txt ", ios::app);
     myfile.open ("particle1_travelpath.xlsx");
     ofstream out33("particle1_travelpath.xlsx", ios::app);
     myfile.open ("particle2_travelpath.xlsx");
     ofstream out4("particle2_travelpath.xlsx", ios::app);
     myfile.open ("particle3_travelpath.xlsx");
     ofstream out5("particle3_travelpath.xlsx", ios::app);
     myfile.open ("particle4_travelpath.xlsx");
     ofstream out6("particle4_travelpath.xlsx", ios::app);
     myfile.open ("particle5_travelpath.xlsx");
     ofstream out7("particle5_travelpath.xlsx", ios::app);
     srand (time(NULL));
  

     for(int i=0; i<5; i++){
         double srk=0, srk1=0, radiOfNanoParticle1=0, radiOfNanoParticle2=0, delt=0, randomAngle=0, dmin=0, centerToPointDistance=0, ux=0, uy=0, nextStepDistance=0;
         double radius12=0,distance11=0,l11=0,l12=0,li=0,T=0,theta11=0,theta12=0,theta13=0,arcx=0,arcy=0;
         double dummy=0;
         t=0,ypx1=0, ypy1=0, ypx=0, ypy=0, ypx2=0, ypy2=0;  
         
         radiOfNanoParticle1 = radiOfNanoParticle;
         radiOfNanoParticle2 = radiOfNanoParticle1*1e-3;
         
         double D = (((1.3807e-23)*(310))/(6*PI*(2.5e-3)*(radiOfNanoParticle1*1e-9)))*(1e12);
         ypy1 = (sizeY - 0) * ( (double)rand() / (double)RAND_MAX ) + 0;

         if(i==0){
            out33<<t<<" "<<ypx1<<" "<<ypy1<<endl;
         }
         if(i==1){
            out4<<t<<" "<<ypx1<<" "<<ypy1<<endl;
         }
         if(i==2){
            out5<<t<<" "<<ypx1<<" "<<ypy1<<endl;
         }
         if(i==3){
            out6<<t<<" "<<ypx1<<" "<<ypy1<<endl;
         }
         if(i==4){
            out7<<t<<" "<<ypx1<<" "<<ypy1<<endl;
         }
       
         while (t<10000){
	            delt=bulkDelt;
                ux=0, uy=0;
                int sConstantTerm = 0, sForceCounter=0, sForceCounter1=1;
                randomAngle = (360 - 0) * ( (double)rand() / (double)RAND_MAX ) + 0;

                //code for velocity calculation start
                for(int sDistance=0; sDistance<count; sDistance++){
	                sak[sConstantTerm] = ypx1-aPointx[sDistance];
	                sbk[sConstantTerm] = ypy1-aPointy[sDistance];
	                srk1 = sak[sConstantTerm]*sak[sConstantTerm]+sbk[sConstantTerm]*sbk[sConstantTerm];
	                srk = sqrt(srk1);
	                salpha[sConstantTerm] = (1/(4*PI*(2.5)))*(-(0.5*(log((srk*srk)+(0.25))))+((0.25)/((srk*srk)+0.25)));
	                sbeta[sConstantTerm] = (1/(4*PI*(2.5)))*(1/((srk*srk)+(0.25)));
	                sConstantTerm = sConstantTerm + 1;
                }
                for(int sDistance1=0; sDistance1<count; sDistance1++){
	                ux = ux + (salpha[sDistance1] + sbeta[sDistance1]*sak[sDistance1]*sak[sDistance1])*forceValue[sForceCounter] + sbeta[sDistance1]*sak[sDistance1]*sbk[sDistance1]*forceValue[sForceCounter1];
	                uy = uy + sbeta[sDistance1]*sak[sDistance1]*sbk[sDistance1]*forceValue[sForceCounter] + (salpha[sDistance1] + sbeta[sDistance1]*sbk[sDistance1]*sbk[sDistance1])*forceValue[sForceCounter1];
	                sForceCounter = sForceCounter + 2;
	                sForceCounter1 = sForceCounter1 + 2;
                }  //code for velocity calculation end
  
                // code for calculation of minimum diameter start
  	            for(int centerDistanceCounter3 = 0; centerDistanceCounter3<n; centerDistanceCounter3++){
		            cbeta[centerDistanceCounter3] = ((ypx1 - x[centerDistanceCounter3])*(ypx1 - x[centerDistanceCounter3])  + (ypy1 - y[centerDistanceCounter3])*(ypy1 - y[centerDistanceCounter3]));
		            cbeta[centerDistanceCounter3] = sqrt (cbeta[centerDistanceCounter3])- r - radiOfNanoParticle2;
		        }
		        dmin = cbeta[0];
		        for (int centerDistanceCounter3 = 0; centerDistanceCounter3<n; centerDistanceCounter3++){
			         if(dmin>cbeta[centerDistanceCounter3]){
			            dmin = cbeta[centerDistanceCounter3];
		             }
		        }
		        double v = sqrt((ux*ux) + (uy*uy));


  
               // code for advancement of nano drug particle start		 
               ypx=ypx1;
               ypy=ypy1;
               ypx1 = ypx1 + (ux*delt) + (sqrt(4*D*delt))*cos( randomAngle* PI / 180.0 );
               ypy1 = ypy1 + (uy*delt) + (sqrt(4*D*delt))*sin( randomAngle* PI / 180.0 );
  
  
               //code for out of tissue boundary start
               while(ypx1<0){
	                 ypx1=ypx;
	                 ypy1=ypy;
	                 randomAngle = (360 - 0) * ( (double)rand() / (double)RAND_MAX ) + 0;
	                 ypx1 = ypx1 + (ux*delt) + (sqrt(4*D*delt))*cos( randomAngle* PI / 180.0 );
                     ypy1 = ypy1 + (uy*delt) + (sqrt(4*D*delt))*sin( randomAngle* PI / 180.0 );
                     if(ypx1>0){
		                goto skip01;
	                 }	   
                 }
                 skip01:
                if(ypy1>sizeY || ypy1<0){
	               ypx2=ypx1;
	               ypy2=ypy1;
	               ypx1 = ypx1;
	               ypy1 = ypy;
 	               for (int centerDistanceCounter = 0; centerDistanceCounter<n; centerDistanceCounter++){
		                centerToPointDistance = sqrt((ypx1 - x[centerDistanceCounter])*(ypx1 - x[centerDistanceCounter])  + (ypy1 - y[centerDistanceCounter])*(ypy1 - y[centerDistanceCounter]));
                        if(centerToPointDistance<(r+radiOfNanoParticle2)){
			               //reflective boundary condition 1
			               radius12=0,distance11=0,l11=0,l12=0,li=0,T=0,theta11=0,theta12=0,theta13=0,arcx=0,arcy=0;
			               radius12 = sqrt(((x[centerDistanceCounter]-ypx2)*(x[centerDistanceCounter]-ypx2))+((y[centerDistanceCounter]-ypy2)*(y[centerDistanceCounter]-ypy2)));
			               distance11 = radius12-r;
			               theta11 = abs(atan((ypy1-ypy2)/(ypx1-ypx2)) - atan((y[centerDistanceCounter]-ypy2)/(x[centerDistanceCounter]-ypx2)));
			               l11 = radius12*cos ( theta11* PI / 180.0 ) - distance11;
			               arcx = (distance11*ypx1 + l11*ypx2)/(distance11 + l11);
			               arcy = (distance11*ypy1 + l11*ypy2)/(distance11 + l11);
			               theta12 = atan((y[centerDistanceCounter]-arcy)/(x[centerDistanceCounter]-arcx));
			               li = sqrt(((ypx1-arcx)*(ypx1-arcx)) + ((ypy1-arcy)*(ypy1-arcy)));
			               theta13 = abs(atan((ypy1-arcy)/(ypx1-arcx)) - atan((y[centerDistanceCounter]-arcy)/(x[centerDistanceCounter]-arcx)));
			               l12 = li*cos ( theta13* PI / 180.0 );
			               T = 2*l12;
			               if(ypx1==arcx && ypx1==x[centerDistanceCounter] && ypy2< arcy){
				    	      ypx1 = ypx1;
					          ypy1 = ypy1 - T;
					          goto skip1;
				           }
				    	   if(ypx1==arcx &&  ypx1==x[centerDistanceCounter] && ypy2> arcy){
					          ypx1 = ypx1;
					          ypy1 = ypy1 + T;
					          goto skip1;
					       }
					       if(ypx2<arcx && ypy1==arcy && ypy1==y[centerDistanceCounter]){
					          ypx1 = ypx1 - T;
					          ypy1 = ypy1;
					          goto skip1;
					       }
					       if(ypx2>arcx && ypy1==arcy && ypy1==y[centerDistanceCounter]){
					          ypx1 = ypx1 + T;
					          ypy1 = ypy1;
					          goto skip1;
					       }
					       if(arcx>x[centerDistanceCounter]){
					          ypx1 = ypx1 + T*cos ( theta12* PI / 180.0 );
					          ypy1 = ypy1 + T*sin ( theta12* PI / 180.0 );
					       }
					       if(arcx<x[centerDistanceCounter]){
					          ypx1 = ypx1 - T*cos ( theta12* PI / 180.0 );
					          ypy1 = ypy1 - T*sin ( theta12* PI / 180.0 );
					       }
					       skip1:
					       dummy=dummy+0;  
                       }
	              }
              }	   // code for out of tissue boundary end
  
              nextStepDistance = sqrt((ypx1 - ypx)*(ypx1 - ypx)  + (ypy1 - ypy)*(ypy1 - ypy));
  
  
              //code for near solid boundary
              if(dmin<nextStepDistance){
	             delt = (((dmin+(D/v)) - sqrt(((dmin+(D/v))*(dmin+(D/v)))-(dmin)*(dmin)))/v);
		         if(dmin<=4*radiOfNanoParticle2){
		            delt = 0.001;
	             }
                 randomAngle = (360 - 0) * ( (double)rand() / (double)RAND_MAX ) + 0;
                 ypx1=ypx;
                 ypy1=ypy;
                 ypx1 = ypx1 + (ux*delt) + (sqrt(4*D*delt))*cos( randomAngle* PI / 180.0 );
                 ypy1 = ypy1 + (uy*delt) + (sqrt(4*D*delt))*sin( randomAngle* PI / 180.0 );
 
              }



              if(i==0){
                 out33<<t<<" "<<ypx1<<" "<<ypy1<<endl;
              }
              if(i==1){
                 out4<<t<<" "<<ypx1<<" "<<ypy1<<endl;
              }
              if(i==2){
              out5<<t<<" "<<ypx1<<" "<<ypy1<<endl;
              }
              if(i==3){
                 out6<<t<<" "<<ypx1<<" "<<ypy1<<endl;
              }
              if(i==4){
                 out7<<t<<" "<<ypx1<<" "<<ypy1<<endl;
              }
              if(ypx1>sizeX){
				 goto finish;  
			  }
			  
              t = t + delt;
 
  
             // code for cellular uptake or reflection start
             for (int centerDistanceCounter = 0; centerDistanceCounter<n; centerDistanceCounter++){
		          centerToPointDistance = sqrt((ypx1 - x[centerDistanceCounter])*(ypx1 - x[centerDistanceCounter])  + (ypy1 - y[centerDistanceCounter])*(ypy1 - y[centerDistanceCounter]));
		          if(centerToPointDistance<=(r+radiOfNanoParticle2)){
			         double randompick = 0;
			         randompick = (1 - 0) * ( (double)rand() / (double)RAND_MAX ) + 0;
			         if(randompick <= probabilityOfGettingCapture ){
			            out3<< "drug molecule is captured"<<endl<<"center of the cell which captured the drug molecule is"<<endl;
			            out3<<x[centerDistanceCounter]<<endl<<y[centerDistanceCounter]<<endl;
			            t = t + delt;
			            out3<<"time required = "<<t<<" sec"<<endl;
			            goto finish;
		             }
                     // reflective boundary conditions
		             if (randompick > probabilityOfGettingCapture ){
			             radius12=0,distance11=0,l11=0,l12=0,li=0,T=0,theta11=0,theta12=0,theta13=0,arcx=0,arcy=0;
			             radius12 = sqrt(((x[centerDistanceCounter]-ypx)*(x[centerDistanceCounter]-ypx))+((y[centerDistanceCounter]-ypy)*(y[centerDistanceCounter]-ypy)));
			             distance11 = radius12-r;
			             theta11 = abs(atan((ypy1-ypy)/(ypx1-ypx)) - atan((y[centerDistanceCounter]-ypy)/(x[centerDistanceCounter]-ypx)));
			             l11 = radius12*cos ( theta11* PI / 180.0 ) - distance11;
			             arcx = (distance11*ypx1 + l11*ypx)/(distance11 + l11);
			             arcy = (distance11*ypy1 + l11*ypy)/(distance11 + l11);
			             theta12 = atan((y[centerDistanceCounter]-arcy)/(x[centerDistanceCounter]-arcx));
			             li = sqrt(((ypx1-arcx)*(ypx1-arcx)) + ((ypy1-arcy)*(ypy1-arcy)));
			             theta13 = abs(atan((ypy1-arcy)/(ypx1-arcx)) - atan((y[centerDistanceCounter]-arcy)/(x[centerDistanceCounter]-arcx)));
			             l12 = li*cos ( theta13* PI / 180.0 );
			             T = 2*l12;
			             if(ypx1==arcx && ypx1==x[centerDistanceCounter] && ypy< arcy){
					     	ypx1 = ypx1;
						    ypy1 = ypy1 - T;
						    goto skip11;
					     }
					     if(ypx1==arcx && ypx1==x[centerDistanceCounter] && ypy> arcy){
						    ypx1 = ypx1;
						    ypy1 = ypy1 + T;
						    goto skip11;
					     }
					     if(ypx<arcx && ypy1==arcy && ypy1==y[centerDistanceCounter]){
						    ypx1 = ypx1 - T;
						    ypy1 = ypy1;
						    goto skip11;
					     }
					     if(ypx>arcx && ypy1==arcy && ypy1==y[centerDistanceCounter] ){
						    ypx1 = ypx1 + T;
						    ypy1 = ypy1;
						    goto skip11;
					     }
					     if(arcx>x[centerDistanceCounter]){
						    ypx1 = ypx1 + T*cos ( theta12* PI / 180.0 );
						    ypy1 = ypy1 + T*sin ( theta12* PI / 180.0 );
					     }
					     if(arcx<x[centerDistanceCounter]){
						    ypx1 = ypx1 - T*cos ( theta12* PI / 180.0 );
						    ypy1 = ypy1 - T*sin ( theta12* PI / 180.0 );
					     }
					     skip11:
					     dummy=dummy; 
		             }
		                   
                     if(i==0){
                        out33<<t<<" "<<ypx1<<" "<<ypy1<<endl;
                     }
                     if(i==1){
                        out4<<t<<" "<<ypx1<<" "<<ypy1<<endl;
                     }
                     if(i==2){
                        out5<<t<<" "<<ypx1<<" "<<ypy1<<endl;
                     }
                     if(i==3){
                        out6<<t<<" "<<ypx1<<" "<<ypy1<<endl;
                     }
                     if(i==4){
                        out7<<t<<" "<<ypx1<<" "<<ypy1<<endl;
                     }
	                  		                   
                     t = t+delt;
                     if(ypx1>sizeX){
				        goto finish;  
			         }
                           
	              } //if loop

                         
                  if(t>10000){
	                 goto finish2;
                  }
             } //for loop  // code for cellular uptake or reflection end

             finish2: 
             dummy=dummy;
             if(t>timechecker){
				 
                cout<<"Particle No."<<i+1<<endl<<"time of simulation passed = "<<t<<endl;
                timechecker = timechecker + 10;
			}  
         } //while loop  
         
     finish:
     dummy=dummy;
     timechecker=10; 

     } //inner loop

// code for advancement of nano drug particle end

	
  myfile.close();	
}


  
