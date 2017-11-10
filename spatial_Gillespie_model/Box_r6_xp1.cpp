// Spatial gillespie as used in Zeng and Holmes J Neurophysiol 103:1798-1808, 2010
// modified for a geometry of a box instead of a spine for use in model 
// verification by Li and Holmes 2017.  Here CaMKII has 6 subunits in a ring.
// Original code by Shangyou Zeng, modified by William Holmes

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "zeng_canngr6.h"

using namespace std;

double const D_CaMKII_C1 = 0.002;
double const D_CaMKII_C2 = 0;
double const D_Ca_C = 0.223; 
double const D_CaM00_C = 0.223/4;
double const D_CaM01_C = 0.223/4;
double const D_CaM10_C = 0.223/4;
double const D_CaM11_C = 0.223/4;
double const D_CaM02_C = 0.223/4;
double const D_CaM20_C = 0.223/4;
double const D_CaM12_C = 0.223/4;
double const D_CaM21_C = 0.223/4;
double const D_CaM22_C = 0.223/4;
//
struct voxel box[10][10][10];
void initial();
struct CaMKII *creat(int n);
double transition_lamda();
double box_diffusion_lamda(int i, int j, int k);
double voxel_chemical_lamda(struct voxel grid);
void reaction_box(double choose2, int p, int q, int r);
void chemical_reaction(struct voxel *p1, double choose2, int index);
int CaMKII_number(struct voxel grid);
double ran2(long *idum);
void reaction();
struct CaMKII * delete_CaMKII(struct voxel grid, struct CaMKII *point);
struct CaMKII * insert_CaMKII(struct voxel grid, struct CaMKII *point);

double reaction_lamda;
double random_num;
long *idum;
int x;
int y;
int z;
//

int main(){
  int i, j, k, m, n, num, index, flag;
  struct CaMKII *p;
  long l;
  double time, Ca_lamda, Ca_lamda1, total_Ca, total_Ca1;
  FILE *fp1;
  fp1=fopen("Box_r6-xp1.dat","wb"); 

  int Ca, Ca1, Ca2, Ca4, CaMKII, CaMKIICaM00, CaMKIICaM10, CaMKIICaM01, CaMKIICaM11, CaMKIICaM20, CaMKIICaM02, CaMKIICaM12, CaMKIICaM21, CaMKIICaM22;
  int Trapped00, Trapped01, Trapped10, Trapped11, Trapped02, Trapped20, Trapped12, Trapped21, Trapped22, Auton, Capped;
  int CaM00, CaM00_1, CaM01, CaM10, CaM11, CaM20, CaM02, CaM12, CaM21, CaM22;
  int Trapped0, Trapped1, Trapped2, Trapped3;
  int bound, trapped;
  
  l = 1000;
  idum = &l;
  time = 0;
  num = 0;
  Ca = 0;
  Ca1 = 0;
  Ca2 = 0;
  Ca4 = 0;
  CaMKII = 0;
  CaMKIICaM00 = 0;
  CaMKIICaM01 = 0;
  CaMKIICaM10 = 0;
  CaMKIICaM11 = 0;
  CaMKIICaM20 = 0;
  CaMKIICaM02 = 0;
  CaMKIICaM12 = 0;
  CaMKIICaM21 = 0;
  CaMKIICaM22 = 0;
  Trapped00 = 0;
  Trapped01 = 0;
  Trapped10 = 0;
  Trapped11 = 0;
  Trapped02 = 0;
  Trapped20 = 0;
  Trapped12 = 0;
  Trapped21 = 0;
  Trapped22 = 0;
  Auton = 0;
  Capped = 0;
  CaM00 = 0;
  CaM00_1 = 0;
  CaM01 = 0;
  CaM10 = 0;
  CaM11 = 0;
  CaM02 = 0;
  CaM20 = 0;
  CaM12 = 0;
  CaM21 = 0;
  CaM22 = 0;
  Trapped0 = 0;
  Trapped1 = 0;
  Trapped2 = 0;
  Trapped3 = 0;
  bound = 0;
  trapped = 0;
  flag = 1;
  
  initial();

  reaction_lamda = transition_lamda();
  Ca_lamda1 = 0;
  total_Ca = 0;
  while(time < 30000){
    total_Ca1 = total_Ca;

//  may have to change this later

    time += log(1.0/ran2(idum)) / reaction_lamda;
//    cout << "time = " << time << " " << reaction_lamda <<endl;
    if(time > (num * 1)){
      for(i=0; i<10; i++)
	for(j=0; j<10; j++)
           for (k=0; k<10; k++){
  	     Ca += box[i][j][k].Ca;
	     Ca1 += box[i][j][k].Ca;
	     Ca2 += box[i][j][k].Ca;
	     CaM00 += box[i][j][k].CaM00;
	     CaM00_1 += box[i][j][k].CaM00;
	     CaM01 += box[i][j][k].CaM01;
             CaM10 += box[i][j][k].CaM10;
	     CaM11 += box[i][j][k].CaM11;
	     CaM02 += box[i][j][k].CaM02;
	     CaM20 += box[i][j][k].CaM20;
	     CaM12 += box[i][j][k].CaM12;
	     CaM21 += box[i][j][k].CaM21;
	     CaM22 += box[i][j][k].CaM22;
	  
	     p = box[i][j][k].head;
	     while(p != NULL){
	       for(m=0; m<1; m++)
	         for(n=0; n<6; n++){
		   if(p->subunit[m][n] == 0){
		     CaMKII++;
		   }else if(p->subunit[m][n] == 1){
		     CaMKIICaM00++;
                     bound++;
		   }else if(p->subunit[m][n] == 2){
		     CaMKIICaM10++;
                     bound++;
		   }else if(p->subunit[m][n] == 3){
		     CaMKIICaM01++;
                     bound++;
		   }else if(p->subunit[m][n] == 4){
		     CaMKIICaM11++;
                     bound++;
		   }else if(p->subunit[m][n] == 5){
		     CaMKIICaM20++;
                     bound++;
		   }else if(p->subunit[m][n] == 6){
		     CaMKIICaM02++;
                     bound++;
		   }else if(p->subunit[m][n] == 7){
		     CaMKIICaM21++;
                     bound++;
		   }else if(p->subunit[m][n] == 8){
		     CaMKIICaM12++;
                     bound++;
		   }else if(p->subunit[m][n] == 9){
		     CaMKIICaM22++;
                     bound++;
		   }else if(p->subunit[m][n] == 10){
		     Trapped00++;
                     trapped++;
		   }else if(p->subunit[m][n] == 11){
		     Trapped10++;
                     trapped++;
		   }else if(p->subunit[m][n] == 12){
		     Trapped01++;
                     trapped++;
		   }else if(p->subunit[m][n] == 13){
		     Trapped11++;
                     trapped++;
		   }else if(p->subunit[m][n] == 14){
		     Trapped20++;
                     trapped++;
		   }else if(p->subunit[m][n] == 15){
		     Trapped02++;
                     trapped++;
		   }else if(p->subunit[m][n] == 16){
		     Trapped21++;
                     trapped++;
		   }else if(p->subunit[m][n] == 17){
		     Trapped12++;
                     trapped++;
		   }else if(p->subunit[m][n] == 18){
		     Trapped22++;
                     trapped++;
		   }else if(p->subunit[m][n] == 19){
		     Auton++;
		   }else if(p->subunit[m][n] == 20){
		     Capped++;
	 	   }
	        }  // n
	     p = p->next;
	  }  // while(p
	}  // k
        
      
       fprintf(fp1, "%f %d %f %d %d %d  %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d  %d %d %d %d  %d %d %d %d %d %d %d %d %d %d\n", time, Ca2, Ca2/7.5, Ca, Ca1, CaMKII, CaMKIICaM00, CaMKIICaM10, CaMKIICaM01, CaMKIICaM11, CaMKIICaM20, CaMKIICaM02, CaMKIICaM21, CaMKIICaM12, CaMKIICaM22, Trapped00, Trapped10, Trapped01, Trapped11, Trapped20, Trapped02, Trapped21, Trapped12, Trapped22, bound, trapped, Auton, Capped, CaM00, CaM10, CaM01, CaM11, CaM20, CaM02, CaM21, CaM12, CaM22, CaM00_1);

       num++;

       Ca = 0;
       Ca1 = 0;
       Ca2 = 0;
       CaMKII = 0;
       CaMKIICaM00 = 0;
       CaMKIICaM01 = 0;
       CaMKIICaM10 = 0;
       CaMKIICaM11 = 0;
       CaMKIICaM20 = 0;
       CaMKIICaM02 = 0;
       CaMKIICaM12 = 0;
       CaMKIICaM21 = 0;
       CaMKIICaM22 = 0;
       Trapped00 = 0;
       Trapped01 = 0;
       Trapped10 = 0;
       Trapped11 = 0;
       Trapped02 = 0;
       Trapped20 = 0;
       Trapped12 = 0;
       Trapped21 = 0;
       Trapped22 = 0;
       Auton = 0;
       Capped = 0;
       CaM00 = 0;
       CaM00_1 = 0;
       CaM01 = 0;
       CaM10 = 0;
       CaM11 = 0;
       CaM02 = 0;
       CaM20 = 0;
       CaM12 = 0;
       CaM21 = 0;
       CaM22 = 0;
//
       Trapped0 = 0;
       Trapped1 = 0;
       Trapped2 = 0;
       Trapped3 = 0;
	bound = 0;
	trapped = 0;

       cout << "time = " << time << " " << reaction_lamda << endl;
    }
    reaction();
  }
    
  return 0;
}

void initial(){
  int i, j, k, m, inter;
  double random_number;
  struct CaMKII *head;
  /******************************************************************************/

//   initialize  diffusion in box

  for(i=0; i<10; i++)
    for(j=0; j<10; j++)
      for(k=0; k<10; k++){
         for(m=0; m<6; m++){
//	box[i][j][k].D_CaMKII[m] = D_CaMKII_C1 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	box[i][j][k].D_CaMKII[m] = D_CaMKII_C2 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	box[i][j][k].D_Ca[m] = D_Ca_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	box[i][j][k].D_CaM00[m] = D_CaM00_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
        box[i][j][k].D_CaM10[m] = D_CaM10_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	box[i][j][k].D_CaM01[m] = D_CaM01_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	box[i][j][k].D_CaM11[m] = D_CaM11_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	box[i][j][k].D_CaM02[m] = D_CaM02_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	box[i][j][k].D_CaM20[m] = D_CaM20_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	box[i][j][k].D_CaM21[m] = D_CaM21_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	box[i][j][k].D_CaM12[m] = D_CaM12_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	box[i][j][k].D_CaM22[m] = D_CaM22_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
        }
      }
    

// Set boundaries  top top=0  m=4
  for(i=0; i<10; i++)
    for(j=0; j<10; j++){
      box[i][j][0].D_CaMKII[4] = 0;
      box[i][j][0].D_Ca[4] =  0;
      box[i][j][0].D_CaM00[4] = 0;
      box[i][j][0].D_CaM10[4] = 0;
      box[i][j][0].D_CaM01[4] = 0;
      box[i][j][0].D_CaM11[4] = 0;
      box[i][j][0].D_CaM02[4] = 0;
      box[i][j][0].D_CaM20[4] = 0;
      box[i][j][0].D_CaM21[4] = 0;
      box[i][j][0].D_CaM12[4] = 0;
      box[i][j][0].D_CaM22[4] = 0;
    }
     

  // bottom bottom = 9  m=5
 
  for(i=0; i<10; i++)
    for(j=0; j<10; j++){
      box[i][j][9].D_CaMKII[5] = 0;
      box[i][j][9].D_Ca[5] = 0;
      box[i][j][9].D_CaM00[5] = 0;
      box[i][j][9].D_CaM10[5] = 0;
      box[i][j][9].D_CaM01[5] = 0;
      box[i][j][9].D_CaM11[5] = 0;
      box[i][j][9].D_CaM02[5] = 0;
      box[i][j][9].D_CaM20[5] = 0;
      box[i][j][9].D_CaM21[5] = 0;
      box[i][j][9].D_CaM12[5] = 0;
      box[i][j][9].D_CaM22[5] = 0;
    }
    

//set the boudary  i=0 side m=0
  for(j=0; j<10; j++)
     for(k=0; k<10; k++){
        box[0][j][k].D_CaMKII[0] = 0;
        box[0][j][k].D_Ca[0] = 0;
        box[0][j][k].D_CaM00[0] = 0;
	box[0][j][k].D_CaM10[0] = 0;
	box[0][j][k].D_CaM01[0] = 0;
	box[0][j][k].D_CaM11[0] = 0;
	box[0][j][k].D_CaM02[0] = 0;
	box[0][j][k].D_CaM20[0] = 0;
	box[0][j][k].D_CaM21[0] = 0;
	box[0][j][k].D_CaM12[0] = 0;
	box[0][j][k].D_CaM22[0] = 0;
    }
//  J=0 side  m=1
  
  for(i=0; i<10; i++)
     for(k=0; k<10; k++){
        box[i][0][k].D_CaMKII[1] = 0;
        box[i][0][k].D_Ca[1] = 0;
        box[i][0][k].D_CaM00[1] = 0;
	box[i][0][k].D_CaM01[1] = 0;
	box[i][0][k].D_CaM10[1] = 0;
	box[i][0][k].D_CaM11[1] = 0;
	box[i][0][k].D_CaM20[1] = 0;
	box[i][0][k].D_CaM02[1] = 0;
	box[i][0][k].D_CaM21[1] = 0;
	box[i][0][k].D_CaM12[1] = 0;
	box[i][0][k].D_CaM22[1] = 0;    
  }

// i=9 side  m=2
  for(j=0; j<10; j++)
     for(k=0; k<10; k++){
        box[9][j][k].D_CaMKII[2] = 0;
        box[9][j][k].D_Ca[2] = 0;
        box[9][j][k].D_CaM00[2] = 0;
	box[9][j][k].D_CaM10[2] = 0;
	box[9][j][k].D_CaM01[2] = 0;
	box[9][j][k].D_CaM11[2] = 0;
	box[9][j][k].D_CaM20[2] = 0;
	box[9][j][k].D_CaM02[2] = 0;
	box[9][j][k].D_CaM12[2] = 0;
	box[9][j][k].D_CaM21[2] = 0;
	box[9][j][k].D_CaM22[2] = 0;    
  }
//    J=9 side m=3 
  for(i=0; i<10; i++)
     for(k=0; k<10; k++){
        box[i][9][k].D_CaMKII[3] = 0;
        box[i][9][k].D_Ca[3] = 0;
        box[i][9][k].D_CaM00[3] = 0;
	box[i][9][k].D_CaM01[3] = 0;
	box[i][9][k].D_CaM10[3] = 0;
	box[i][9][k].D_CaM11[3] = 0;
	box[i][9][k].D_CaM20[3] = 0;
	box[i][9][k].D_CaM02[3] = 0;
	box[i][9][k].D_CaM21[3] = 0;
	box[i][9][k].D_CaM12[3] = 0;
	box[i][9][k].D_CaM22[3] = 0;    
   }
  /***************************************************************************************/
  //initialize the concentration
  for(i=0; i<10; i++)
    for(j=0; j<10; j++)
      for(k=0; k<10; k++){
      box[i][j][k].CaMKII_num = 0;
      box[i][j][k].Ca_pump = 0;
      box[i][j][k].Ca = 6;
      box[i][j][k].CaM00 = 3;
	  box[i][j][k].CaM10 = 0;
	  box[i][j][k].CaM01 = 0;
	  box[i][j][k].CaM11 = 0;
	  box[i][j][k].CaM02 = 0;
	  box[i][j][k].CaM20 = 0;
	  box[i][j][k].CaM21 = 0;
	  box[i][j][k].CaM12 = 0;
	  box[i][j][k].CaM22 = 0;
    }

  
  /*************************************************************/
  //initialize transition rate
  for(i=0; i<10; i++)
    for(j=0; j<10; j++)
       for(k=0; k<10; k++){
	box[i][j][k].CaM00_to_CaM10 = 5.02857 / (6.02 * pow(10, 5) * pow(0.1, 3));
	box[i][j][k].CaM10_to_CaM00 = 52.8 / 1000.0;
	box[i][j][k].CaM00_to_CaM01 = 95 / (6.02 * pow(10, 5) * pow(0.1, 3));
	box[i][j][k].CaM01_to_CaM00 = 3750.0 / 1000.0;
	box[i][j][k].CaM10_to_CaM11 = 95 / (6.02 * pow(10, 5) * pow(0.1, 3));
	box[i][j][k].CaM11_to_CaM10 = 3750.0 / 1000.0;
	box[i][j][k].CaM10_to_CaM20 = 5.02857 / (6.02 * pow(10, 5) * pow(0.1, 3));
	box[i][j][k].CaM20_to_CaM10 = 7.84 / 1000.0;
	box[i][j][k].CaM01_to_CaM11 = 5.02857 / (6.02 * pow(10, 5) * pow(0.1, 3));
	box[i][j][k].CaM11_to_CaM01 = 52.8 / 1000.0;
	box[i][j][k].CaM01_to_CaM02 = 95 / (6.02 * pow(10, 5) * pow(0.1, 3));
	box[i][j][k].CaM02_to_CaM01 = 750 / 1000.0;
	box[i][j][k].CaM11_to_CaM21 = 5.02857 / (6.02 * pow(10, 5) * pow(0.1, 3));
	box[i][j][k].CaM21_to_CaM11 = 7.84 / 1000.0;
	box[i][j][k].CaM11_to_CaM12 = 95 / (6.02 * pow(10, 5) * pow(0.1, 3));
	box[i][j][k].CaM12_to_CaM11 = 750.0 / 1000.0;
	box[i][j][k].CaM20_to_CaM21 = 95 / (6.02 * pow(10, 5) * pow(0.1, 3));
	box[i][j][k].CaM21_to_CaM20 = 3750.0 / 1000.0;
	box[i][j][k].CaM02_to_CaM12 = 5.02857 / (6.02 * pow(10, 5) * pow(0.1, 3));
	box[i][j][k].CaM12_to_CaM02 = 52.8 / 1000.0;
	box[i][j][k].CaM21_to_CaM22 = 95 / (6.02 * pow(10, 5) * pow(0.1, 3));
	box[i][j][k].CaM22_to_CaM21 = 750.0 / 1000.0;
	box[i][j][k].CaM12_to_CaM22 = 5.02857 / (6.02 * pow(10, 5) * pow(0.1, 3));
	box[i][j][k].CaM22_to_CaM12 = 7.84 / 1000.0;

	box[i][j][k].CaM00_to_CaMKIICaM00 = 0.019 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM00_to_CaM00 = 5.586 / 1000.0;
        box[i][j][k].CaM10_to_CaMKIICaM10 = 0.295 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM10_to_CaM10 = 6.195 / 1000.0;
        box[i][j][k].CaM01_to_CaMKIICaM01 = 0.11 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM01_to_CaM01 = 3.234 / 1000.0;
        box[i][j][k].CaM11_to_CaMKIICaM11 = 1.65 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM11_to_CaM11 = 3.465 / 1000.0;
        box[i][j][k].CaM20_to_CaMKIICaM20 = 4.6 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM20_to_CaM20 = 6.9 / 1000.0;
        box[i][j][k].CaM02_to_CaMKIICaM02 = 0.6 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM02_to_CaM02 = 1.764 / 1000.0;
        box[i][j][k].CaM12_to_CaMKIICaM12 = 9.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM12_to_CaM12 = 1.995 / 1000.0;
        box[i][j][k].CaM21_to_CaMKIICaM21 = 26 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM21_to_CaM21 = 3.9 / 1000.0;
        box[i][j][k].CaM22_to_CaMKIICaM22 = 77.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM22_to_CaM22 = 1.1625 / 1000.0;

	box[i][j][k].CaMKIICaM00_to_CaMKIICaM10 = 44 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM10_to_CaMKIICaM00 = 33 / 1000.0; 
        box[i][j][k].CaMKIICaM00_to_CaMKIICaM01 = 76 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM01_to_CaMKIICaM00 = 300 / 1000.0;
        box[i][j][k].CaMKIICaM10_to_CaMKIICaM11 = 76 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM11_to_CaMKIICaM10 = 300 / 1000.0;
        box[i][j][k].CaMKIICaM10_to_CaMKIICaM20 = 44 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM20_to_CaMKIICaM10 = 4.9 / 1000.0;
        box[i][j][k].CaMKIICaM01_to_CaMKIICaM11 = 44 / (6.02 * pow(10, 5) * pow(0.1, 3));
	box[i][j][k].CaMKIICaM11_to_CaMKIICaM01 = 33 / 1000.0;
        box[i][j][k].CaMKIICaM01_to_CaMKIICaM02 = 76 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM02_to_CaMKIICaM01 = 60 / 1000.0;
        box[i][j][k].CaMKIICaM11_to_CaMKIICaM21 = 44 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM21_to_CaMKIICaM11 = 4.9 / 1000.0;
        box[i][j][k].CaMKIICaM11_to_CaMKIICaM12 = 76 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM12_to_CaMKIICaM11 = 60 / 1000.0;
        box[i][j][k].CaMKIICaM20_to_CaMKIICaM21 = 76 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM21_to_CaMKIICaM20 = 300 / 1000.0;
        box[i][j][k].CaMKIICaM02_to_CaMKIICaM12 = 44 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM12_to_CaMKIICaM02 = 33 / 1000.0;
        box[i][j][k].CaMKIICaM21_to_CaMKIICaM22 = 76 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM22_to_CaMKIICaM21 = 60 / 1000.0;
        box[i][j][k].CaMKIICaM12_to_CaMKIICaM22 = 44 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].CaMKIICaM22_to_CaMKIICaM12 = 4.9 / 1000.0;

//

        box[i][j][k].CaMKIICaM00_to_Trapped00 = 0 / 1000.0;
        box[i][j][k].Trapped00_to_CaMKIICaM00 = 0.003 / 1000.0;
        box[i][j][k].CaMKIICaM10_to_Trapped10 = 0.5 / 1000.0;
        box[i][j][k].Trapped10_to_CaMKIICaM10 = 0.003 / 1000.0;
        box[i][j][k].CaMKIICaM01_to_Trapped01 = 0.5 / 1000.0;
        box[i][j][k].Trapped01_to_CaMKIICaM01 = 0.003 / 1000.0;
        box[i][j][k].CaMKIICaM11_to_Trapped11 = 2.5 / 1000.0;
        box[i][j][k].Trapped11_to_CaMKIICaM11 = 0.003 / 1000.0;
        box[i][j][k].CaMKIICaM20_to_Trapped20 = 1/ 1000.0;
        box[i][j][k].Trapped20_to_CaMKIICaM20 = 0.003 / 1000.0;
        box[i][j][k].CaMKIICaM02_to_Trapped02 = 1 / 1000.0;
        box[i][j][k].Trapped02_to_CaMKIICaM02 = 0.003 / 1000.0;
        box[i][j][k].CaMKIICaM12_to_Trapped12 = 5 / 1000.0;
        box[i][j][k].Trapped12_to_CaMKIICaM12 = 0.003 / 1000.0;
        box[i][j][k].CaMKIICaM21_to_Trapped21 = 5 / 1000.0;
        box[i][j][k].Trapped21_to_CaMKIICaM21 = 0.003 / 1000.0;
        box[i][j][k].CaMKIICaM22_to_Trapped22 = 10 / 1000.0;
        box[i][j][k].Trapped22_to_CaMKIICaM22 = 0.003 / 1000.0;
  
        box[i][j][k].Trapped00_to_Auton = 0.167589 / 1000.0;
        box[i][j][k].Auton_to_Trapped00 = 0.019 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped10_to_Auton = 0.03717 / 1000.0;
        box[i][j][k].Auton_to_Trapped10 = 0.295 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped01_to_Auton = 0.019404 / 1000.0;
        box[i][j][k].Auton_to_Trapped01 = 0.11 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped11_to_Auton = 0.004158 / 1000.0;
        box[i][j][k].Auton_to_Trapped11 = 1.65 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped20_to_Auton = 0.0207 / 1000.0;
        box[i][j][k].Auton_to_Trapped20 = 4.6 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped02_to_Auton = 0.005292 / 1000.0;
        box[i][j][k].Auton_to_Trapped02 = 0.6 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped12_to_Auton = 0.001197 / 1000.0;
        box[i][j][k].Auton_to_Trapped12 = 9.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped21_to_Auton = 0.00234 / 1000.0;
        box[i][j][k].Auton_to_Trapped21 = 26 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped22_to_Auton = 0.000348755 / 1000.0;
        box[i][j][k].Auton_to_Trapped22 = 77.5 / (6.02 * pow(10, 5) * pow(0.1, 3));

        box[i][j][k].Trapped00_to_Trapped10 = 44 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped10_to_Trapped00 = 6.6 / 1000.0;
        box[i][j][k].Trapped00_to_Trapped01 = 76 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped01_to_Trapped00 = 60 / 1000.0;
        box[i][j][k].Trapped10_to_Trapped11 = 76 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped11_to_Trapped10 = 60 / 1000.0;
        box[i][j][k].Trapped10_to_Trapped20 = 44 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped20_to_Trapped10 = 2.45 / 1000.0;
        box[i][j][k].Trapped01_to_Trapped11 = 44 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped11_to_Trapped01 = 6.6 / 1000.0;
        box[i][j][k].Trapped01_to_Trapped02 = 76 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped02_to_Trapped01 = 30 / 1000.0;
        box[i][j][k].Trapped11_to_Trapped21 = 44 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped21_to_Trapped11 = 2.45 / 1000.0;
        box[i][j][k].Trapped11_to_Trapped12 = 76 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped12_to_Trapped11 = 30 / 1000.0;
        box[i][j][k].Trapped20_to_Trapped21 = 76 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped21_to_Trapped20 = 60 / 1000.0;
        box[i][j][k].Trapped02_to_Trapped12 = 44 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped12_to_Trapped02 = 6.6 / 1000.0;
        box[i][j][k].Trapped21_to_Trapped22 = 76 / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped22_to_Trapped21 = 30 / 1000.0;
        box[i][j][k].Trapped12_to_Trapped22 = 44  / (6.02 * pow(10, 5) * pow(0.1, 3));
        box[i][j][k].Trapped22_to_Trapped12 = 2.45 / 1000.0;

        box[i][j][k].Auton_to_CaMKII = 0.003 / 1000.0;
        box[i][j][k].Auton_to_Capped = 0.1 / 1000.0;
        box[i][j][k].Capped_to_Auton = 0.01 / 1000.0;   
		
    }
  /*****************************************************************************/
  //initialize CaMKII
 
  for(i=0; i<10; i++)
    for(j=0; j<10; j++)
       for(k=0; k<10; k++)
         box[i][j][k].head = NULL;

//  initialize CaMKII number here, 83 or 40
//  initialize CaMKII number here, one per voxel
    
    m=0;
    for(i=0; i<10; i++)
      for(j=0; j<10; j++)
        for(k=0; k<10; k++){
          box[i][j][k].head = creat(m);
          m++;
       }

/*
   for(i=0; i<83; i++){
      random_number = 100 * ran2(idum);
      inter = (int) random_number;
      j = inter / 10;
      k = inter % 10;
      random_number2 = 100 * ran2(idum);
      inter2 = (int) random_number2;
      m = inter2 / 10;
      if(box[m][j][k].head == NULL){
        box[m][j][k].head = creat(i);
      }else{
        head = creat(i);
        head->next = box[m][j][k].head;
        box[m][j][k].head = head;
      }
   }
*/

}

struct CaMKII * creat(int n){
  struct CaMKII * head;
  int j, k;
  
  head = (struct CaMKII *) malloc(sizeof(struct CaMKII));
  for(j=0; j<1; j++)
    for(k=0; k<6; k++)
      head->subunit[j][k] = 0;
  head->index = n;
  
  return(head);
}

double transition_lamda(){
  int i, j, k;
  double lamda;
  lamda = 0;

  for(i=0; i<10; i++)
    for(j=0; j<10; j++)
       for(k=0; k<10; k++){
         box[i][j][k].diffusion_lamda = box_diffusion_lamda(i, j, k);
         box[i][j][k].chemical_lamda = voxel_chemical_lamda(box[i][j][k]);
         lamda += box[i][j][k].diffusion_lamda + box[i][j][k].chemical_lamda;
    }
	
  return lamda;
}

double voxel_chemical_lamda(struct voxel grid){
  double lamda;
  struct CaMKII *p1;
  int i, j;
	
  lamda = 0;
  
  lamda += grid.CaM00 * grid.Ca * grid.CaM00_to_CaM10;
  lamda += grid.CaM10 * grid.CaM10_to_CaM00;
  lamda += grid.CaM00 * grid.Ca * grid.CaM00_to_CaM01;
  lamda += grid.CaM01 * grid.CaM01_to_CaM00;
  lamda += grid.CaM10 * grid.Ca * grid.CaM10_to_CaM11;
  lamda += grid.CaM11 * grid.CaM11_to_CaM10;
  lamda += grid.CaM10 * grid.Ca * grid.CaM10_to_CaM20;
  lamda += grid.CaM20 * grid.CaM20_to_CaM10;
  lamda += grid.CaM01 * grid.Ca * grid.CaM01_to_CaM11;
  lamda += grid.CaM11 * grid.CaM11_to_CaM01;
  lamda += grid.CaM01 * grid.Ca * grid.CaM01_to_CaM02;
  lamda += grid.CaM02 * grid.CaM02_to_CaM01;
  lamda += grid.CaM11 * grid.Ca * grid.CaM11_to_CaM21;
  lamda += grid.CaM21 * grid.CaM21_to_CaM11;
  lamda += grid.CaM11 * grid.Ca * grid.CaM11_to_CaM12;
  lamda += grid.CaM12 * grid.CaM12_to_CaM11;
  lamda += grid.CaM20 * grid.Ca * grid.CaM20_to_CaM21;
  lamda += grid.CaM21 * grid.CaM21_to_CaM20;
  lamda += grid.CaM02 * grid.Ca * grid.CaM02_to_CaM12;
  lamda += grid.CaM12 * grid.CaM12_to_CaM02;
  lamda += grid.CaM21 * grid.Ca * grid.CaM21_to_CaM22;
  lamda += grid.CaM22 * grid.CaM22_to_CaM21;
  lamda += grid.CaM12 * grid.Ca * grid.CaM12_to_CaM22;
  lamda += grid.CaM22 * grid.CaM22_to_CaM12;


  //chemical reaction related to CaMKII
  p1 = grid.head;
  
  while(p1 != NULL){
    for(i=0; i<1; i++)
      for(j=0; j<6; j++){
	if(p1->subunit[i][j] ==0){
	  lamda += grid.CaM00 * grid.CaM00_to_CaMKIICaM00;
	  lamda += grid.CaM10 * grid.CaM10_to_CaMKIICaM10;
	  lamda += grid.CaM01 * grid.CaM01_to_CaMKIICaM01;
	  lamda += grid.CaM11 * grid.CaM11_to_CaMKIICaM11;
	  lamda += grid.CaM20 * grid.CaM20_to_CaMKIICaM20;
	  lamda += grid.CaM02 * grid.CaM02_to_CaMKIICaM02;
	  lamda += grid.CaM12 * grid.CaM12_to_CaMKIICaM12;
	  lamda += grid.CaM21 * grid.CaM21_to_CaMKIICaM21;
	  lamda += grid.CaM22 * grid.CaM22_to_CaMKIICaM22;
	}else if(p1->subunit[i][j] == 1){
	  lamda += grid.CaMKIICaM00_to_CaM00;
	  lamda += grid.Ca * grid.CaMKIICaM00_to_CaMKIICaM01;
	  lamda += grid.Ca * grid.CaMKIICaM00_to_CaMKIICaM10;

	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
         lamda += grid.CaMKIICaM00_to_Trapped00;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20))
  || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
         (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
*/
//	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
//		 lamda += 0.4 * grid.CaMKIICaM00_to_Trapped00;
//
	}else if(p1->subunit[i][j] == 2){
	  lamda += grid.CaMKIICaM10_to_CaMKIICaM00;
	  lamda += grid.CaMKIICaM10_to_CaM10;
	  lamda += grid.Ca * grid.CaMKIICaM10_to_CaMKIICaM11;
          lamda += grid.Ca * grid.CaMKIICaM10_to_CaMKIICaM20;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
         lamda += grid.CaMKIICaM10_to_Trapped10;
/*
//	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
//		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20))
                                             || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
*/
//	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20))  // ||
//		 lamda += 0.4 * grid.CaMKIICaM10_to_Trapped10;
//
	}else if(p1->subunit[i][j] == 3){
	  lamda += grid.CaMKIICaM01_to_CaMKIICaM00;
	  lamda += grid.CaMKIICaM01_to_CaM01;
	  lamda += grid.Ca * grid.CaMKIICaM01_to_CaMKIICaM11;
	  lamda += grid.Ca * grid.CaMKIICaM01_to_CaMKIICaM02;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
         lamda += grid.CaMKIICaM01_to_Trapped01;
/*	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20))
                                              || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
         (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
*/
//  	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20)) 
//		 lamda += 0.4 * grid.CaMKIICaM01_to_Trapped01;
//
	}else if(p1->subunit[i][j] == 4){
	  lamda += grid.CaMKIICaM11_to_CaMKIICaM01;
          lamda += grid.CaMKIICaM11_to_CaMKIICaM10;
	  lamda += grid.CaMKIICaM11_to_CaM11;
	  lamda += grid.Ca * grid.CaMKIICaM11_to_CaMKIICaM21;
	  lamda += grid.Ca * grid.CaMKIICaM11_to_CaMKIICaM12;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
         lamda += grid.CaMKIICaM11_to_Trapped11;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] ==1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM11_to_Trapped11;
*/
//
//	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20))
// one way only     ||
//                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
//		 lamda += 0.4 * grid.CaMKIICaM11_to_Trapped11;
//
	}else if(p1->subunit[i][j] == 5){
	  lamda += grid.CaMKIICaM20_to_CaMKIICaM10;
	  lamda += grid.Ca * grid.CaMKIICaM20_to_CaMKIICaM21;
	  lamda += grid.CaMKIICaM20_to_CaM20;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
         lamda += grid.CaMKIICaM20_to_Trapped20;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM20_to_Trapped20;
*/
//
//	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20))
//   ||
//                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
//		 lamda += 0.4 * grid.CaMKIICaM20_to_Trapped20;
//
	}else if(p1->subunit[i][j] == 6){
	  lamda += grid.CaMKIICaM02_to_CaMKIICaM01;
	  lamda += grid.Ca * grid.CaMKIICaM02_to_CaMKIICaM12;
	  lamda += grid.CaMKIICaM02_to_CaM02;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
         lamda += grid.CaMKIICaM02_to_Trapped02;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM02_to_Trapped02;
*/
//
//	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20))
// ||
//                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
//		 lamda += 0.4 * grid.CaMKIICaM02_to_Trapped02;
//
	}else if(p1->subunit[i][j] == 7){
	  lamda += grid.CaMKIICaM21_to_CaMKIICaM11;
          lamda += grid.CaMKIICaM21_to_CaMKIICaM20;
	  lamda += grid.Ca * grid.CaMKIICaM21_to_CaMKIICaM22;
	  lamda += grid.CaMKIICaM21_to_CaM21;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
         lamda += grid.CaMKIICaM21_to_Trapped21;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM21_to_Trapped21;
*/
//
//	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20)) 
// ||
//                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
//		 lamda += 0.4 * grid.CaMKIICaM21_to_Trapped21;
//
	}else if(p1->subunit[i][j] == 8){
	  lamda += grid.CaMKIICaM12_to_CaMKIICaM11;
          lamda += grid.CaMKIICaM12_to_CaMKIICaM02;
	  lamda += grid.Ca * grid.CaMKIICaM12_to_CaMKIICaM22;
	  lamda += grid.CaMKIICaM12_to_CaM12;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
         lamda += grid.CaMKIICaM12_to_Trapped12;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM12_to_Trapped12;
*/
//
//	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20))
// ||
//                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
//		 lamda += 0.4 * grid.CaMKIICaM12_to_Trapped12;
//
	}else if(p1->subunit[i][j] == 9){
	  lamda += grid.CaMKIICaM22_to_CaMKIICaM21;
          lamda += grid.CaMKIICaM22_to_CaMKIICaM12;
	  lamda += grid.CaMKIICaM22_to_CaM22;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
         lamda += grid.CaMKIICaM22_to_Trapped22;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM22_to_Trapped22;
*/
//
//	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20))
//  ||
//                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
//		 lamda += 0.4 * grid.CaMKIICaM22_to_Trapped22;
//
	}else if(p1->subunit[i][j] == 10){
	  lamda += grid.Trapped00_to_CaMKIICaM00;
	  lamda += grid.Ca * grid.Trapped00_to_Trapped01;
	  lamda += grid.Ca * grid.Trapped00_to_Trapped10;
	  lamda += grid.Trapped00_to_Auton;
	 
	 }else if(p1->subunit[i][j] == 11){
	  lamda += grid.Trapped10_to_CaMKIICaM10;
	  lamda += grid.Trapped10_to_Trapped00;
          lamda += grid.Ca * grid.Trapped10_to_Trapped11;
	  lamda += grid.Ca * grid.Trapped10_to_Trapped20;
	  lamda += grid.Trapped10_to_Auton;

	}else if(p1->subunit[i][j] == 12){
	  lamda += grid.Trapped01_to_CaMKIICaM01;
	  lamda += grid.Trapped01_to_Trapped00;
          lamda += grid.Ca * grid.Trapped01_to_Trapped11;
	  lamda += grid.Ca * grid.Trapped01_to_Trapped02;
	  lamda += grid.Trapped01_to_Auton;

	}else if(p1->subunit[i][j] == 13){
	   lamda += grid.Trapped11_to_CaMKIICaM11;
	   lamda += grid.Trapped11_to_Trapped01;
	   lamda += grid.Trapped11_to_Trapped10;
	   lamda += grid.Ca * grid.Trapped11_to_Trapped21;
	   lamda += grid.Ca * grid.Trapped11_to_Trapped12;
	   lamda += grid.Trapped11_to_Auton;

	}else if(p1->subunit[i][j] == 14){
	   lamda += grid.Trapped20_to_CaMKIICaM20;
	   lamda += grid.Trapped20_to_Trapped10;
	   lamda += grid.Ca * grid.Trapped20_to_Trapped21;
	   lamda += grid.Trapped20_to_Auton;

	}else if(p1->subunit[i][j] == 15){
	   lamda += grid.Trapped02_to_CaMKIICaM02;
	   lamda += grid.Trapped02_to_Trapped01;
	   lamda += grid.Ca * grid.Trapped02_to_Trapped12;
	   lamda += grid.Trapped02_to_Auton;

	}else if(p1->subunit[i][j] == 16){
	   lamda += grid.Trapped21_to_CaMKIICaM21;
	   lamda += grid.Trapped21_to_Trapped11;
	   lamda += grid.Trapped21_to_Trapped20;
	   lamda += grid.Ca * grid.Trapped21_to_Trapped22;
	   lamda += grid.Trapped21_to_Auton;

	}else if(p1->subunit[i][j] == 17){
	   lamda += grid.Trapped12_to_CaMKIICaM12;
	   lamda += grid.Trapped12_to_Trapped11;
	   lamda += grid.Trapped12_to_Trapped02;
	   lamda += grid.Ca * grid.Trapped12_to_Trapped22;
	   lamda += grid.Trapped12_to_Auton;

	}else if(p1->subunit[i][j] == 18){
	   lamda += grid.Trapped22_to_CaMKIICaM22;
	   lamda += grid.Trapped22_to_Trapped21;
	   lamda += grid.Trapped22_to_Trapped12;
	   lamda += grid.Trapped22_to_Auton;

	}else if(p1->subunit[i][j] == 19){
	  lamda += grid.CaM00 * grid.Auton_to_Trapped00;
	  lamda += grid.CaM01 * grid.Auton_to_Trapped01;
	  lamda += grid.CaM10 * grid.Auton_to_Trapped10;
	  lamda += grid.CaM11 * grid.Auton_to_Trapped11;
	  lamda += grid.CaM20 * grid.Auton_to_Trapped20;
	  lamda += grid.CaM02 * grid.Auton_to_Trapped02;
	  lamda += grid.CaM21 * grid.Auton_to_Trapped21;
	  lamda += grid.CaM12 * grid.Auton_to_Trapped12;
	  lamda += grid.CaM22 * grid.Auton_to_Trapped22;
	  lamda += grid.Auton_to_CaMKII;

	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
	     lamda += grid.Auton_to_Capped;
/*
	   if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	     lamda += 0.4 * grid.Auton_to_Capped;
*/
//
//	   if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20))
//  ||
//                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
//	     lamda += 0.4 * grid.Auton_to_Capped;
//
	}else{
	  lamda += grid.Capped_to_Auton;
	}
	   
      }
    p1 = p1->next;
  }
  return lamda;
}

double box_diffusion_lamda(int p, int q, int r){
  double lamda;
  int k;
  int CaMKII_num;
  lamda = 0;
  
  CaMKII_num = CaMKII_number(box[p][q][r]);
  for(k=0; k<6; k++){
        lamda = lamda + CaMKII_num * box[p][q][r].D_CaMKII[k];
        lamda = lamda + box[p][q][r].Ca * box[p][q][r].D_Ca[k];

        lamda = lamda + box[p][q][r].CaM00 * box[p][q][r].D_CaM00[k];
        lamda = lamda + box[p][q][r].CaM10 * box[p][q][r].D_CaM10[k];
        lamda = lamda + box[p][q][r].CaM01 * box[p][q][r].D_CaM01[k];
        lamda = lamda + box[p][q][r].CaM11 * box[p][q][r].D_CaM11[k];
        lamda = lamda + box[p][q][r].CaM20 * box[p][q][r].D_CaM20[k];
	lamda = lamda + box[p][q][r].CaM02 * box[p][q][r].D_CaM02[k];
	lamda = lamda + box[p][q][r].CaM21 * box[p][q][r].D_CaM21[k];
	lamda = lamda + box[p][q][r].CaM12 * box[p][q][r].D_CaM12[k];
	lamda = lamda + box[p][q][r].CaM22 * box[p][q][r].D_CaM22[k];
  }
  return lamda;
}

int CaMKII_number(struct voxel grid){
  struct CaMKII *p1;
  int n;
  
  n = 0;
  p1 = grid.head;
  while(p1 != NULL){
    n++;
    p1 = p1->next;
  }

  return n;
}

void reaction(){
  int i, j, k;
  double choose, choose1;

  random_num = ran2(idum) * reaction_lamda;
  choose = 0;

  for(i=0; i<10; i++){
    for(j=0; j<10; j++){
       for(k=0; k<10; k++){
      choose1 = choose;
      choose = choose +  box[i][j][k].diffusion_lamda +  box[i][j][k].chemical_lamda;
      if((choose > random_num) && (choose1 <= random_num)){
	reaction_box(choose1, i, j, k);
	goto loop0;
        }
       }
     }
   }


 loop0: k = 0;
// cout<< " reacion   " << choose << " " << random_num <<" "<< reaction_lamda <<endl;    
}

void reaction_box(double choose2, int p, int q, int r ){
  int k;
  struct CaMKII *p1;
  double choose, choose1, v1, v2, v3, v4;

  choose = choose2;

  choose1 = choose;

  for(k=0; k<6; k++){
    p1 = box[p][q][r].head;
    while(p1 != NULL){
      choose1 = choose;
      choose += box[p][q][r].D_CaMKII[k];
      if((choose > random_num) && (choose1 <= random_num)){
	box[p][q][r].head = delete_CaMKII(box[p][q][r], p1);
	v1 = box[p][q][r].diffusion_lamda;
	v2 = box[p][q][r].chemical_lamda;
	box[p][q][r].diffusion_lamda = box_diffusion_lamda(p, q, r);
	box[p][q][r].chemical_lamda = voxel_chemical_lamda(box[p][q][r]);

	if(k == 0){	  
	  box[p-1][q][r].head = insert_CaMKII(box[p-1][q][r], p1);
	  v3 = box[p-1][q][r].diffusion_lamda;
	  v4 = box[p-1][q][r].chemical_lamda;
	  box[p-1][q][r].diffusion_lamda = box_diffusion_lamda(p-1, q, r);
	  box[p-1][q][r].chemical_lamda = voxel_chemical_lamda(box[p-1][q][r]);
	  reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p-1][q][r].diffusion_lamda +
	    box[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
	}else if(k == 1){
	  box[p][q-1][r].head = insert_CaMKII(box[p][q-1][r], p1);
	  v3 = box[p][q-1][r].diffusion_lamda;
	  v4 = box[p][q-1][r].chemical_lamda;
	  box[p][q-1][r].diffusion_lamda = box_diffusion_lamda(p, q-1, r);
	  box[p][q-1][r].chemical_lamda = voxel_chemical_lamda(box[p][q-1][r]);
	  reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q-1][r].diffusion_lamda +
	    box[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
	}else if(k == 2){
	  box[p+1][q][r].head = insert_CaMKII(box[p+1][q][r], p1);
	  v3 = box[p+1][q][r].diffusion_lamda;
	  v4 = box[p+1][q][r].chemical_lamda;
	  box[p+1][q][r].diffusion_lamda = box_diffusion_lamda(p+1, q, r);
	  box[p+1][q][r].chemical_lamda = voxel_chemical_lamda(box[p+1][q][r]);
	  reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p+1][q][r].diffusion_lamda +
	    box[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
	}else if(k == 3){
	  box[p][q+1][r].head = insert_CaMKII(box[p][q+1][r], p1);
	  v3 = box[p][q+1][r].diffusion_lamda;
	  v4 = box[p][q+1][r].chemical_lamda;
	  box[p][q+1][r].diffusion_lamda = box_diffusion_lamda(p, q+1, r);
	  box[p][q+1][r].chemical_lamda = voxel_chemical_lamda(box[p][q+1][r]);
	  reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q+1][r].diffusion_lamda +
	    box[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
	}else if(k == 4){
	  box[p][q][r-1].head = insert_CaMKII(box[p][q][r-1], p1);
	  v3 = box[p][q][r-1].diffusion_lamda;
	  v4 = box[p][q][r-1].chemical_lamda;
	  box[p][q][r-1].diffusion_lamda = box_diffusion_lamda(p, q, r-1);
	  box[p][q][r-1].chemical_lamda = voxel_chemical_lamda(box[p][q][r-1]);
	  reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r-1].diffusion_lamda +
	    box[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}else if(k == 5){
	  box[p][q][r+1].head = insert_CaMKII(box[p][q][r+1], p1);
	  v3 = box[p][q][r+1].diffusion_lamda;
	  v4 = box[p][q][r+1].chemical_lamda;
	  box[p][q][r+1].diffusion_lamda = box_diffusion_lamda(p, q, r+1);
	  box[p][q][r+1].chemical_lamda = voxel_chemical_lamda(box[p][q][r+1]);
	  reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r+1].diffusion_lamda +
	    box[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
        }
	goto loop1;
      }
      p1 = p1->next;
    }

    choose1 = choose;
    choose += box[p][q][r].Ca * box[p][q][r].D_Ca[k];
    if((choose > random_num) && (choose1 <= random_num)){
      box[p][q][r].Ca--;
      v1 = box[p][q][r].diffusion_lamda;
      v2 = box[p][q][r].chemical_lamda;
      box[p][q][r].diffusion_lamda = box_diffusion_lamda(p, q, r);
      box[p][q][r].chemical_lamda = voxel_chemical_lamda(box[p][q][r]);

      if(k == 0){
	box[p-1][q][r].Ca++;
	v3 = box[p-1][q][r].diffusion_lamda;
	v4 = box[p-1][q][r].chemical_lamda;
	box[p-1][q][r].diffusion_lamda = box_diffusion_lamda(p-1, q, r);
	box[p-1][q][r].chemical_lamda = voxel_chemical_lamda(box[p-1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p-1][q][r].diffusion_lamda +
	  box[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	box[p][q-1][r].Ca++;
	v3 = box[p][q-1][r].diffusion_lamda;
	v4 = box[p][q-1][r].chemical_lamda;
	box[p][q-1][r].diffusion_lamda = box_diffusion_lamda(p, q-1, r);
	box[p][q-1][r].chemical_lamda = voxel_chemical_lamda(box[p][q-1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q-1][r].diffusion_lamda +
	  box[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	box[p+1][q][r].Ca++;
	v3 = box[p+1][q][r].diffusion_lamda;
	v4 = box[p+1][q][r].chemical_lamda;
	box[p+1][q][r].diffusion_lamda = box_diffusion_lamda(p+1, q, r);
	box[p+1][q][r].chemical_lamda = voxel_chemical_lamda(box[p+1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p+1][q][r].diffusion_lamda +
	  box[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	box[p][q+1][r].Ca++;
	v3 = box[p][q+1][r].diffusion_lamda;
	v4 = box[p][q+1][r].chemical_lamda;
	box[p][q+1][r].diffusion_lamda = box_diffusion_lamda(p, q+1, r);
	box[p][q+1][r].chemical_lamda = voxel_chemical_lamda(box[p][q+1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q+1][r].diffusion_lamda +
	  box[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	box[p][q][r-1].Ca++;
	v3 = box[p][q][r-1].diffusion_lamda;
	v4 = box[p][q][r-1].chemical_lamda;
	box[p][q][r-1].diffusion_lamda = box_diffusion_lamda(p, q, r-1);
	box[p][q][r-1].chemical_lamda = voxel_chemical_lamda(box[p][q][r-1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r-1].diffusion_lamda +
	  box[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 5){
	box[p][q][r+1].Ca++;
	v3 = box[p][q][r+1].diffusion_lamda;
	v4 = box[p][q][r+1].chemical_lamda;
	box[p][q][r+1].diffusion_lamda = box_diffusion_lamda(p, q, r+1);
	box[p][q][r+1].chemical_lamda = voxel_chemical_lamda(box[p][q][r+1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r+1].diffusion_lamda +
	  box[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

    choose1 = choose;
    choose += box[p][q][r].CaM00 * box[p][q][r].D_CaM00[k];
    if((choose > random_num) && (choose1 <= random_num)){
      box[p][q][r].CaM00--;
      v1 = box[p][q][r].diffusion_lamda;
      v2 = box[p][q][r].chemical_lamda;
      box[p][q][r].diffusion_lamda = box_diffusion_lamda(p, q, r);
      box[p][q][r].chemical_lamda = voxel_chemical_lamda(box[p][q][r]);
      if(k == 0){
	box[p-1][q][r].CaM00++;
	v3 = box[p-1][q][r].diffusion_lamda;
	v4 = box[p-1][q][r].chemical_lamda;
	box[p-1][q][r].diffusion_lamda = box_diffusion_lamda(p-1, q, r);
	box[p-1][q][r].chemical_lamda = voxel_chemical_lamda(box[p-1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p-1][q][r].diffusion_lamda +
	  box[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	box[p][q-1][r].CaM00++;
	v3 = box[p][q-1][r].diffusion_lamda;
	v4 = box[p][q-1][r].chemical_lamda;
	box[p][q-1][r].diffusion_lamda = box_diffusion_lamda(p, q-1, r);
	box[p][q-1][r].chemical_lamda = voxel_chemical_lamda(box[p][q-1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q-1][r].diffusion_lamda +
	  box[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	box[p+1][q][r].CaM00++;
	v3 = box[p+1][q][r].diffusion_lamda;
	v4 = box[p+1][q][r].chemical_lamda;
	box[p+1][q][r].diffusion_lamda = box_diffusion_lamda(p+1, q, r);
	box[p+1][q][r].chemical_lamda = voxel_chemical_lamda(box[p+1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p+1][q][r].diffusion_lamda +
	  box[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	box[p][q+1][r].CaM00++;
	v3 = box[p][q+1][r].diffusion_lamda;
	v4 = box[p][q+1][r].chemical_lamda;
	box[p][q+1][r].diffusion_lamda = box_diffusion_lamda(p, q+1, r);
	box[p][q+1][r].chemical_lamda = voxel_chemical_lamda(box[p][q+1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q+1][r].diffusion_lamda +
	  box[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	box[p][q][r-1].CaM00++;
	v3 = box[p][q][r-1].diffusion_lamda;
	v4 = box[p][q][r-1].chemical_lamda;
	box[p][q][r-1].diffusion_lamda = box_diffusion_lamda(p, q, r-1);
	box[p][q][r-1].chemical_lamda = voxel_chemical_lamda(box[p][q][r-1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r-1].diffusion_lamda +
	  box[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 5){
	box[p][q][r+1].CaM00++;
	v3 = box[p][q][r+1].diffusion_lamda;
	v4 = box[p][q][r+1].chemical_lamda;
	box[p][q][r+1].diffusion_lamda = box_diffusion_lamda(p, q, r+1);
	box[p][q][r+1].chemical_lamda = voxel_chemical_lamda(box[p][q][r+1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r+1].diffusion_lamda +
	  box[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += box[p][q][r].CaM10 * box[p][q][r].D_CaM10[k];
    if((choose > random_num) && (choose1 <= random_num)){
      box[p][q][r].CaM10--;
      v1 = box[p][q][r].diffusion_lamda;
      v2 = box[p][q][r].chemical_lamda;
      box[p][q][r].diffusion_lamda = box_diffusion_lamda(p, q, r);
      box[p][q][r].chemical_lamda = voxel_chemical_lamda(box[p][q][r]);
      if(k == 0){
	box[p-1][q][r].CaM10++;
	v3 = box[p-1][q][r].diffusion_lamda;
	v4 = box[p-1][q][r].chemical_lamda;
	box[p-1][q][r].diffusion_lamda = box_diffusion_lamda(p-1, q, r);
	box[p-1][q][r].chemical_lamda = voxel_chemical_lamda(box[p-1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p-1][q][r].diffusion_lamda +
	  box[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	box[p][q-1][r].CaM10++;
	v3 = box[p][q-1][r].diffusion_lamda;
	v4 = box[p][q-1][r].chemical_lamda;
	box[p][q-1][r].diffusion_lamda = box_diffusion_lamda(p, q-1, r);
	box[p][q-1][r].chemical_lamda = voxel_chemical_lamda(box[p][q-1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q-1][r].diffusion_lamda +
	  box[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	box[p+1][q][r].CaM10++;
	v3 = box[p+1][q][r].diffusion_lamda;
	v4 = box[p+1][q][r].chemical_lamda;
	box[p+1][q][r].diffusion_lamda = box_diffusion_lamda(p+1, q, r);
	box[p+1][q][r].chemical_lamda = voxel_chemical_lamda(box[p+1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p+1][q][r].diffusion_lamda +
	  box[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	box[p][q+1][r].CaM10++;
	v3 = box[p][q+1][r].diffusion_lamda;
	v4 = box[p][q+1][r].chemical_lamda;
	box[p][q+1][r].diffusion_lamda = box_diffusion_lamda(p, q+1, r);
	box[p][q+1][r].chemical_lamda = voxel_chemical_lamda(box[p][q+1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q+1][r].diffusion_lamda +
	  box[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	box[p][q][r-1].CaM10++;
	v3 = box[p][q][r-1].diffusion_lamda;
	v4 = box[p][q][r-1].chemical_lamda;
	box[p][q][r-1].diffusion_lamda = box_diffusion_lamda(p, q, r-1);
	box[p][q][r-1].chemical_lamda = voxel_chemical_lamda(box[p][q][r-1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r-1].diffusion_lamda +
	  box[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 5){
	box[p][q][r+1].CaM10++;
	v3 = box[p][q][r+1].diffusion_lamda;
	v4 = box[p][q][r+1].chemical_lamda;
	box[p][q][r+1].diffusion_lamda = box_diffusion_lamda(p, q, r+1);
	box[p][q][r+1].chemical_lamda = voxel_chemical_lamda(box[p][q][r+1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r+1].diffusion_lamda +
	  box[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += box[p][q][r].CaM01 * box[p][q][r].D_CaM01[k];
    if((choose > random_num) && (choose1 <= random_num)){
      box[p][q][r].CaM01--;
      v1 = box[p][q][r].diffusion_lamda;
      v2 = box[p][q][r].chemical_lamda;
      box[p][q][r].diffusion_lamda = box_diffusion_lamda(p, q, r);
      box[p][q][r].chemical_lamda = voxel_chemical_lamda(box[p][q][r]);
      if(k == 0){
	box[p-1][q][r].CaM01++;
	v3 = box[p-1][q][r].diffusion_lamda;
	v4 = box[p-1][q][r].chemical_lamda;
	box[p-1][q][r].diffusion_lamda = box_diffusion_lamda(p-1, q, r);
	box[p-1][q][r].chemical_lamda = voxel_chemical_lamda(box[p-1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p-1][q][r].diffusion_lamda +
	  box[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	box[p][q-1][r].CaM01++;
	v3 = box[p][q-1][r].diffusion_lamda;
	v4 = box[p][q-1][r].chemical_lamda;
	box[p][q-1][r].diffusion_lamda = box_diffusion_lamda(p, q-1, r);
	box[p][q-1][r].chemical_lamda = voxel_chemical_lamda(box[p][q-1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q-1][r].diffusion_lamda +
	  box[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	box[p+1][q][r].CaM01++;
	v3 = box[p+1][q][r].diffusion_lamda;
	v4 = box[p+1][q][r].chemical_lamda;
	box[p+1][q][r].diffusion_lamda = box_diffusion_lamda(p+1, q, r);
	box[p+1][q][r].chemical_lamda = voxel_chemical_lamda(box[p+1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p+1][q][r].diffusion_lamda +
	  box[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	box[p][q+1][r].CaM01++;
	v3 = box[p][q+1][r].diffusion_lamda;
	v4 = box[p][q+1][r].chemical_lamda;
	box[p][q+1][r].diffusion_lamda = box_diffusion_lamda(p, q+1, r);
	box[p][q+1][r].chemical_lamda = voxel_chemical_lamda(box[p][q+1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q+1][r].diffusion_lamda +
	  box[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	box[p][q][r-1].CaM01++;
	v3 = box[p][q][r-1].diffusion_lamda;
	v4 = box[p][q][r-1].chemical_lamda;
	box[p][q][r-1].diffusion_lamda = box_diffusion_lamda(p, q, r-1);
	box[p][q][r-1].chemical_lamda = voxel_chemical_lamda(box[p][q][r-1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r-1].diffusion_lamda +
	  box[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 5){
	box[p][q][r+1].CaM01++;
	v3 = box[p][q][r+1].diffusion_lamda;
	v4 = box[p][q][r+1].chemical_lamda;
	box[p][q][r+1].diffusion_lamda = box_diffusion_lamda(p, q, r+1);
	box[p][q][r+1].chemical_lamda = voxel_chemical_lamda(box[p][q][r+1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r+1].diffusion_lamda +
	  box[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += box[p][q][r].CaM11 * box[p][q][r].D_CaM11[k];
    if((choose > random_num) && (choose1 <= random_num)){
      box[p][q][r].CaM11--;
      v1 = box[p][q][r].diffusion_lamda;
      v2 = box[p][q][r].chemical_lamda;
      box[p][q][r].diffusion_lamda = box_diffusion_lamda(p, q, r);
      box[p][q][r].chemical_lamda = voxel_chemical_lamda(box[p][q][r]);
      if(k == 0){
	box[p-1][q][r].CaM11++;
	v3 = box[p-1][q][r].diffusion_lamda;
	v4 = box[p-1][q][r].chemical_lamda;
	box[p-1][q][r].diffusion_lamda = box_diffusion_lamda(p-1, q, r);
	box[p-1][q][r].chemical_lamda = voxel_chemical_lamda(box[p-1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p-1][q][r].diffusion_lamda +
	  box[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	box[p][q-1][r].CaM11++;
	v3 = box[p][q-1][r].diffusion_lamda;
	v4 = box[p][q-1][r].chemical_lamda;
	box[p][q-1][r].diffusion_lamda = box_diffusion_lamda(p, q-1, r);
	box[p][q-1][r].chemical_lamda = voxel_chemical_lamda(box[p][q-1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q-1][r].diffusion_lamda +
	  box[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	box[p+1][q][r].CaM11++;
	v3 = box[p+1][q][r].diffusion_lamda;
	v4 = box[p+1][q][r].chemical_lamda;
	box[p+1][q][r].diffusion_lamda = box_diffusion_lamda(p+1, q, r);
	box[p+1][q][r].chemical_lamda = voxel_chemical_lamda(box[p+1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p+1][q][r].diffusion_lamda +
	  box[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	box[p][q+1][r].CaM11++;
	v3 = box[p][q+1][r].diffusion_lamda;
	v4 = box[p][q+1][r].chemical_lamda;
	box[p][q+1][r].diffusion_lamda = box_diffusion_lamda(p, q+1, r);
	box[p][q+1][r].chemical_lamda = voxel_chemical_lamda(box[p][q+1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q+1][r].diffusion_lamda +
	  box[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	box[p][q][r-1].CaM11++;
	v3 = box[p][q][r-1].diffusion_lamda;
	v4 = box[p][q][r-1].chemical_lamda;
	box[p][q][r-1].diffusion_lamda = box_diffusion_lamda(p, q, r-1);
	box[p][q][r-1].chemical_lamda = voxel_chemical_lamda(box[p][q][r-1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r-1].diffusion_lamda +
	  box[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 5){
	box[p][q][r+1].CaM11++;
	v3 = box[p][q][r+1].diffusion_lamda;
	v4 = box[p][q][r+1].chemical_lamda;
	box[p][q][r+1].diffusion_lamda = box_diffusion_lamda(p, q, r+1);
	box[p][q][r+1].chemical_lamda = voxel_chemical_lamda(box[p][q][r+1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r+1].diffusion_lamda +
	  box[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += box[p][q][r].CaM20 * box[p][q][r].D_CaM20[k];
    if((choose > random_num) && (choose1 <= random_num)){
      box[p][q][r].CaM20--;
      v1 = box[p][q][r].diffusion_lamda;
      v2 = box[p][q][r].chemical_lamda;
      box[p][q][r].diffusion_lamda = box_diffusion_lamda(p, q, r);
      box[p][q][r].chemical_lamda = voxel_chemical_lamda(box[p][q][r]);
      if(k == 0){
	box[p-1][q][r].CaM20++;
	v3 = box[p-1][q][r].diffusion_lamda;
	v4 = box[p-1][q][r].chemical_lamda;
	box[p-1][q][r].diffusion_lamda = box_diffusion_lamda(p-1, q, r);
	box[p-1][q][r].chemical_lamda = voxel_chemical_lamda(box[p-1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p-1][q][r].diffusion_lamda +
	  box[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	box[p][q-1][r].CaM20++;
	v3 = box[p][q-1][r].diffusion_lamda;
	v4 = box[p][q-1][r].chemical_lamda;
	box[p][q-1][r].diffusion_lamda = box_diffusion_lamda(p, q-1, r);
	box[p][q-1][r].chemical_lamda = voxel_chemical_lamda(box[p][q-1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q-1][r].diffusion_lamda +
	  box[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	box[p+1][q][r].CaM20++;
	v3 = box[p+1][q][r].diffusion_lamda;
	v4 = box[p+1][q][r].chemical_lamda;
	box[p+1][q][r].diffusion_lamda = box_diffusion_lamda(p+1, q, r);
	box[p+1][q][r].chemical_lamda = voxel_chemical_lamda(box[p+1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p+1][q][r].diffusion_lamda +
	  box[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	box[p][q+1][r].CaM20++;
	v3 = box[p][q+1][r].diffusion_lamda;
	v4 = box[p][q+1][r].chemical_lamda;
	box[p][q+1][r].diffusion_lamda = box_diffusion_lamda(p, q+1, r);
	box[p][q+1][r].chemical_lamda = voxel_chemical_lamda(box[p][q+1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q+1][r].diffusion_lamda +
	  box[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	box[p][q][r-1].CaM20++;
	v3 = box[p][q][r-1].diffusion_lamda;
	v4 = box[p][q][r-1].chemical_lamda;
	box[p][q][r-1].diffusion_lamda = box_diffusion_lamda(p, q, r-1);
	box[p][q][r-1].chemical_lamda = voxel_chemical_lamda(box[p][q][r-1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r-1].diffusion_lamda +
	  box[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 5){
	box[p][q][r+1].CaM20++;
	v3 = box[p][q][r+1].diffusion_lamda;
	v4 = box[p][q][r+1].chemical_lamda;
	box[p][q][r+1].diffusion_lamda = box_diffusion_lamda(p, q, r+1);
	box[p][q][r+1].chemical_lamda = voxel_chemical_lamda(box[p][q][r+1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r+1].diffusion_lamda +
	  box[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += box[p][q][r].CaM02 * box[p][q][r].D_CaM02[k];
    if((choose > random_num) && (choose1 <= random_num)){
      box[p][q][r].CaM02--;
      v1 = box[p][q][r].diffusion_lamda;
      v2 = box[p][q][r].chemical_lamda;
      box[p][q][r].diffusion_lamda = box_diffusion_lamda(p, q, r);
      box[p][q][r].chemical_lamda = voxel_chemical_lamda(box[p][q][r]);
      if(k == 0){
	box[p-1][q][r].CaM02++;
	v3 = box[p-1][q][r].diffusion_lamda;
	v4 = box[p-1][q][r].chemical_lamda;
	box[p-1][q][r].diffusion_lamda = box_diffusion_lamda(p-1, q, r);
	box[p-1][q][r].chemical_lamda = voxel_chemical_lamda(box[p-1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p-1][q][r].diffusion_lamda +
	  box[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	box[p][q-1][r].CaM02++;
	v3 = box[p][q-1][r].diffusion_lamda;
	v4 = box[p][q-1][r].chemical_lamda;
	box[p][q-1][r].diffusion_lamda = box_diffusion_lamda(p, q-1, r);
	box[p][q-1][r].chemical_lamda = voxel_chemical_lamda(box[p][q-1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q-1][r].diffusion_lamda +
	  box[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	box[p+1][q][r].CaM02++;
	v3 = box[p+1][q][r].diffusion_lamda;
	v4 = box[p+1][q][r].chemical_lamda;
	box[p+1][q][r].diffusion_lamda = box_diffusion_lamda(p+1, q, r);
	box[p+1][q][r].chemical_lamda = voxel_chemical_lamda(box[p+1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p+1][q][r].diffusion_lamda +
	  box[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	box[p][q+1][r].CaM02++;
	v3 = box[p][q+1][r].diffusion_lamda;
	v4 = box[p][q+1][r].chemical_lamda;
	box[p][q+1][r].diffusion_lamda = box_diffusion_lamda(p, q+1, r);
	box[p][q+1][r].chemical_lamda = voxel_chemical_lamda(box[p][q+1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q+1][r].diffusion_lamda +
	  box[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	box[p][q][r-1].CaM02++;
	v3 = box[p][q][r-1].diffusion_lamda;
	v4 = box[p][q][r-1].chemical_lamda;
	box[p][q][r-1].diffusion_lamda = box_diffusion_lamda(p, q, r-1);
	box[p][q][r-1].chemical_lamda = voxel_chemical_lamda(box[p][q][r-1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r-1].diffusion_lamda +
	  box[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 5){
	box[p][q][r+1].CaM02++;
	v3 = box[p][q][r+1].diffusion_lamda;
	v4 = box[p][q][r+1].chemical_lamda;
	box[p][q][r+1].diffusion_lamda = box_diffusion_lamda(p, q, r+1);
	box[p][q][r+1].chemical_lamda = voxel_chemical_lamda(box[p][q][r+1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r+1].diffusion_lamda +
	  box[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += box[p][q][r].CaM21 * box[p][q][r].D_CaM21[k];
    if((choose > random_num) && (choose1 <= random_num)){
      box[p][q][r].CaM21--;
      v1 = box[p][q][r].diffusion_lamda;
      v2 = box[p][q][r].chemical_lamda;
      box[p][q][r].diffusion_lamda = box_diffusion_lamda(p, q, r);
      box[p][q][r].chemical_lamda = voxel_chemical_lamda(box[p][q][r]);
      if(k == 0){
	box[p-1][q][r].CaM21++;
	v3 = box[p-1][q][r].diffusion_lamda;
	v4 = box[p-1][q][r].chemical_lamda;
	box[p-1][q][r].diffusion_lamda = box_diffusion_lamda(p-1, q, r);
	box[p-1][q][r].chemical_lamda = voxel_chemical_lamda(box[p-1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p-1][q][r].diffusion_lamda +
	  box[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	box[p][q-1][r].CaM21++;
	v3 = box[p][q-1][r].diffusion_lamda;
	v4 = box[p][q-1][r].chemical_lamda;
	box[p][q-1][r].diffusion_lamda = box_diffusion_lamda(p, q-1, r);
	box[p][q-1][r].chemical_lamda = voxel_chemical_lamda(box[p][q-1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q-1][r].diffusion_lamda +
	  box[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	box[p+1][q][r].CaM21++;
	v3 = box[p+1][q][r].diffusion_lamda;
	v4 = box[p+1][q][r].chemical_lamda;
	box[p+1][q][r].diffusion_lamda = box_diffusion_lamda(p+1, q, r);
	box[p+1][q][r].chemical_lamda = voxel_chemical_lamda(box[p+1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p+1][q][r].diffusion_lamda +
	  box[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	box[p][q+1][r].CaM21++;
	v3 = box[p][q+1][r].diffusion_lamda;
	v4 = box[p][q+1][r].chemical_lamda;
	box[p][q+1][r].diffusion_lamda = box_diffusion_lamda(p, q+1, r);
	box[p][q+1][r].chemical_lamda = voxel_chemical_lamda(box[p][q+1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q+1][r].diffusion_lamda +
	  box[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	box[p][q][r-1].CaM21++;
	v3 = box[p][q][r-1].diffusion_lamda;
	v4 = box[p][q][r-1].chemical_lamda;
	box[p][q][r-1].diffusion_lamda = box_diffusion_lamda(p, q, r-1);
	box[p][q][r-1].chemical_lamda = voxel_chemical_lamda(box[p][q][r-1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r-1].diffusion_lamda +
	  box[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 5){
	box[p][q][r+1].CaM21++;
	v3 = box[p][q][r+1].diffusion_lamda;
	v4 = box[p][q][r+1].chemical_lamda;
	box[p][q][r+1].diffusion_lamda = box_diffusion_lamda(p, q, r+1);
	box[p][q][r+1].chemical_lamda = voxel_chemical_lamda(box[p][q][r+1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r+1].diffusion_lamda +
	  box[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += box[p][q][r].CaM12 * box[p][q][r].D_CaM12[k];
    if((choose > random_num) && (choose1 <= random_num)){
      box[p][q][r].CaM12--;
      v1 = box[p][q][r].diffusion_lamda;
      v2 = box[p][q][r].chemical_lamda;
      box[p][q][r].diffusion_lamda = box_diffusion_lamda(p, q, r);
      box[p][q][r].chemical_lamda = voxel_chemical_lamda(box[p][q][r]);
      if(k == 0){
	box[p-1][q][r].CaM12++;
	v3 = box[p-1][q][r].diffusion_lamda;
	v4 = box[p-1][q][r].chemical_lamda;
	box[p-1][q][r].diffusion_lamda = box_diffusion_lamda(p-1, q, r);
	box[p-1][q][r].chemical_lamda = voxel_chemical_lamda(box[p-1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p-1][q][r].diffusion_lamda +
	  box[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	box[p][q-1][r].CaM12++;
	v3 = box[p][q-1][r].diffusion_lamda;
	v4 = box[p][q-1][r].chemical_lamda;
	box[p][q-1][r].diffusion_lamda = box_diffusion_lamda(p, q-1, r);
	box[p][q-1][r].chemical_lamda = voxel_chemical_lamda(box[p][q-1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q-1][r].diffusion_lamda +
	  box[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	box[p+1][q][r].CaM12++;
	v3 = box[p+1][q][r].diffusion_lamda;
	v4 = box[p+1][q][r].chemical_lamda;
	box[p+1][q][r].diffusion_lamda = box_diffusion_lamda(p+1, q, r);
	box[p+1][q][r].chemical_lamda = voxel_chemical_lamda(box[p+1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p+1][q][r].diffusion_lamda +
	  box[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	box[p][q+1][r].CaM12++;
	v3 = box[p][q+1][r].diffusion_lamda;
	v4 = box[p][q+1][r].chemical_lamda;
	box[p][q+1][r].diffusion_lamda = box_diffusion_lamda(p, q+1, r);
	box[p][q+1][r].chemical_lamda = voxel_chemical_lamda(box[p][q+1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q+1][r].diffusion_lamda +
	  box[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	box[p][q][r-1].CaM12++;
	v3 = box[p][q][r-1].diffusion_lamda;
	v4 = box[p][q][r-1].chemical_lamda;
	box[p][q][r-1].diffusion_lamda = box_diffusion_lamda(p, q, r-1);
	box[p][q][r-1].chemical_lamda = voxel_chemical_lamda(box[p][q][r-1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r-1].diffusion_lamda +
	  box[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 5){
	box[p][q][r+1].CaM12++;
	v3 = box[p][q][r+1].diffusion_lamda;
	v4 = box[p][q][r+1].chemical_lamda;
	box[p][q][r+1].diffusion_lamda = box_diffusion_lamda(p, q, r+1);
	box[p][q][r+1].chemical_lamda = voxel_chemical_lamda(box[p][q][r+1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r+1].diffusion_lamda +
	  box[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += box[p][q][r].CaM22 * box[p][q][r].D_CaM22[k];
    if((choose > random_num) && (choose1 <= random_num)){
      box[p][q][r].CaM22--;
      v1 = box[p][q][r].diffusion_lamda;
      v2 = box[p][q][r].chemical_lamda;
      box[p][q][r].diffusion_lamda = box_diffusion_lamda(p, q, r);
      box[p][q][r].chemical_lamda = voxel_chemical_lamda(box[p][q][r]);
      if(k == 0){
	box[p-1][q][r].CaM22++;
	v3 = box[p-1][q][r].diffusion_lamda;
	v4 = box[p-1][q][r].chemical_lamda;
	box[p-1][q][r].diffusion_lamda = box_diffusion_lamda(p-1, q, r);
	box[p-1][q][r].chemical_lamda = voxel_chemical_lamda(box[p-1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p-1][q][r].diffusion_lamda +
	  box[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	box[p][q-1][r].CaM22++;
	v3 = box[p][q-1][r].diffusion_lamda;
	v4 = box[p][q-1][r].chemical_lamda;
	box[p][q-1][r].diffusion_lamda = box_diffusion_lamda(p, q-1, r);
	box[p][q-1][r].chemical_lamda = voxel_chemical_lamda(box[p][q-1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q-1][r].diffusion_lamda +
	  box[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	box[p+1][q][r].CaM22++;
	v3 = box[p+1][q][r].diffusion_lamda;
	v4 = box[p+1][q][r].chemical_lamda;
	box[p+1][q][r].diffusion_lamda = box_diffusion_lamda(p+1, q, r);
	box[p+1][q][r].chemical_lamda = voxel_chemical_lamda(box[p+1][q][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p+1][q][r].diffusion_lamda +
	  box[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	box[p][q+1][r].CaM22++;
	v3 = box[p][q+1][r].diffusion_lamda;
	v4 = box[p][q+1][r].chemical_lamda;
	box[p][q+1][r].diffusion_lamda = box_diffusion_lamda(p, q+1, r);
	box[p][q+1][r].chemical_lamda = voxel_chemical_lamda(box[p][q+1][r]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q+1][r].diffusion_lamda +
	  box[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	box[p][q][r-1].CaM22++;
	v3 = box[p][q][r-1].diffusion_lamda;
	v4 = box[p][q][r-1].chemical_lamda;
	box[p][q][r-1].diffusion_lamda = box_diffusion_lamda(p, q, r-1);
	box[p][q][r-1].chemical_lamda = voxel_chemical_lamda(box[p][q][r-1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r-1].diffusion_lamda +
	  box[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 5){
	box[p][q][r+1].CaM22++;
	v3 = box[p][q][r+1].diffusion_lamda;
	v4 = box[p][q][r+1].chemical_lamda;
	box[p][q][r+1].diffusion_lamda = box_diffusion_lamda(p, q, r+1);
	box[p][q][r+1].chemical_lamda = voxel_chemical_lamda(box[p][q][r+1]);
	reaction_lamda = reaction_lamda + box[p][q][r].diffusion_lamda + box[p][q][r].chemical_lamda + box[p][q][r+1].diffusion_lamda +
	  box[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }
 } 
//  not sure about the {} here
  x = p;
  y = q;
  z = r;
  chemical_reaction(&box[p][q][r], choose, 1);
 loop1: k = 0;
// cout<< " r_l ch " << reaction_lamda  << " " << choose << endl;    
}


void chemical_reaction(struct voxel *point, double choose2, int index){
  int i, j;
  double choose, choose1, v1, v2;
  struct CaMKII *p1;

  choose = choose2;

// cout<< " ch choice " << choose << " " << random_num << endl;    
  choose1 = choose;
  choose += point->CaM00 * point->Ca * point->CaM00_to_CaM10;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM00--;
    point->Ca--;
    point->CaM10++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM10 * point->CaM10_to_CaM00;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM00++;
    point->Ca++;
    point->CaM10--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM00 * point->Ca * point->CaM00_to_CaM01;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM00--;
    point->Ca--;
    point->CaM01++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM01 * point->CaM01_to_CaM00;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM00++;
    point->Ca++;
    point->CaM01--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM10 * point->Ca * point->CaM10_to_CaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM10--;
    point->Ca--;
    point->CaM11++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM11 * point->CaM11_to_CaM10;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM10++;
    point->Ca++;
    point->CaM11--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM10 * point->Ca * point->CaM10_to_CaM20;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM10--;
    point->Ca--;
    point->CaM20++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM20 * point->CaM20_to_CaM10;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM10++;
    point->Ca++;
    point->CaM20--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM01 * point->Ca * point->CaM01_to_CaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM01--;
    point->Ca--;
    point->CaM11++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM11 * point->CaM11_to_CaM01;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM01++;
    point->Ca++;
    point->CaM11--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM01 * point->Ca * point->CaM01_to_CaM02;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM01--;
    point->Ca--;
    point->CaM02++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM02 * point->CaM02_to_CaM01;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM01++;
    point->Ca++;
    point->CaM02--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM11 * point->Ca * point->CaM11_to_CaM21;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM11--;
    point->Ca--;
    point->CaM21++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM21 * point->CaM21_to_CaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM11++;
    point->Ca++;
    point->CaM21--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM11 * point->Ca * point->CaM11_to_CaM12;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM11--;
    point->Ca--;
    point->CaM12++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM12 * point->CaM12_to_CaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM11++;
    point->Ca++;
    point->CaM12--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM20 * point->Ca * point->CaM20_to_CaM21;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM20--;
    point->Ca--;
    point->CaM21++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM21 * point->CaM21_to_CaM20;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM20++;
    point->Ca++;
    point->CaM21--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM02 * point->Ca * point->CaM02_to_CaM12;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM02--;
    point->Ca--;
    point->CaM12++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM12 * point->CaM12_to_CaM02;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM02++;
    point->Ca++;
    point->CaM12--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM21 * point->Ca * point->CaM21_to_CaM22;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM21--;
    point->Ca--;
    point->CaM22++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM22 * point->CaM22_to_CaM21;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM21++;
    point->Ca++;
    point->CaM22--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM12 * point->Ca * point->CaM12_to_CaM22;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM12--;
    point->Ca--;
    point->CaM22++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM22 * point->CaM22_to_CaM12;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM12++;
    point->Ca++;
    point->CaM22--;
    goto loop5;
  }
  
//------------------------------------------------------------------ 

  p1 = point->head;
  while(p1 != NULL){
    for(i=0; i<1; i++){
      for(j=0; j<6; j++){

	if(p1->subunit[i][j] ==0){
		
	  choose1 = choose;
	  choose += point->CaM00 * point->CaM00_to_CaMKIICaM00;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM00--;
	    p1->subunit[i][j] = 1;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM10 * point->CaM10_to_CaMKIICaM10;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM10--;
	    p1->subunit[i][j] = 2;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM01 * point->CaM01_to_CaMKIICaM01;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM01--;
	    p1->subunit[i][j] = 3;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM11 * point->CaM11_to_CaMKIICaM11;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM11--;
	    p1->subunit[i][j] = 4;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM20 * point->CaM20_to_CaMKIICaM20;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM20--;
	    p1->subunit[i][j] = 5;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM02 * point->CaM02_to_CaMKIICaM02;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM02--;
	    p1->subunit[i][j] = 6;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM21 * point->CaM21_to_CaMKIICaM21;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM21--;
	    p1->subunit[i][j] = 7;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM12 * point->CaM12_to_CaMKIICaM12;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM12--;
	    p1->subunit[i][j] = 8;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM22 * point->CaM22_to_CaMKIICaM22;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM22--;
	    p1->subunit[i][j] = 9;
	    goto loop5;
	  }
 
	}else if(p1->subunit[i][j] == 1){

	  choose1 = choose;
	  choose += point->CaMKIICaM00_to_CaM00;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM00++;
	    p1->subunit[i][j] = 0;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM00_to_CaMKIICaM10;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 2;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM00_to_CaMKIICaM01;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 3;
	    goto loop5;
	  }

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
     	 choose += point->CaMKIICaM00_to_Trapped00;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 10;
	    goto loop5;
	  }
/*
	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
         (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
         choose += 0.4 * point->CaMKIICaM00_to_Trapped00;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 10;
	    goto loop5;
	  }
*/
	}else if(p1->subunit[i][j] == 2){

	  choose1 = choose;
	  choose += point->CaMKIICaM10_to_CaMKIICaM00;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 1;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM10_to_CaM10;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM10++;
	    p1->subunit[i][j] = 0;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM10_to_CaMKIICaM11;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 4;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM10_to_CaMKIICaM20;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 5;
	    goto loop5;
	  }
	  
          choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
     	 choose += point->CaMKIICaM10_to_Trapped10;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 11;
	    goto loop5;
	  }

/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
	  choose1 = choose;
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
         (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
  
    	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20)) // ||
  	     choose += 0.4 * point->CaMKIICaM10_to_Trapped10;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 11;
	    goto loop5;
	  }
 */ 
	}else if(p1->subunit[i][j] == 3){

	  choose1 = choose;
	  choose += point->CaMKIICaM01_to_CaMKIICaM00;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 1;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM01_to_CaM01;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM01++;
	    p1->subunit[i][j] = 0;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM01_to_CaMKIICaM11;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 4;
	    goto loop5;
	  }

      choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM01_to_CaMKIICaM02;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 6;
	    goto loop5;
	  }

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
     	 choose += point->CaMKIICaM01_to_Trapped01;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 12;
	    goto loop5;
	  }

/*	  choose1 = choose;
  	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
         (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
  	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20))  //||
	     choose += 0.4 * point->CaMKIICaM01_to_Trapped01;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 12;
	    goto loop5;
	  }
*/
  
	}else if(p1->subunit[i][j] == 4){
	  
	  choose1 = choose;
	  choose += point->CaMKIICaM11_to_CaMKIICaM01;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 3;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM11_to_CaMKIICaM10;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 2;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose1 = point->CaMKIICaM11_to_CaM11;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM11++;
	    p1->subunit[i][j] = 0;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose1 = point->Ca * point->CaMKIICaM11_to_CaMKIICaM21;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 7;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose1 = point->Ca * point->CaMKIICaM11_to_CaMKIICaM12;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 8;
	    goto loop5;
	  }

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
         choose += point->CaMKIICaM11_to_Trapped11;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 13;
	    goto loop5;
	  }

/*
	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	     choose += 0.4 * point->CaMKIICaM11_to_Trapped11;
//
//	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20))
// ||
//                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
//	     choose += 0.4 * point->CaMKIICaM11_to_Trapped11;
//
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 13;
	    goto loop5;
	  }
*/
	  
	}else if(p1->subunit[i][j] == 5){

	  choose1 = choose;
	  choose += point->CaMKIICaM20_to_CaMKIICaM10;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 2;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM20_to_CaMKIICaM21;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 7;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM20_to_CaM20;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM20++;
	    p1->subunit[i][j] = 0;
	    goto loop5;
	  }

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
        choose += point->CaMKIICaM20_to_Trapped20;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 14;
	    goto loop5;
	  }

/*
	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM20_to_Trapped20;
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20))
//  ||
//                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM20_to_Trapped20;
//
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 14;
	    goto loop5;
	  }
*/

	}else if(p1->subunit[i][j] == 6){
	
	  choose1 = choose;
	  choose += point->CaMKIICaM02_to_CaMKIICaM01;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 3;
	  point->Ca++;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM02_to_CaMKIICaM12;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 8;
	  point->Ca--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM02_to_CaM02;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 0;
	  point->CaM02++;
	  goto loop5;
	  }

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
        choose += point->CaMKIICaM02_to_Trapped02;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 15;
	    goto loop5;
	  }
/*

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM02_to_Trapped02;
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20))
//  ||
//                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM02_to_Trapped02;
//
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 15;
	    goto loop5;
	  }
*/
	  
	}else if(p1->subunit[i][j] == 7){
		
	  choose1 = choose;
	  choose += point->CaMKIICaM21_to_CaMKIICaM11;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 4;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM21_to_CaMKIICaM20;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 5;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM21_to_CaMKIICaM22;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 9;
	    goto loop5;
	  }

      choose1 = choose;
	  choose += point->CaMKIICaM21_to_CaM21;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM21++;
	    p1->subunit[i][j] = 0;
	    goto loop5;
	  }

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
        choose += point->CaMKIICaM21_to_Trapped21;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 16;
	    goto loop5;
	  }

/*
	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM21_to_Trapped21;
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20))
//  ||
//                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM21_to_Trapped21;
//
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 16;
	    goto loop5;
	  }
*/
	}else if(p1->subunit[i][j] == 8){

	  choose1 = choose;
	  choose += point->CaMKIICaM12_to_CaMKIICaM11;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 4;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM12_to_CaMKIICaM02;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 6;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM12_to_CaMKIICaM22;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 9;
	    goto loop5;
	  }

      choose1 = choose;
	  choose += point->CaMKIICaM12_to_CaM12;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM12++;
	    p1->subunit[i][j] = 0;
	    goto loop5;
	  }

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
        choose += point->CaMKIICaM12_to_Trapped12;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 17;
	    goto loop5;
	  }

/*
	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM12_to_Trapped12;
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20))
// ||
//                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM12_to_Trapped12;
//
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 17;
	    goto loop5;
	  }
*/
	}else if(p1->subunit[i][j] == 9){
      
	  choose1 = choose;
	  choose += point->CaMKIICaM22_to_CaMKIICaM21;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 7;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM22_to_CaMKIICaM12;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 8;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM22_to_CaM22;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM22++;
	    p1->subunit[i][j] = 0;
	    goto loop5;
	  }

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
        choose += point->CaMKIICaM22_to_Trapped22;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 18;
	    goto loop5;
	  }

/*
	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM22_to_Trapped22;
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20))
//  ||
//                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM22_to_Trapped22;
//
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 18;
	    goto loop5;
	  }
*/
	}else if(p1->subunit[i][j] == 10){

      choose1 = choose;
	  choose += point->Trapped00_to_CaMKIICaM00;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 1;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped00_to_Trapped01;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 12;
		 point->Ca--;
	     goto loop5;
	   }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped00_to_Trapped10;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 11;
		 point->Ca--;
	     goto loop5;
	   }

	  choose1 = choose;
	  choose += point->Trapped00_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM00++;
	     goto loop5;
	   }

      }else if(p1->subunit[i][j] == 11){
      choose1 = choose;
	  choose += point->Trapped10_to_CaMKIICaM10;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 2;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped10_to_Trapped00;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 10;
	  goto loop5;
	  }

      choose1 = choose;
	  choose += point->Ca * point->Trapped10_to_Trapped11;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 13;
		 point->Ca--;
	     goto loop5;
	   }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped10_to_Trapped20;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 14;
		 point->Ca--;
	     goto loop5;
	   }

      choose1 = choose;
	  choose += point->Trapped10_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM10++;
	     goto loop5;
	   }

	  }else if(p1->subunit[i][j] == 12){
      choose1 = choose;
	  choose += point->Trapped01_to_CaMKIICaM01;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 3;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped01_to_Trapped00;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 10;
	  goto loop5;
	  }

      choose1 = choose;
	  choose += point->Ca * point->Trapped01_to_Trapped11;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 13;
		 point->Ca--;
	     goto loop5;
	   }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped01_to_Trapped02;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 15;
		 point->Ca--;
	     goto loop5;
	   }

      choose1 = choose;
	  choose += point->Trapped01_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM01++;
	     goto loop5;
	   }

	  }else if(p1->subunit[i][j] == 13){
      choose1 = choose;
	  choose += point->Trapped11_to_CaMKIICaM11;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 4;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped11_to_Trapped10;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 11;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped11_to_Trapped01;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 12;
	  goto loop5;
	  }

      choose1 = choose;
	  choose += point->Ca * point->Trapped11_to_Trapped21;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 16;
		 point->Ca--;
	     goto loop5;
	   }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped11_to_Trapped12;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 17;
		 point->Ca--;
	     goto loop5;
	   }

	  choose1 = choose;
	  choose += point->Trapped11_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM11++;
	     goto loop5;
	   }

	  }else if(p1->subunit[i][j] == 14){
      choose1 = choose;
	  choose += point->Trapped20_to_CaMKIICaM20;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 5;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped20_to_Trapped10;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 11;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped20_to_Trapped21;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 16;
		 point->Ca--;
	     goto loop5;
	   }
	  
	  choose1 = choose;
	  choose += point->Trapped20_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM20++;
	     goto loop5;
	   }

	  }else if(p1->subunit[i][j] == 15){
      choose1 = choose;
	  choose += point->Trapped02_to_CaMKIICaM02;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 6;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped02_to_Trapped01;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 12;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped02_to_Trapped12;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 17;
		 point->Ca--;
	     goto loop5;
	   }
	  
	  choose1 = choose;
	  choose += point->Trapped02_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM02++;
	     goto loop5;
	   }

	  }else if(p1->subunit[i][j] == 16){
      choose1 = choose;
	  choose += point->Trapped21_to_CaMKIICaM21;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 7;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped21_to_Trapped20;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 14;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped21_to_Trapped11;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 13;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped21_to_Trapped22;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 18;
		 point->Ca--;
	     goto loop5;
	   }
	  
	  choose1 = choose;
	  choose += point->Trapped21_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM21++;
	     goto loop5;
	   }

	  }else if(p1->subunit[i][j] == 17){
      choose1 = choose;
	  choose += point->Trapped12_to_CaMKIICaM12;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 8;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped12_to_Trapped02;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 15;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped12_to_Trapped11;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 13;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped12_to_Trapped22;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 18;
		 point->Ca--;
	     goto loop5;
	   }

	  choose1 = choose;
	  choose += point->Trapped21_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM21++;
	     goto loop5;
	   }

	  }else if(p1->subunit[i][j] == 18){
      choose1 = choose;
	  choose += point->Trapped22_to_CaMKIICaM22;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 9;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped22_to_Trapped12;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 17;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped22_to_Trapped21;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 16;
	  goto loop5;
	  }
	  
	  choose1 = choose;
	  choose += point->Trapped22_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM22++;
	     goto loop5;
	   }

	  }else if(p1->subunit[i][j] == 19){
          choose1 = choose;
	  choose += point->CaM00 * point->Auton_to_Trapped00;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 10;
	  point->CaM00--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM10 * point->Auton_to_Trapped10;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 11;
	  point->CaM10--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM01 * point->Auton_to_Trapped01;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 12;
	  point->CaM01--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM11 * point->Auton_to_Trapped11;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 13;
	  point->CaM11--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM20 * point->Auton_to_Trapped20;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 14;
	  point->CaM20--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM02 * point->Auton_to_Trapped02;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 15;
	  point->CaM02--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM21 * point->Auton_to_Trapped21;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 16;
	  point->CaM21--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM12 * point->Auton_to_Trapped12;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 17;
	  point->CaM12--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM22 * point->Auton_to_Trapped22;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 18;
	  point->CaM22--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Auton_to_CaMKII;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 0;
	  goto loop5;
	  }

          choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) ||
		 (p1->subunit[i][(j+5)%6] == 2) || (p1->subunit[i][(j+5)%6] == 3) ||
		 (p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
		 (p1->subunit[i][(j+5)%6] == 18))
   	    choose += point->Auton_to_Capped;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 20;
	     goto loop5;
	   }

/*
	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20)
                                              || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))

	    choose += 0.4 * point->Auton_to_Capped;
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20))
//  ||
//                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->Auton_to_Capped;
//
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 20;
	    goto loop5;
	  }
*/

	  }else{
          choose1 = choose;
	  choose += point->Capped_to_Auton;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 19;
	  goto loop5;
	  }
	  }
	  }
    }
    p1 = p1->next;
  }
// cout<< " no choice " << choose << " " << random_num << endl;    

 loop5: v1 = point->diffusion_lamda;
  v2 = point->chemical_lamda;
  if(index == 1){
    point->diffusion_lamda = box_diffusion_lamda(x, y, z);
    point->chemical_lamda = voxel_chemical_lamda(box[x][y][z]);
  }
  reaction_lamda += point->diffusion_lamda + point->chemical_lamda - v1 - v2;
}   


struct CaMKII * delete_CaMKII(struct voxel grid, struct CaMKII *point){
  struct CaMKII *p1, *p2, *head;

  head = grid.head;
  if(head == NULL) cout << "empty list" << endl;
  p1 = head;
  while(point->index != p1->index && p1->next != NULL){
    p2 = p1;
    p1 = p1->next;
  }
  if(point->index == p1->index){
    if(p1 == head) head = p1->next;
    else p2->next = p1->next;
  }

  return head;
}

struct CaMKII * insert_CaMKII(struct voxel grid, struct CaMKII *point){
 
  point->next = grid.head;
  return point;
}
	
double ran2(long *idum){
  long const IM1 = 2147483563;
  long const IM2 = 2147483399;
  double const AM = 1.0/IM1;
  long const IMM1 = IM1-1;
  int const IA1 = 40014;
  int const IA2 = 40692;
  int const IQ1 = 53668;
  int const IQ2 = 52774;
  int const IR1 = 12211;
  int const IR2 = 3791;
  int const NTAB = 32;
  long const NDIV = (1 + IMM1/NTAB);
  double const EPS = 1.2e-7;
  double const RNMX = 1.0 - EPS;

  int j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  double temp;
  
  if (*idum <= 0){
    if ( -(*idum) < 1) *idum = 1;
    else *idum = -(*idum);
    idum2 = (*idum);
    for ( j = NTAB + 7; j >= 0; j--){
      k = (*idum)/ IQ1;
      *idum = IA1 * (*idum - k*IQ1) - k * IR1;
      if ( *idum < 0) *idum += IM1;
      if ( j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }
  
  k = (*idum) / IQ1;
  *idum = IA1 * (*idum - k*IQ1) - k * IR1;
  if ( *idum < 0) *idum += IM1;
  k = idum2 / IQ2;
  idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
  if (idum2 < 0) idum2 += IM2;
  j = iy/ NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  
  if ( iy < 1) iy += IMM1;
  if (( temp = AM*iy) > RNMX){
    return RNMX;
  }
  else{
    return temp;
  }
}
