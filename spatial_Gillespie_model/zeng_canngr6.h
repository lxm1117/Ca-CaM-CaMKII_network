struct voxel{
  double diffusion_lamda;
  double chemical_lamda;

  double Ca_pump; //the Ca2+ pump 
  double pump_area; //the surface pump area 
  double Ca_influx; //the Ca2+ influx molecular number by NMDA

  /* D[0] is diffusion coefficient to voxel [i-1,j,k]
     D[1] is diffusion coefficient to voxel [i,j-1,k]
     D[2] is diffusion coefficient to voxel [i+1,j,k]
     D[3] is diffusion coefficient to voxel [i,j+1,k]
     D[4] is diffusion coefficient to voxel [i,j,k-1]
     D[5] is diffusion coefficient to voxel [i,j,k+1]
  */
  
  double D_CaMKII[6];
  double D_Ca[6]; //diffusion coefficient of Ca2+ in six direction
  double D_CaM00[6];
  double D_CaM01[6];
  double D_CaM10[6];
  double D_CaM11[6];
  double D_CaM02[6];
  double D_CaM20[6];
  double D_CaM12[6];
  double D_CaM21[6];
  double D_CaM22[6];
  double D_CB00[6];
  double D_CB10[6];
  double D_CB01[6];
  double D_CB11[6];
  double D_CB20[6];
  double D_CB02[6];
  double D_CB21[6];
  double D_CB12[6];
  double D_CB22[6];

  //concentration of different molecules
  int CaMKII_num;
  int Ca;
  int CaM00;
  int CaM01;
  int CaM10;
  int CaM11;
  int CaM02;
  int CaM20;
  int CaM12;
  int CaM21;
  int CaM22;  

  int CaN;
  int CaNCaM00;
  int CaNCaM01;
  int CaNCaM10;
  int CaNCaM11;
  int CaNCaM02;
  int CaNCaM20;
  int CaNCaM21;
  int CaNCaM12;
  int CaNCaM22;

  int Ng;
  int NgCaM00;
  int NgCaM01;
  int NgCaM10;
  int NgCaM11;
  int NgCaM02;
  int NgCaM20;
  int NgCaM21;
  int NgCaM12;
  int NgCaM22;
  
  int CB00;
  int CB01;
  int CB10;
  int CB11;
  int CB02;
  int CB20;
  int CB12;
  int CB21;
  int CB22;

  //transition rate
  double CaM00_to_CaM10;
  double CaM10_to_CaM00;
  double CaM00_to_CaM01;
  double CaM01_to_CaM00;
  double CaM10_to_CaM11;
  double CaM11_to_CaM10;
  double CaM10_to_CaM20;
  double CaM20_to_CaM10;
  double CaM01_to_CaM11;
  double CaM11_to_CaM01;
  double CaM01_to_CaM02;
  double CaM02_to_CaM01;
  double CaM11_to_CaM21;
  double CaM21_to_CaM11;
  double CaM11_to_CaM12;
  double CaM12_to_CaM11;
  double CaM20_to_CaM21;
  double CaM21_to_CaM20;
  double CaM02_to_CaM12;
  double CaM12_to_CaM02;
  double CaM21_to_CaM22;
  double CaM22_to_CaM21;
  double CaM12_to_CaM22;
  double CaM22_to_CaM12;

  double CaM00_to_CaMKIICaM00;
  double CaMKIICaM00_to_CaM00;
  double CaM10_to_CaMKIICaM10;
  double CaMKIICaM10_to_CaM10;
  double CaM01_to_CaMKIICaM01;
  double CaMKIICaM01_to_CaM01;
  double CaM11_to_CaMKIICaM11;
  double CaMKIICaM11_to_CaM11;
  double CaM20_to_CaMKIICaM20;
  double CaMKIICaM20_to_CaM20;
  double CaM02_to_CaMKIICaM02;
  double CaMKIICaM02_to_CaM02;
  double CaM12_to_CaMKIICaM12;
  double CaMKIICaM12_to_CaM12;
  double CaM21_to_CaMKIICaM21;
  double CaMKIICaM21_to_CaM21;
  double CaM22_to_CaMKIICaM22;
  double CaMKIICaM22_to_CaM22;

  double CaMKIICaM00_to_CaMKIICaM10;
  double CaMKIICaM10_to_CaMKIICaM00;
  double CaMKIICaM00_to_CaMKIICaM01;
  double CaMKIICaM01_to_CaMKIICaM00;
  double CaMKIICaM10_to_CaMKIICaM11;
  double CaMKIICaM11_to_CaMKIICaM10;
  double CaMKIICaM10_to_CaMKIICaM20;
  double CaMKIICaM20_to_CaMKIICaM10;
  double CaMKIICaM01_to_CaMKIICaM11;
  double CaMKIICaM11_to_CaMKIICaM01;
  double CaMKIICaM01_to_CaMKIICaM02;
  double CaMKIICaM02_to_CaMKIICaM01;
  double CaMKIICaM11_to_CaMKIICaM21;
  double CaMKIICaM21_to_CaMKIICaM11;
  double CaMKIICaM11_to_CaMKIICaM12;
  double CaMKIICaM12_to_CaMKIICaM11;
  double CaMKIICaM20_to_CaMKIICaM21;
  double CaMKIICaM21_to_CaMKIICaM20;
  double CaMKIICaM02_to_CaMKIICaM12;
  double CaMKIICaM12_to_CaMKIICaM02;
  double CaMKIICaM21_to_CaMKIICaM22;
  double CaMKIICaM22_to_CaMKIICaM21;
  double CaMKIICaM12_to_CaMKIICaM22;
  double CaMKIICaM22_to_CaMKIICaM12;

  double CaM00_to_CaNCaM00;
  double CaNCaM00_to_CaM00;
  double CaM10_to_CaNCaM10;
  double CaNCaM10_to_CaM10;
  double CaM01_to_CaNCaM01;
  double CaNCaM01_to_CaM01;
  double CaM11_to_CaNCaM11;
  double CaNCaM11_to_CaM11;
  double CaM20_to_CaNCaM20;
  double CaNCaM20_to_CaM20;
  double CaM02_to_CaNCaM02;
  double CaNCaM02_to_CaM02;
  double CaM12_to_CaNCaM12;
  double CaNCaM12_to_CaM12;
  double CaM21_to_CaNCaM21;
  double CaNCaM21_to_CaM21;
  double CaM22_to_CaNCaM22;
  double CaNCaM22_to_CaM22;

  double CaNCaM00_to_CaNCaM10;
  double CaNCaM10_to_CaNCaM00;
  double CaNCaM00_to_CaNCaM01;
  double CaNCaM01_to_CaNCaM00;
  double CaNCaM10_to_CaNCaM11;
  double CaNCaM11_to_CaNCaM10;
  double CaNCaM10_to_CaNCaM20;
  double CaNCaM20_to_CaNCaM10;
  double CaNCaM01_to_CaNCaM11;
  double CaNCaM11_to_CaNCaM01;
  double CaNCaM01_to_CaNCaM02;
  double CaNCaM02_to_CaNCaM01;
  double CaNCaM11_to_CaNCaM21;
  double CaNCaM21_to_CaNCaM11;
  double CaNCaM11_to_CaNCaM12;
  double CaNCaM12_to_CaNCaM11;
  double CaNCaM20_to_CaNCaM21;
  double CaNCaM21_to_CaNCaM20;
  double CaNCaM02_to_CaNCaM12;
  double CaNCaM12_to_CaNCaM02;
  double CaNCaM21_to_CaNCaM22;
  double CaNCaM22_to_CaNCaM21;
  double CaNCaM12_to_CaNCaM22;
  double CaNCaM22_to_CaNCaM12;
  
  double CaM00_to_NgCaM00;
  double NgCaM00_to_CaM00;
  double CaM10_to_NgCaM10;
  double NgCaM10_to_CaM10;
  double CaM01_to_NgCaM01;
  double NgCaM01_to_CaM01;
  double CaM11_to_NgCaM11;
  double NgCaM11_to_CaM11;
  double CaM20_to_NgCaM20;
  double NgCaM20_to_CaM20;
  double CaM02_to_NgCaM02;
  double NgCaM02_to_CaM02;
  double CaM12_to_NgCaM12;
  double NgCaM12_to_CaM12;
  double CaM21_to_NgCaM21;
  double NgCaM21_to_CaM21;
  double CaM22_to_NgCaM22;
  double NgCaM22_to_CaM22;

  double NgCaM00_to_NgCaM10;
  double NgCaM10_to_NgCaM00;
  double NgCaM00_to_NgCaM01;
  double NgCaM01_to_NgCaM00;
  double NgCaM10_to_NgCaM11;
  double NgCaM11_to_NgCaM10;
  double NgCaM10_to_NgCaM20;
  double NgCaM20_to_NgCaM10;
  double NgCaM01_to_NgCaM11;
  double NgCaM11_to_NgCaM01;
  double NgCaM01_to_NgCaM02;
  double NgCaM02_to_NgCaM01;
  double NgCaM11_to_NgCaM21;
  double NgCaM21_to_NgCaM11;
  double NgCaM11_to_NgCaM12;
  double NgCaM12_to_NgCaM11;
  double NgCaM20_to_NgCaM21;
  double NgCaM21_to_NgCaM20;
  double NgCaM02_to_NgCaM12;
  double NgCaM12_to_NgCaM02;
  double NgCaM21_to_NgCaM22;
  double NgCaM22_to_NgCaM21;
  double NgCaM12_to_NgCaM22;
  double NgCaM22_to_NgCaM12;

  double CaMKIICaM00_to_Trapped00;
  double Trapped00_to_CaMKIICaM00;
  double CaMKIICaM10_to_Trapped10;
  double Trapped10_to_CaMKIICaM10;
  double CaMKIICaM01_to_Trapped01;
  double Trapped01_to_CaMKIICaM01;
  double CaMKIICaM11_to_Trapped11;
  double Trapped11_to_CaMKIICaM11;
  double CaMKIICaM20_to_Trapped20;
  double Trapped20_to_CaMKIICaM20;
  double CaMKIICaM02_to_Trapped02;
  double Trapped02_to_CaMKIICaM02;
  double CaMKIICaM12_to_Trapped12;
  double Trapped12_to_CaMKIICaM12;
  double CaMKIICaM21_to_Trapped21;
  double Trapped21_to_CaMKIICaM21;
  double CaMKIICaM22_to_Trapped22;
  double Trapped22_to_CaMKIICaM22;
  
  double Trapped00_to_Auton;
  double Auton_to_Trapped00;
  double Trapped10_to_Auton;
  double Auton_to_Trapped10;
  double Trapped01_to_Auton;
  double Auton_to_Trapped01;
  double Trapped11_to_Auton;
  double Auton_to_Trapped11;
  double Trapped20_to_Auton;
  double Auton_to_Trapped20;
  double Trapped02_to_Auton;
  double Auton_to_Trapped02;
  double Trapped12_to_Auton;
  double Auton_to_Trapped12;
  double Trapped21_to_Auton;
  double Auton_to_Trapped21;
  double Trapped22_to_Auton;
  double Auton_to_Trapped22;

  double Trapped00_to_Trapped10;
  double Trapped10_to_Trapped00;
  double Trapped00_to_Trapped01;
  double Trapped01_to_Trapped00;
  double Trapped10_to_Trapped11;
  double Trapped11_to_Trapped10;
  double Trapped10_to_Trapped20;
  double Trapped20_to_Trapped10;
  double Trapped01_to_Trapped11;
  double Trapped11_to_Trapped01;
  double Trapped01_to_Trapped02;
  double Trapped02_to_Trapped01;
  double Trapped11_to_Trapped21;
  double Trapped21_to_Trapped11;
  double Trapped11_to_Trapped12;
  double Trapped12_to_Trapped11;
  double Trapped20_to_Trapped21;
  double Trapped21_to_Trapped20;
  double Trapped02_to_Trapped12;
  double Trapped12_to_Trapped02;
  double Trapped21_to_Trapped22;
  double Trapped22_to_Trapped21;
  double Trapped12_to_Trapped22;
  double Trapped22_to_Trapped12;

  double Auton_to_CaMKII;
  double Auton_to_Capped;
  double Capped_to_Auton;

  //transition rate related to CB
  double CB00_to_CB10;
  double CB10_to_CB00;
  double CB00_to_CB01;
  double CB01_to_CB00;
  double CB10_to_CB11;
  double CB11_to_CB10;
  double CB10_to_CB20;
  double CB20_to_CB10;
  double CB01_to_CB11;
  double CB11_to_CB01;
  double CB01_to_CB02;
  double CB02_to_CB01;
  double CB11_to_CB21;
  double CB21_to_CB11;
  double CB11_to_CB12;
  double CB12_to_CB11;
  double CB20_to_CB21;
  double CB21_to_CB20;
  double CB02_to_CB12;
  double CB12_to_CB02;
  double CB12_to_CB22;
  double CB22_to_CB12;
  double CB21_to_CB22;
  double CB22_to_CB21;

  //point of linked list to store CaMKII
  
  struct CaMKII *head;
};

struct CaMKII{
  int index;
  int subunit[1][6]; /* each CaMKII has 12 subunits. Each subunit has 21 states. [0] is for CaMKII, [1] is for CaMKIICaM00;
			[2] is for CaMKIICaM10, [3] is for CaMKIICaM01, [4] is for CaMKIICaM11, [5] is for CaMKIICaM20,
			[6] is for CaMKIICaM02, [7] is for CaMKIICaM21, [8] is for CaMKIICaM12, [9] is for CaMKIICaM22, [10] is for 
			Trapped00, [11] is for Trapped10, [12] is for Trapped01, [13] is for Trapped11, [14] is for Trapped20, 
			[15] is for Trapped02, [16] is for Trapped21, [17] is for Trapped12, [18] is for Trapped22, [19] is for Autonomous,
			[20] is for Capped.
		     */

  struct CaMKII *next;
};
