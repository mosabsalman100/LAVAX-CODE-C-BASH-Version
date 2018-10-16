# LAVAX-CODE-C-BASH-Version

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define nSIA 1 
//#define NA 1536
#define NRUN 30000
#define TIME_STEP 0.002
#define ND 1 // number of diffsion path 1d, 2d or, 3d
#define THERMO 1 // at how many step show the results.
#define DUMP_CONFIG_FILE "/home2/mosab/2017.03.22_Uniaxialstress/3d_zr_Vacancy_and_Int/tets_iron/deformation/ackland_potential_model/Using_Our_own_codes/perfect_1"
#define LAMMPS_EXE "mpirun -np 12 /home/bin/lammps_voronoi"			// lammps execution command
#define INP_LAMMPS_REF       "/home2/mosab/2017.03.22_Uniaxialstress/3d_zr_Vacancy_and_Int/tets_iron/deformation/ackland_potential_model/Using_Our_own_codes/perfect_1/resource/input_diffusion_MSD"	// reference input file (to be read and then be modified)
#define NLINE_INP_LAMMPS_REF 35						// the number of lines in the reference input file
#define POT_FILE1 "/home2/mosab/2017.03.22_Uniaxialstress/3d_zr_Vacancy_and_Int/tets_iron/deformation/ackland_potential_model/Using_Our_own_codes/perfect_1/resource/Fe_ackland.eam.fs"		// potential file name
#define RES_FILE1 "/home2/mosab/2017.03.22_Uniaxialstress/3d_zr_Vacancy_and_Int/tets_iron/deformation/ackland_potential_model/Using_Our_own_codes/perfect_1/resource/data.after_defortmation_event_perfect_1sia"

//void func_create_MSD_files(char dump[],char dstar[], double x[][1], double y[][1],double z[][1],double late_s[3][2],double atom_confB[][5],double atom_confA[][5]);
void func_create_MSD_files(char dump[],char dstar[]);
void func_create_input_files(char fold1[], char proj1[], int Temp);


int main()
{	
	int N_TIMING,s,i,n,ns,nr,NA,Temp;
	NA=1537;
	char fname_dump_config[300];
	char fname1[300],RUN_FOLDER[300];
	char cline[300],command1[300],proj1[300];
	double atom_config[NA][8];
	double lattice_size[3][2];
	char fname_dump_star[300];
	char fname_MSD[300];
	double MSDx[n][1],MSDy[n][1],MSDz[n][1],MSD[n][1];
	double atom_configB[n][5],atom_configA[n][5];
	char dump[300],dstar[300]; 
	double x[n][1]; 
	double y[n][1];
	double z[n][1];
	double late_s[3][2];
	double atom_confB[n][5];
	double atom_confA[n][5];


	N_TIMING=NRUN/THERMO;
	FILE *fr,*fw;
	


                sprintf(command1,"mkdir MSD_config\n");
                system(command1);
		
                sprintf(command1,"mkdir dump_config\n");
                system(command1);
	        sprintf(command1,"mkdir MSD_AVG_config\n");
                system(command1);

	for(Temp=350;Temp<650;Temp=Temp+50)
	{
	printf("#Temp: %d\n",Temp);//pka_ene);
		sprintf(proj1,"Temp%d_%dSIA",Temp,nSIA);
		func_create_input_files(RUN_FOLDER, proj1, Temp);
		sprintf(command1,"date");
		system(command1);
		sprintf(command1,"%s < input.%s > output.%s\n",LAMMPS_EXE,proj1,proj1);
		system(command1);


	 for(s=0;s<N_TIMING+1;s++)
        {	
	        sprintf(fname_dump_config,"%s/config.dump.%d.%d",DUMP_CONFIG_FILE,Temp,nSIA);
	        if( (fr=fopen(fname_dump_config,"r"))==NULL )   {  printf("error in atom config file: %s\n",fname_dump_config);  exit(1);  }
		if(s>0){for(i=0;i<((NA+9)*s);i++) {fgets(cline,500,fr);}}
		fgets(cline,500,fr);
 	        fscanf(fr,"%d\n",&ns);
		//printf("N_timing=%d  s=%d",ns,s);
           	if((s*THERMO)!=ns)       {  printf("error: mismatch in ns number\n");  exit(1);  } 
		fgets(cline,500,fr);
 	        fscanf(fr,"%d\n",&NA);
		//printf("NA=%d",NA);
		fgets(cline,500,fr);

		for(nr=0;nr<3;nr++){
              				for(i=0;i<2;i++)        {fscanf(fr,"%lf",&lattice_size[nr][i]);}
					//printf("[nr][i]=%f  ",lattice_size[nr][i]);
					fgets(cline,500,fr);
					}
		fgets(cline,500,fr);
		//fgets(cline,500,fr);
				for(n=0;n<NA;n++){
              				for(i=0;i<8;i++)       { fscanf(fr,"%lf",&atom_config[n][i]);}//printf("[n][i]=%lf\n  ",atom_config[n][i]);}
					fgets(cline,500,fr);
						}
						fclose(fr);

		sprintf(fname1,"dump_config.%d.%dK.%dSIA",s*THERMO,Temp,nSIA);
       		if( (fw=fopen(fname1,"w"))==NULL )      {  printf("error: cannot open dump_config");  exit(1);  }
		fprintf(fw,"ITEM: TIMESTEP\n");
		fprintf(fw,"%d\n",ns);
		fprintf(fw,"ITEM: NUMBER OF ATOMS\n");
		fprintf(fw,"%d\n",NA);
		fprintf(fw,"ITEM: BOX BOUNDS pp pp pp\n");

		for(nr=0;nr<3;nr++){
              				for(i=0;i<2;i++)        {fprintf(fw,"%lf",lattice_size[nr][i]);}
					fprintf(fw,"\n");
					}

		fprintf(fw,"ITEM: ATOMS id type x y z ix iy iz \n");


     	  	 for(n=0;n<NA;n++)
     		   {
			for(i=0;i<8;i++) {fprintf(fw,"%lf    ",atom_config[n][i]);}
			fprintf(fw,"\n");
			}
	        fclose(fw);
	
			}
			}	
	//func_create_MSD_files(dump,dstar, x, y,z,late_s,atom_confB,atom_confA);
	func_create_MSD_files(dump,dstar);
	
	
        sprintf(command1,"mv %s/dump_config.* %s/dump_config",DUMP_CONFIG_FILE,DUMP_CONFIG_FILE);
        system(command1);
	sprintf(command1,"mv %s/MSDavg_config.* %s/MSD_AVG_config",DUMP_CONFIG_FILE,DUMP_CONFIG_FILE);
        system(command1);
	sprintf(command1,"mv %s/MSD_config.* %s/MSD_config",DUMP_CONFIG_FILE,DUMP_CONFIG_FILE);
        system(command1);
	return 0;
}

void func_create_input_files(char fold1[], char proj1[], int Temp)
{
	int i,j,k;
	char cline1[500],fname1[300];
	FILE *fw,*fr;

    	sprintf(fname1,"input.%s",proj1);
      	if( (fw=fopen(fname1,"w"))==NULL )      {  printf("error in open %s file\n",fname1);  exit(1);  }
     	fprintf(fw,"# Temp=%d    Number of SIA=%d\n",Temp,nSIA);

   	if( (fr=fopen(INP_LAMMPS_REF,"r"))==NULL )      {  printf("error in open INP_LAMMPS_REF file\n");  exit(1);  }

      	for(i=0;i<NLINE_INP_LAMMPS_REF;i++)
      	{
          	fgets(cline1,500,fr);
              	if(i==1)        fprintf(fw,"variable T equal       %d\n",Temp);
              	else if(i==2)   fprintf(fw,"variable s equal        %d\n",nSIA);
              	else if(i==3)   fprintf(fw,"variable pot_file1 string    %s\n",POT_FILE1);
              	else if(i==4)   fprintf(fw,"variable res_file1 string    %s\n",RES_FILE1);


             	else    fprintf(fw,"%s",cline1);
      	} 
     	fclose(fr);
     	fclose(fw);
}

void func_create_MSD_files(char dump[],char dstar[])

{
int N_TIMING,s,i,n,ns,nr,NA,j,Temp;
	NA=1537;
	char fname_dump_config[300];
	char fname1[300];
	char cline[300];
	double lattice_size[3][2];
	char fname_dump_star[300];
	char fname_MSD[300],fname_MSDavg[300];
	double MSD[NA],MSDavg[N_TIMING],MSD_AVG[N_TIMING],Diff_Coeff[N_TIMING];
	double atom_configB[NA][8],atom_configA[NA][8];
	double x[NA]; 
	double y[NA];
	double z[NA];
	double x2[NA]; 
	double y2[NA];
	double z2[NA];
	double late_s[3][2];
	double atom_confB[NA][8];
	double atom_confA[NA][8];
	FILE *fr,*fw,*fr1,*fw1,*fw2;

	N_TIMING=NRUN/THERMO;
	for(Temp=350;Temp<650;Temp=Temp+50)
	{	

	sprintf(fname_dump_config,"dump_config.0.%dK.%dSIA",Temp,nSIA);
        if( (fr=fopen(fname_dump_config,"r"))==NULL )   {  printf("error in opening atom config file: %s\n",fname_dump_config);  exit(1);}
                
         for(i=0;i<(9);i++) {fgets(cline,500,fr);}
       	 for(n=0;n<NA;n++){
//         for(i=0;i<(9);i++) {fgets(cline,500,fr);}
         for(i=0;i<8;i++)       { fscanf(fr,"%lf",&atom_configA[n][i]);}
                            	    fgets(cline,500,fr);}
                                 fclose(fr);

	for(s=1;s<N_TIMING;s++)
	{
		sprintf(fname_dump_star,"dump_config.%d.%dK.%dSIA",s*THERMO,Temp,nSIA);
                if( (fr1=fopen(fname_dump_star,"r"))==NULL )      {  printf("error: cannot open dump_config");  exit(1);  }
  
      for(j=0;j<(9);j++) {fgets(cline,500,fr1);}
	for(n=0;n<NA;n++){
				        for(j=0;j<8;j++)       { fscanf(fr1,"%lf",&atom_configB[n][j]);}//printf("[n][j]=%lf\n  ",atom_config[n][i]);}
					}
					fclose(fr1);	
		for(n=0;n<NA;n++){			
			
				x2[n]=pow(((atom_configB[n][2]+atom_configB[n][5])-(atom_configA[n][2]+atom_configA[n][5])),2);
				y2[n]=pow(((atom_configB[n][3]+atom_configB[n][6])-(atom_configA[n][3]+atom_configA[n][6])),2);
				z2[n]=pow(((atom_configB[n][4]+atom_configB[n][7])-(atom_configA[n][4]+atom_configA[n][7])),2);

				MSD[n]=sqrt(x2[n]+y2[n]+z2[n]);
				}				
			
		sprintf(fname_MSD,"MSD_config.%d.%dK.%dSIA",s*THERMO,Temp,nSIA);
       		if( (fw1=fopen(fname_MSD,"w"))==NULL )      {  printf("error: cannot open dump_config");  exit(1);  }
	

		fprintf(fw1,"id_original  ID_current   x2   y2   z2   MSD \n");
   	        for(n=0;n<NA;n++)
     		   {
			fprintf(fw1,"%lf  %lf   %lf   %lf   %lf   %lf   ",atom_configA[n][0],atom_configB[n][0],x2[n],y2[n],z2[n],MSD[n]);
			fprintf(fw1,"\n");
			MSDavg[s]+=MSD[n];
			}
	        fclose(fw1);
	}
		

    sprintf(fname_MSDavg,"MSDavg_config.%dK.%dSIA",Temp,nSIA);
                if( (fw2=fopen(fname_MSDavg,"w"))==NULL )      {  printf("error: cannot open p_config");  exit(1);  }


                fprintf(fw2,"Step    MSD_avg	Diffusion_Coefficient \n");
                for(s=0;s<N_TIMING;s++)
			{	
				MSD_AVG[s]=MSDavg[s]/NA;
				Diff_Coeff[s]=MSD_AVG[s]/(2*ND*TIME_STEP*(s*THERMO));
				fprintf(fw2,"%d     %lf	      %lf ",s*THERMO,MSD_AVG[s],Diff_Coeff[s]);
				fprintf(fw2,"\n");
			}
			fclose(fw2);
	
}
}
