/*---------------------------------------------------*/
/* module  : FonctionDemo10.c                         */
/* auteur  : Mignotte Max                            */
/* date    : --/--/2010                              */              
/* langage : C                                       */
/* labo    : DIRO                                    */
/*---------------------------------------------------*/

/*------------------------------------------------*/
/* FICHIERS INCLUS -------------------------------*/
/*------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "FonctionDemo10.h"

/*------------------------------------------------*/
/* FONCTIONS -------------------------------------*/
/*------------------------------------------------*/
/*---------------------------------------------------------*/
/*  Alloue de la memoire pour une matrice 1d de float      */
/*---------------------------------------------------------*/
float* fmatrix_allocate_1d(int hsize)
 {
  float* matrix;

  matrix=(float*)malloc(sizeof(float)*hsize); 
  if (matrix==NULL) printf("probleme d'allocation memoire");

  return matrix; 
 }

/*----------------------------------------------------------*/
/*  Alloue de la memoire pour une matrice 2d de float       */
/*----------------------------------------------------------*/
float** fmatrix_allocate_2d(int vsize,int hsize)
 {
  int i;
  float** matrix;
  float *imptr;

  matrix=(float**)malloc(sizeof(float*)*vsize);
  if (matrix==NULL) printf("probleme d'allocation memoire");

  imptr=(float*)malloc(sizeof(float)*hsize*vsize);
  if (imptr==NULL) printf("probleme d'allocation memoire");
 
  for(i=0;i<vsize;i++,imptr+=hsize) matrix[i]=imptr;
  return matrix;
 }

/*----------------------------------------------------------*/
/* Libere la memoire de la matrice 1d de float              */
/*----------------------------------------------------------*/
void free_fmatrix_1d(float* pmat)
 { 
  free(pmat); 
 }

//----------------------------------------------------------*/
/* Libere la memoire de la matrice 2d de float              */
/*----------------------------------------------------------*/
void free_fmatrix_2d(float** pmat)
 { 
  free(pmat[0]);
  free(pmat);
 }

/*----------------------------*/
/* -LECTURE/SAUVEGARDE IMAGE- */
/*----------------------------*/
/*----------------------------------------------------------*/
/* Chargement de l'image de nom <name> (en pgm)             */
/*----------------------------------------------------------*/
float** LoadImagePgm(char* name,int *length,int *width)
 {
  int i,j,k;
  unsigned char var;
  char buff[NBCHAR];
  float** mat;
  char stringTmp1[NBCHAR],stringTmp2[NBCHAR];
  int ta1,ta2,ta3;
  FILE *fic;

  /*-----nom du fichier pgm-----*/
  strcpy(buff,name);
  strcat(buff,".pgm");
  printf("---> Ouverture de %s",buff);

  /*----ouverture du fichier----*/
  fic=fopen(buff,"r");
  if (fic==NULL)
    { printf("\n- Grave erreur a l'ouverture de %s  -\n",buff);
      exit(-1); }

  /*--recuperation de l'entete--*/
  fgets(stringTmp1,100,fic);
  fgets(stringTmp2,100,fic);
  fscanf(fic,"%d %d",&ta1,&ta2);
  fscanf(fic,"%d\n",&ta3);

  /*--affichage de l'entete--*/
  printf("\n\n--Entete--");
  printf("\n----------");
  printf("\n%s%s%d %d \n%d\n",stringTmp1,stringTmp2,ta1,ta2,ta3);

  *length=ta1;
  *width=ta2;
  mat=fmatrix_allocate_2d(*length,*width);
   
  /*--chargement dans la matrice--*/
     for(i=0;i<*length;i++)
      for(j=0;j<*width;j++)  
        { fread(&var,1,1,fic);
          mat[i][j]=var; }

   /*---fermeture du fichier---*/
  fclose(fic);

  return(mat);
 }


/*----------------------------------------------------------*/
/* Sauvegarde de l'image de nom <name> au format pgm        */
/*----------------------------------------------------------*/
void SaveImagePgm(char* name,float** mat,int length,int width)
 {
  int i,j,k;
  char buff[NBCHAR];
  FILE* fic;
  time_t tm;

  /*--extension--*/
  strcpy(buff,name);
  strcat(buff,".pgm");

  /*--ouverture fichier--*/
  fic=fopen(buff,"w");
    if (fic==NULL) 
        { printf(" Probleme dans la sauvegarde de %s",buff); 
          exit(-1); }
  printf("\n Sauvegarde de %s au format pgm\n",name);

  /*--sauvegarde de l'entete--*/
  fprintf(fic,"P5");
  if ((ctime(&tm))==NULL) fprintf(fic,"\n#\n");
  else fprintf(fic,"\n# IMG Module, %s",ctime(&tm));
  fprintf(fic,"%d %d",width,length);
  fprintf(fic,"\n255\n");

  /*--enregistrement--*/
     for(i=0;i<length;i++)
      for(j=0;j<width;j++) 
        fprintf(fic,"%c",(char)mat[i][j]);
   
  /*--fermeture fichier--*/
   fclose(fic); 
 } 

/*--------------*/
/* FOURIER FFT--*/
/*--------------*/
/*------------------------------------------------*/
/*  FOURN ----------------------------------------*/
/*------------------------------------------------*/
void fourn(float data[], unsigned long nn[], int ndim, int isign)
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	float tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for (ntot=1,idim=1;idim<=ndim;idim++)
		ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
						tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}


/*----------------------------------------------------------*/
/* FFTDD                                                    */
/*----------------------------------------------------------*/
void FFTDD(float** mtxR,float** mtxI,int lgth, int wdth)
{
 int i,j;
 int posx,posy;

 float* data;
 float* ImgFreqR;
 float* ImgFreqI;
 unsigned long* nn;

 /*allocation memoire*/
 data=(float*)malloc(sizeof(float)*(2*wdth*lgth)+1);
 ImgFreqR=(float*)malloc(sizeof(float)*(wdth*lgth));
 ImgFreqI=(float*)malloc(sizeof(float)*(wdth*lgth));
 nn=(unsigned long*)malloc(sizeof(unsigned long)*(FFT2D+1)); 

 /*Remplissage de nn*/
 nn[1]=lgth; nn[2]=wdth;

 /*Remplissage de data*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
   { data[2*(i*lgth+j)+1]=mtxR[i][j];
     data[2*(i*lgth+j)+2]=mtxI[i][j]; }

 /*FFTDD*/
 fourn(data,nn,FFT2D,FFT);

 /*Remplissage*/
 for(i=0;i<(wdth*lgth);i++)
  { ImgFreqR[i]=data[(2*i)+1];
    ImgFreqI[i]=data[(2*i)+2];  }

 /*Conversion en Matrice*/
 for(i=0;i<(wdth*lgth);i++)
  { posy=(int)(i/wdth);
    posx=(int)(i%wdth);

    mtxR[posy][posx]=ImgFreqR[i];
    mtxI[posy][posx]=ImgFreqI[i];}

 /*Liberation memoire*/
 free(data);
 free(ImgFreqR);
 free(ImgFreqI);
 free(nn);
}


/*----------------------------------------------------------*/
/* IFFTDD                                                   */
/*----------------------------------------------------------*/
void IFFTDD(float** mtxR,float**  mtxI,int lgth,int wdth)
{
 int i,j;
 int posx,posy;

 float* data;
 float* ImgFreqR;
 float* ImgFreqI;
 unsigned long* nn;

 /*allocation memoire*/
 data=(float*)malloc(sizeof(float)*(2*wdth*lgth)+1);
 ImgFreqR=(float*)malloc(sizeof(float)*(wdth*lgth));
 ImgFreqI=(float*)malloc(sizeof(float)*(wdth*lgth));
 nn=(unsigned long*)malloc(sizeof(unsigned long)*(FFT2D+1));

 /*Remplissage de nn*/
 nn[1]=lgth; nn[2]=wdth;

 /*Remplissage de data*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
   { data[2*(i*lgth+j)+1]=mtxR[i][j];
     data[2*(i*lgth+j)+2]=mtxI[i][j]; }

 /*FFTDD*/
 fourn(data,nn,FFT2D,IFFT);

 /*Remplissage*/
 for(i=0;i<(wdth*lgth);i++)
  { ImgFreqR[i]=data[(2*i)+1];
    ImgFreqI[i]=data[(2*i)+2]; }

 /*Conversion en Matrice*/
 for(i=0;i<(wdth*lgth);i++)
  { posy=(int)(i/wdth);
    posx=(int)(i%wdth);

   mtxR[posy][posx]=ImgFreqR[i]/(wdth*lgth);  
   mtxI[posy][posx]=ImgFreqI[i]/(wdth*lgth); }

 /*Liberation memoire*/
 free(data);
 free(ImgFreqR);
 free(ImgFreqI);
 free(nn);
}

/*----------------------------------------------------------*/
/* FFT1D                                                    */
/*----------------------------------------------------------*/
void FFT1D(float* vctR,float* vctI,int lgth)
{
 int i;
 float* data;
 unsigned long* nn;

 /*allocation memoire*/
 data=(float*)malloc(sizeof(float)*((2*lgth)+1));
 nn=(unsigned long*)malloc(sizeof(unsigned long)*2); 
 nn[1]=lgth; 

 /*Remplissage de data*/
 for(i=0;i<lgth;i++)  
   { data[(2*i)+1]=vctR[i];
     data[(2*i)+2]=vctI[i]; 
    }

 /*FFTDD*/
 fourn(data,nn,1,FFT);

 /*Resultat*/
 for(i=0;i<lgth;i++)
  { vctR[i]=data[(2*i)+1];
    vctI[i]=data[(2*i)+2];  }

 /*Liberation memoire*/
 free(data);
 free(nn);
}


/*----------------------------------------------------------*/
/* IFFT1D                                                   */
/*----------------------------------------------------------*/
void IFFT1D(float* vctR,float* vctI,int lgth)
{
 int i;
 float* data;
 unsigned long* nn;

 /*allocation memoire*/
 data=(float*)malloc(sizeof(float)*((2*lgth)+1));
 nn=(unsigned long*)malloc(sizeof(unsigned long)*2); 
 nn[1]=lgth; 

 /*Remplissage de data*/
 for(i=0;i<lgth;i++)  
   { data[(2*i)+1]=vctR[i];
     data[(2*i)+2]=vctI[i]; }

 /*FFTDD*/
 fourn(data,nn,1,IFFT);

 /*Resultat*/
 for(i=0;i<lgth;i++)
  { vctR[i]=data[(2*i)+1];
    vctI[i]=data[(2*i)+2];  }

 /*Liberation memoire*/
 free(data);
 free(nn);
}

/*--------------*/
/* FOURIER DFT--*/
/*--------------*/
//----------------------------------------------------------
// dft                                         
//---------------------------------------------------------- 
void dft(float* datar,float* datai,float* dataR,float* dataI,int sz,int inverse)
{
 int x,y;
 float a,ca,sa;
 float pi2=(inverse)?2.0*PI:-2.0*PI;
 float invs=1.0/(float)sz;

 for(y=0;y<sz;y++) 
    { dataR[y]=0;
      dataI[y]=0;

      for(x=0;x<sz;x++) 
        { a=pi2*y*x*invs;
          ca=cos(a);
          sa=sin(a);
          dataR[y]+=datar[x]*ca-datai[x]*sa;
          dataI[y]+=datar[x]*sa+datai[x]*ca; }

       if(!inverse) 
         { dataR[y]*=invs; 
           dataI[y]*=invs; }
    }
}
                
//----------------------------------------------------------
// DFT                                         
//---------------------------------------------------------- 
void DFTDD(float** datar,float** datai,int lgth,int wdth)
{
 int i,j;

 //Memory Allocation
 float* Vctr=fmatrix_allocate_1d(lgth);
 float* Vcti=fmatrix_allocate_1d(lgth);
 float* VctR=fmatrix_allocate_1d(lgth);
 float* VctI=fmatrix_allocate_1d(lgth);

 float* Vecr=fmatrix_allocate_1d(wdth);
 float* Veci=fmatrix_allocate_1d(wdth);
 float* VecR=fmatrix_allocate_1d(wdth);
 float* VecI=fmatrix_allocate_1d(wdth);

 float** ImgR=fmatrix_allocate_2d(lgth,wdth);
 float** ImgI=fmatrix_allocate_2d(lgth,wdth);
 float** dataR=fmatrix_allocate_2d(lgth,wdth);
 float** dataI=fmatrix_allocate_2d(lgth,wdth);
 

 //Iteration on Columns
 for(j=0;j<wdth;j++) 
    { 
      for(i=0;i<lgth;i++) { Vctr[i]=datar[i][j]; Vcti[i]=datai[i][j]; }
      dft(Vctr,Vcti,VctR,VctI,lgth,0); 
      for(i=0;i<lgth;i++) { ImgR[i][j]=VctR[i];  ImgI[i][j]=VctI[i]; }
    }

 //Iteration on Rows
 for(i=0;i<lgth;i++)
   {
    for(j=0;j<wdth;j++) { Vecr[j]=ImgR[i][j]; Veci[j]=ImgI[i][j]; }
    dft(Vecr,Veci,VecR,VecI,wdth,0); 
    for(j=0;j<wdth;j++) { dataR[i][j]=VecR[j];  dataI[i][j]=VecI[j]; }
   }

 //Transfert
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
   { datar[i][j]=dataR[i][j];
     datai[i][j]=dataI[i][j]; }

 //desallocation memoire
 free_fmatrix_1d(Vecr);
 free_fmatrix_1d(Veci);
 free_fmatrix_1d(VecR);
 free_fmatrix_1d(VecI);

 free_fmatrix_1d(Vctr);
 free_fmatrix_1d(Vcti);
 free_fmatrix_1d(VctR);
 free_fmatrix_1d(VctI);

 free_fmatrix_2d(ImgR);
 free_fmatrix_2d(ImgI); 
}

//----------------------------------------------------------
// DFT                                         
//---------------------------------------------------------- 
void IDFTDD(float** datar,float** datai,int lgth,int wdth)
{
 int i,j;

 //Memory Allocation
 float* Vctr=fmatrix_allocate_1d(lgth);
 float* Vcti=fmatrix_allocate_1d(lgth);
 float* VctR=fmatrix_allocate_1d(lgth);
 float* VctI=fmatrix_allocate_1d(lgth);

 float* Vecr=fmatrix_allocate_1d(wdth);
 float* Veci=fmatrix_allocate_1d(wdth);
 float* VecR=fmatrix_allocate_1d(wdth);
 float* VecI=fmatrix_allocate_1d(wdth);

 float** ImgR=fmatrix_allocate_2d(lgth,wdth);
 float** ImgI=fmatrix_allocate_2d(lgth,wdth);
 float** dataR=fmatrix_allocate_2d(lgth,wdth);
 float** dataI=fmatrix_allocate_2d(lgth,wdth);
 

 //Iteration on Columns
 for(j=0;j<wdth;j++) 
    { 
      for(i=0;i<lgth;i++) { Vctr[i]=datar[i][j]; Vcti[i]=datai[i][j]; }
      dft(Vctr,Vcti,VctR,VctI,lgth,1); 
      for(i=0;i<lgth;i++) { ImgR[i][j]=VctR[i];  ImgI[i][j]=VctI[i]; }
    }

 //Iteration on Rows
 for(i=0;i<lgth;i++)
   {
    for(j=0;j<wdth;j++) { Vecr[j]=ImgR[i][j]; Veci[j]=ImgI[i][j]; }
    dft(Vecr,Veci,VecR,VecI,wdth,1); 
    for(j=0;j<wdth;j++) { dataR[i][j]=VecR[j];  dataI[i][j]=VecI[j]; }
   }

 //Transfert
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
   { datar[i][j]=dataR[i][j];
     datai[i][j]=dataI[i][j]; }

 //desallocation memoire
 free_fmatrix_1d(Vecr);
 free_fmatrix_1d(Veci);
 free_fmatrix_1d(VecR);
 free_fmatrix_1d(VecI);

 free_fmatrix_1d(Vctr);
 free_fmatrix_1d(Vcti);
 free_fmatrix_1d(VctR);
 free_fmatrix_1d(VctI);

 free_fmatrix_2d(ImgR);
 free_fmatrix_2d(ImgI); 
}

//--------------------------//
//--- OUTILS SPECTRAL 2D ---//
//--------------------------//
/*----------------------------------------------------------*/
/* Mod2                                                     */
/*----------------------------------------------------------*/
void Mod(float** matM,float** matR,float** matI,int lgth,int wdth)
{
 int i,j;

 /*Calcul du module*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
 matM[i][j]=sqrt((matR[i][j]*matR[i][j])+(matI[i][j]*matI[i][j]));
}

/*----------------------------------------------------------*/
/* Mult                                                     */
/*----------------------------------------------------------*/
void Mult(float** mat,float coef,int lgth,int wdth)
{
 int i,j;

 /*multiplication*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
    { mat[i][j]*=coef;
      if (mat[i][j]>GREY_LEVEL) mat[i][j]=GREY_LEVEL; }
}

/*----------------------------------------------------------*/
/* Recal                                                    */
/*----------------------------------------------------------*/
void Recal(float** mat,int lgth,int wdth)
{
 int i,j;
 float max,min;

 /*Initialisation*/
 min=mat[0][0];

 /*Recherche du min*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    if (mat[i][j]<min) min=mat[i][j];

 /*plus min*/
   for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    mat[i][j]-=min;

   max=mat[0][0];
 /*Recherche du max*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
    if (mat[i][j]>max) max=mat[i][j];

 /*Recalibre la matrice*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   mat[i][j]*=(GREY_LEVEL/max);      
}

/*----------------------------------------------------------*/
/* Mult 2 matrices complexes                                */
/*----------------------------------------------------------*/
void MultMatrix(float** matRout,float** matIout,float** mat1Rin,float** mat1Iin,
float** mat2Rin,float** mat2Iin,int lgth,int wdth)
{
 int i,j;

 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   { matRout[i][j]=mat1Rin[i][j]*mat2Rin[i][j]-mat1Iin[i][j]*mat2Iin[i][j];
     matIout[i][j]=mat1Rin[i][j]*mat2Iin[i][j]+mat2Rin[i][j]*mat1Iin[i][j]; }
}

/*----------------------------------------------------------*/
/* Matrice complexe au carre                                */
/*----------------------------------------------------------*/
void SquareMatrix(float** matRout,float** matIout,float** matRin,float** matIin,int lgth,int wdth)
{
 int i,j;

 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   { matRout[i][j]=SQUARE(matRin[i][j])-SQUARE(matIin[i][j]);
     matIout[i][j]=2*matRin[i][j]*matIin[i][j]; }
}

/*----------------------------------------------------------*/
/*  CenterImg_                                              */
/*----------------------------------------------------------*/
void CenterImg_(float** mat,int length,int width)
{
  int i, j;
  float** tmp;
  
  tmp=fmatrix_allocate_2d(length,width);
  
  for(i=0;i<length;i++) for(j=0;j<width;j++) 
    tmp[i][j]=0.0;

  IFFTDD(mat,tmp,length,width);

  for(i=0;i<length;i++) for(j=0;j<width;j++) 
    {
    mat[i][j]*=pow(-1,i+j);
    tmp[i][j]*=pow(-1,i+j);
    }

  FFTDD(mat,tmp,length,width);

  free_fmatrix_2d(tmp);
}

/*----------------------------------------------------------*/
/*  CenterImg                                               */
/*----------------------------------------------------------*/
void CenterImg(float** mat,int lgth,int wdth)
{
 int i,j;
 int ci,cj;
 float** mattmp;

 /*Initialisation*/
 ci=(int)(lgth/2);
 cj=(int)(wdth/2);

 /*Allocation memoire*/
 mattmp=fmatrix_allocate_2d(lgth,wdth);

 /*Recadrage*/
 for(i=0;i<ci;i++) for(j=0;j<cj;j++)
 mattmp[ci+i][cj+j]=mat[i][j];

 for(i=ci;i<lgth;i++) for(j=cj;j<wdth;j++)
 mattmp[i-ci][j-cj]=mat[i][j];

 for(i=0;i<ci;i++) for(j=cj;j<wdth;j++)
 mattmp[ci+i][j-cj]=mat[i][j];

 for(i=ci;i<lgth;i++) for(j=0;j<cj;j++)
 mattmp[i-ci][cj+j]=mat[i][j];

 /*Transfert*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
  mat[i][j]=mattmp[i][j];

 /*desallocation memoire*/
 free_fmatrix_2d(mattmp);
}

//--------------------------//
//--- OUTILS SPECTRAL 1D ---//
//--------------------------//
/*----------------------------------------------------------*/
/* ModVct                                                     */
/*----------------------------------------------------------*/
void ModVct(float* vctM,float* vctR,float* vctI,int lgth)
{
 int i,j;

 /*Calcul du module*/
 for(i=0;i<lgth;i++)
 vctM[i]=sqrt((vctR[i]*vctR[i])+(vctI[i]*vctI[i]));
}

/*----------------------------------------------------------*/
/*  CenterVct                                               */
/*----------------------------------------------------------*/
void CenterVct(float* vct,int lgth)
{
 int i;
 int ci;
 float* vcttmp;

 /*Initialisation*/
 ci=(int)(lgth/2);

 /*Allocation memoire*/
 vcttmp=fmatrix_allocate_1d(lgth);

 /*Recadrage*/
 for(i=0;i<ci;i++) vcttmp[ci+i]=vct[i];
 for(i=ci;i<lgth;i++) vcttmp[i-ci]=vct[i];

 /*Transfert*/
 for(i=0;i<lgth;i++) vct[i]=vcttmp[i];

 /*desallocation memoire*/
 free_fmatrix_1d(vcttmp);
}

//-------------------//
//--- DEGRADATION ---//
//-------------------//
//----------------------------------------------------------
//  Gaussian noise  
//----------------------------------------------------------
float gaussian_noise(float var)
{
 float noise,theta;

 //Noise generation 
 noise=sqrt(-2*var*log(1.0-((float)rand()/RAND_MAX)));
 theta=(float)rand()*1.9175345E-4-PI;
 noise=noise*cos(theta);
 return noise;
}

//-----------------------------//
// -LECTURE/SAUVEGARDE SIGNAL -//
//-----------------------------//
/*----------------------------------------------------------*/
/* Chargement du signal de nom <name>.dat                   */
/*----------------------------------------------------------*/
float* LoadSignalDat(char* name,int* length)
 {
  int i;
  float ech1,ech2;
  char buff[NBCHAR];
  float Tech;
  FILE *fic;

  //Allocation
  float** Vct2D=fmatrix_allocate_2d(2,10000000);

  //nom du fichier dat
  strcpy(buff,name);
  strcat(buff,".dat");
  printf("\n  >> Ouverture de %s",buff);

  //ouverture du fichier
  fic=fopen(buff,"r");
  if (fic==NULL)
    { printf("\n- Grave erreur a l'ouverture de %s  -\n",buff);
      exit(-1); }
 
  //Lecture Donn�e & Longueur & Periode Ech
  for(i=0;;i++)
       { fscanf(fic,"%f %f",&ech1,&ech2);
         if (feof(fic)) break;
         //printf("\n[%f::%f]",ech1,ech2);
         Vct2D[0][i]=ech1;
         Vct2D[1][i]=ech2; }

  (*length)=i;
  Tech=Vct2D[0][1];
  Tech=1.0/Tech;
  Tech=(int)Tech;
  printf(" (%d echantillons)",(*length));
  printf("\n  >> Techantillonnage:: %.0f echantillons/seconde",Tech);

  //Chargement
  float* VctFinal=fmatrix_allocate_1d((*length));
  for(i=0;i<(*length);i++) VctFinal[i]=Vct2D[1][i];
 
  //Debug
  //for(i=0;i<(*length);i++) printf("\n [%d:%f]",i,VctFinal[i]);
  //for(i=0;i<(*length);i++) printf("\n [%f:%f]",Vct2D[0][i],Vct2D[1][i]);
 
   //End 
  fclose(fic);
  (*length)=i;
  free_fmatrix_2d(Vct2D);
  return VctFinal;
 }

/*----------------------------------------------------------*/
/* Sauvegarde du signal de nom <name>.dat                   */
/*----------------------------------------------------------*/
void SaveSignalDat(char* name,float* vct,int length)
 {
  int i;
  char buff[NBCHAR];
  FILE* fic;

  /*--extension--*/
  strcpy(buff,name);
  strcat(buff,".dat");

  /*--ouverture fichier--*/
  fic=fopen(buff,"w");
    if (fic==NULL) 
        { printf(" Probleme dans la sauvegarde de %s",buff); 
          exit(-1); }
  printf("\n Sauvegarde de %s au format dat\n",name);

  /*--enregistrement--*/
  for(i=0;i<length;i++) fprintf(fic,"%f %f\n",(float)i,vct[i]);

   /*--fermeture fichier--*/
   fclose(fic);  
 }

/*----------------------------------------------------------*/
/* Sauvegarde du signal de nom <name>.dat                   */
/*----------------------------------------------------------*/
void SaveSignalDat2(char* name,float* vct,int length,float Tech)
 {
  int i;
  char buff[NBCHAR];
  FILE* fic;

  /*--extension--*/
  strcpy(buff,name);
  strcat(buff,".dat");

  /*--ouverture fichier--*/
  fic=fopen(buff,"w");
    if (fic==NULL) 
        { printf(" Probleme dans la sauvegarde de %s",buff); 
          exit(-1); }
  printf("\n Sauvegarde de %s au format dat\n",name);

  /*--enregistrement--*/
  for(i=0;i<length;i++) fprintf(fic,"%f %f\n",(float)i*Tech,vct[i]);

   /*--fermeture fichier--*/
   fclose(fic);  
 }

/*----------------------------------------------------------*/
/* Sauvegarde du signal de nom <name>.dat                   */
/*----------------------------------------------------------*/
void SaveSignalDatWav(char* name,float* vct,int length,int SampRate)
 {
  int i;
  char buffName[NBCHAR];
  char buffSyst[NBCHAR];
  FILE* fic;

  /*--extension--*/
  strcpy(buffName,name);

  /*--ouverture fichier--*/
  fic=fopen("tmp.dat","w");
    if (fic==NULL) 
        { printf(" Probleme dans la sauvegarde de %s",buffName); 
          exit(-1); }
  printf("\n Sauvegarde de %s au format wav",name);

  /*--enregistrement--*/
  fprintf(fic,"; Sample Rate %d",SampRate);
  fprintf(fic,"; Channels 1\n");
  for(i=0;i<length;i++) fprintf(fic,"%f %f\n",(float)i,vct[i]);

   /*--fermeture fichier--*/
   fclose(fic);  

   /*--Conversion Dat -> Wav--*/
   sprintf(buffSyst,"sox tmp.dat %s.wav&",buffName);
   system(buffSyst);
   printf("\n Appel Systeme::sox tmp.dat %s.wav\n\n",buffName); 
 }
