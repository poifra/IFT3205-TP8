/*------------------------------------------------------*/
/* Prog    : Tp10_IFT3205.c                             */
/* Auteur  : Francois Poitras et Charles Langlois       */
/* Date    : --/--/2010                                 */
/* version :                                            */ 
/* langage : C                                          */
/* labo    : DIRO                                       */
/*------------------------------------------------------*/


/*------------------------------------------------*/
/* FICHIERS INCLUS -------------------------------*/
/*------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FonctionDemo10.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/   
/*------------------------------------------------*/
#define NAME_VISUALISER_IMG "./display "
#define NAME_VISUALISER     "./ViewSig.sh "

/*------------------------------------------------*/
/* PROTOTYPE DE FONCTIONS  -----------------------*/   
/*------------------------------------------------*/

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL   -------------------------*/                     
/*------------------------------------------------*/
int main(int argc,char **argv)
{
    int i,j,n;
    char BufSystVisuSig[100];
    int length;

    float*  SignX=LoadSignalDat("SOUND_GoodMorningVietnam",&length);
    float*  SignY=fmatrix_allocate_1d(length);
    
    float*  b=fmatrix_allocate_1d(length);

    for(n=0;n<length;n++) {
        b[n] = (rand()/(float)RAND_MAX) * 2 - 1;
    }
    
    float SamplingRate=11025;

    float am[] = {1,3,4};

    float total = 0;
    int m;
    
    for(n=0;n<length;n++) {

        total = 0;

        m = 1; total += m * sin(2 * PI * am[m - 1] * 440 * n / SamplingRate + 7 * sqrt(m) * sin(2 * PI * 5 * n / SamplingRate));
        m = 2; total += m * sin(2 * PI * am[m - 1] * 440 * n / SamplingRate + 7 * sqrt(m) * sin(2 * PI * 5 * n / SamplingRate));
        m = 3; total += m * sin(2 * PI * am[m - 1] * 440 * n / SamplingRate + 7 * sqrt(m) * sin(2 * PI * 5 * n / SamplingRate));
        
        SignY[n] = total * SignX[n];
    }

    for(n=0;n<length;n++) {
        SignY[n] *= 0.5;
    }
    
    //Sauvegarde
    SaveSignalDatWav("SignalOut9",SignY,length,SamplingRate);
    //SaveSignalDat("SOUND_GoodMorningVietnam1",SignY,length);
 
    //Visu
    //strcpy(BufSystVisuSig,NAME_VISUALISER);
    //strcat(BufSystVisuSig,"SOUND_GoodMorningVietnam1.dat&");
    //printf(" %s",BufSystVisuSig);
    //system(BufSystVisuSig);
       
 
    //==End=========================================================

    //retour sans probleme
    printf("\n C'est fini ... \n\n");
    return 0; 	 
}


