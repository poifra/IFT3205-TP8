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

    //=================================================
    //Question 2.1
    //------------
    //
    // Sample Rate = 11025 ech/sec
    // Frequence Échant= 11025
    //
    //  
    //=================================================

    float*  SignX=LoadSignalDat("SOUND_GoodMorningVietnam",&length);
    float*  SignY=fmatrix_allocate_1d(length);
    
    float*  b=fmatrix_allocate_1d(length);

    for(n=0;n<length;n++) {
        b[n] = (rand()/(float)RAND_MAX) * 2 - 1;
    }
    
    float SamplingRate=11025;

    float R = 0.6;
    float f0 = 0.4;
    
    for(n=0;n<length;n++) {
        SignY[n] = SignX[n];

        if(n >= 1)
            SignY[n] += 2 * R * cos(2 * PI * n * f0) * SignY[n - 1];
        
        if(n >= 2)
            SignY[n] += -SignX[n - 2] + R*R * SignY[n - 2];
    }

    for(n=0;n<length;n++) {
        SignY[n] *= 0.5;
    }
    
    //Sauvegarde
    SaveSignalDatWav("SignalOut7",SignY,length,SamplingRate); 
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


