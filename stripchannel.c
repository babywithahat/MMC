#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sfmisc.h>

#include "fm.h"
#include "env.h"
#include "mix.h"


int filecount(char* file_name)
{
	int lines = 0;
	char ch;

	FILE *fp;

	fp = fopen(file_name,"r");

	if(fp == NULL){
		perror("Error while opening the file.\n");
		printf("\n\tError while opening file: %s\n", file_name);
		exit(EXIT_FAILURE);
	}

	while((ch = fgetc(fp)) != EOF){
		if(ch == '\n'){
			lines++;
		}
	}
	
	fclose(fp);	
	return lines;
}

void file2array(char* file_name, double buffer[])
{
	double val;
	
	FILE *fp;

	fp = fopen(file_name,"r"); //read mode

	if(fp == NULL)
	{
		perror("Error while opening the file.\n");
		printf("\n\tError while opening file: %s\n", file_name);
		exit(EXIT_FAILURE);
	}
	
	int i = 0;
	while(i < filecount(file_name)){
		fscanf(fp, "%lf", &val);
		buffer[i] = val;
		++i;
	}
	
	fclose(fp);
	return ;
}


double durOfTone(int filenum,double startTime){
	
	int i;
	double toneLength = 0;
	double toneEnd = 0;
	
	char file_name[64] ;
	sprintf(file_name, "%dD", filenum) ;

	printf("\tLine counting file %s...", file_name);
	
	int mallocD;
	mallocD = filecount(file_name);
	printf("Done.\n");
	
	double* D;
	D = (double *)Malloc(mallocD * sizeof(double));
	
	printf("\tfile2array file %s...", file_name);
	
	file2array(file_name, D);
	
	printf("Done.\n");
	
	for(i =0;i<mallocD;i++)
		toneLength = toneLength + D[i];	

	Free(D);
	
	toneEnd = startTime + toneLength;
	
	return toneEnd;
}


void allocSetup(int filenum,
		int *mallocCar,
		int *mallocMod,
		int *mallocPD,
		int *mallocD,
		int *mallocAmp,
		int *mallocT,
		int *mallocA,
		int *mallocPan)
{
	char file_name[64];

	sprintf(file_name, "%dCarrier", filenum);
	printf("\tLine counting file %s...", file_name);
	*mallocCar = filecount(file_name);
	printf("%d\n",*mallocCar);
	
	sprintf(file_name, "%dModulator", filenum);
	printf("\tLine counting file %s...", file_name);
	*mallocMod = filecount(file_name);
	printf("%d\n",*mallocMod);

	sprintf(file_name, "%dPD", filenum);
	printf("\tLine counting file %s...", file_name);
	*mallocPD = filecount(file_name);
	printf("%d\n",*mallocPD);

	sprintf(file_name, "%dD", filenum);
	printf("\tLine counting file %s...", file_name);
	*mallocD = filecount(file_name);
	printf("%d\n",*mallocD);
	
	sprintf(file_name, "%dAmp", filenum);
	printf("\tLine counting file %s...", file_name);
	*mallocAmp = filecount(file_name);
	printf("%d\n",*mallocAmp);
	
	sprintf(file_name, "%dD", filenum);
	double* D;
	D = (double *)Malloc(*mallocD * sizeof(double));
	file2array(file_name, D);

	sprintf(file_name, "%dT", filenum);
	printf("\tLine counting file %s...", file_name);
	*mallocT = filecount(file_name);
	printf("%d\n",*mallocT);

	sprintf(file_name, "%dA", filenum);
	printf("\tLine counting file %s...", file_name);
	*mallocA = filecount(file_name);
	printf("%d\n",*mallocA);


	int i;
	double toneLength = 1;
	for(i = 0;i<*mallocD;i++) 
		toneLength = toneLength + D[i];

	*mallocPan = toneLength * 44100.0;
	printf("\tD[0] = %f\n", D[0]);
	printf("\ttoneLength = %f\n", toneLength);
	printf("\tmallocPan = %d\n", *mallocPan);
	
	Free(D);
}


void loadarray(int filenum,
		double Carrier[],
		double Modulator[],
		double PD[],
		double D[],
		double Amp[],
		int *nFMTones,
		double T[],
		double A[],
		int *nPoints,
		double xpan[],
		double ypan[]){
	
	char file_name[64];
	
	sprintf(file_name, "%dCarrier", filenum);
	printf("\tfile2array file %s...\n", file_name);
	file2array(file_name, Carrier);
	
	sprintf(file_name, "%dModulator", filenum);
	printf("\tfile2array file %s...\n", file_name);
	file2array(file_name, Modulator);
	
	sprintf(file_name, "%dPD", filenum);
	printf("\tfile2array file %s...\n", file_name);
	file2array(file_name, PD);

	sprintf(file_name, "%dD", filenum);
	printf("\tfile2array file %s...\n", file_name);
	file2array(file_name, D);

	sprintf(file_name, "%dAmp", filenum);
	printf("\tfile2array file %s...\n", file_name);
	file2array(file_name, Amp);
	
	*nFMTones = filecount(file_name);
	
	int i;
	double toneLength = 0;
	printf("test\n");
	for(i =0;i<*nFMTones;i++)
		toneLength = toneLength + D[i];
	
	printf("test\n");
	//changes Modulator Multiplier to Molulating Frequency
	for(i=0;i<*nFMTones;i++)
		Modulator[i] = Modulator[i] * Carrier[i];

	sprintf(file_name, "%dT", filenum);
	printf("\tfile2array file %s...\n", file_name);
	file2array(file_name, T);

	sprintf(file_name, "%dA", filenum);
	printf("\tfile2array file %s...\n", file_name);
	file2array(file_name, A);
	
	*nPoints = filecount(file_name);
	
	//panning
	//TODO x and y panning for 8 chan
	char xpanTime[64];
	char xpanStep[64];
	sprintf(xpanTime, "%dxpantime", filenum);
	sprintf(xpanStep, "%dxpanstep", filenum);
	
	int xpanPoints = filecount(xpanTime);
	
	char ypanTime[64];
	char ypanStep[64];
	sprintf(ypanTime, "%dypantime", filenum);
	sprintf(ypanStep, "%dypanstep", filenum);
	
	int ypanPoints = filecount(ypanTime);
	
	double *xX;
	double *xY;
	double *yX;
	double *yY;
	xX = (double *)Malloc(xpanPoints * sizeof(double));
	xY = (double *)Malloc(xpanPoints * sizeof(double));
	yX = (double *)Malloc(ypanPoints * sizeof(double));
	yY = (double *)Malloc(ypanPoints * sizeof(double));
	file2array(xpanTime, xX);
	file2array(xpanStep, xY);
	file2array(ypanTime, yX);
	file2array(ypanStep, yY);
	
	int x, xIndx, xrange;
	double yrange, m, m_inc;
	double sr = 44100;
	int samples = toneLength * sr; 
	//linear panning
	double x1, y1, x2, y2;
	
	x1 = xX[0];
	y1 = xY[0];
	xIndx = 0;

	for(i = 1; i < xpanPoints; i++){
		x2 = xX[i];
		y2 = xY[i];

		xrange = (x2 - x1) * samples;
		yrange = (y2 - y1);

		m = 0;
		m_inc = 1.0 / (xrange - 1.0);

		for(x = xIndx ; x < xIndx+xrange; x++){
			xpan[x] = y1 +(m * yrange);
			m = m + m_inc;
			//printf("xpan[x] = %f\n", xpan[x]);
		}

		xIndx = xIndx + xrange;

		y1 = y2;
		x1 = x2;
	}
	
	x1 = yX[0];
	y1 = yY[0];
	xIndx = 0;
	
	for(i = 1; i < ypanPoints; i++){
		x2 = yX[i];
		y2 = yY[i];

		xrange = (x2 - x1) * samples;
		yrange = (y2 - y1);

		m = 0;
		m_inc = 1.0 / (xrange - 1.0);

		for(x = xIndx ; x < xIndx+xrange; x++){
			ypan[x] = y1 +(m * yrange);
			m = m + m_inc;
			//printf("ypan[x] = %f\n", ypan[x]);
		}

		xIndx = xIndx + xrange;

		y1 = y2;
		x1 = x2;
	}
	Free(xX);
	Free(xY);
	Free(yX);
	Free(yY);
}


void stripchannel(double dB,
		double FMstartTime, //in seconds
		double Carrier[],
		double Modulator[], 
		double PD[], 
		double D[],
		double Amp[],
		int *nFMTones,
		double T[],
		double A[],
		int *nPoints,
		int envType,
		double xpan[],
		double ypan[])
{
	//allocate memory
	
	int i;
	double toneDur = 0;
	int sr = 44100;
	double largestD = 1;
	for(i = 0; i < *nFMTones ;i++)
		if(D[i] >= largestD)
			largestD = D[i];
	
	for(i = 0; i < *nFMTones ;i++)
		toneDur = toneDur + D[i]; 
	
	double toneStartTime = FMstartTime;
	int c=0;
	int startc = 0;
	double panCount = 0;
	
	short *Asound;
	short *Bsound;
	short *Csound;
	short *Dsound;
	short *Esound;
	short *Fsound;
	short *Gsound;
	short *Hsound;
	
	int a = 0;
	printf("////////Writing Fm tone...\n");
	for(i = 0 ; i < *nFMTones ; i++ ){
		printf("\tAllocating memory...");
		Asound = (short *)Malloc( largestD * sr * sizeof(short) );
		Bsound = (short *)Malloc( largestD * sr * sizeof(short) );
		Csound = (short *)Malloc( largestD * sr * sizeof(short) );
		Dsound = (short *)Malloc( largestD * sr * sizeof(short) );
		Esound = (short *)Malloc( largestD * sr * sizeof(short) );
		Fsound = (short *)Malloc( largestD * sr * sizeof(short) );
		Gsound = (short *)Malloc( largestD * sr * sizeof(short) );
		Hsound = (short *)Malloc( largestD * sr * sizeof(short) );
		printf("Done.\n");
		
		printf("\t\ti = %d\n", i);
		printf("\t\tD[i] = %f\n ", D[i]);
		fm(Asound, D[i], sr, dB * Amp[i], 
				Carrier[i%*nFMTones], Modulator[i%*nFMTones], PD[i%*nFMTones]);
		fm(Bsound, D[i], sr, dB * Amp[i], 
				Carrier[i%*nFMTones], Modulator[i%*nFMTones], PD[i%*nFMTones]);
		fm(Csound, D[i], sr, dB * Amp[i], 
				Carrier[i%*nFMTones], Modulator[i%*nFMTones], PD[i%*nFMTones]);
		fm(Dsound, D[i], sr, dB * Amp[i], 
				Carrier[i%*nFMTones], Modulator[i%*nFMTones], PD[i%*nFMTones]);
		fm(Esound, D[i], sr, dB * Amp[i], 
				Carrier[i%*nFMTones], Modulator[i%*nFMTones], PD[i%*nFMTones]);
		fm(Fsound, D[i], sr, dB * Amp[i], 
				Carrier[i%*nFMTones], Modulator[i%*nFMTones], PD[i%*nFMTones]);
		fm(Gsound, D[i], sr, dB * Amp[i], 
				Carrier[i%*nFMTones], Modulator[i%*nFMTones], PD[i%*nFMTones]);
		fm(Hsound, D[i], sr, dB * Amp[i], 
				Carrier[i%*nFMTones], Modulator[i%*nFMTones], PD[i%*nFMTones]);
	
		printf("\tRunning adsr2()...");	
		//envelope
		adsr2(Asound, D[i], sr, T, A, *nPoints, envType);
		adsr2(Bsound, D[i], sr, T, A, *nPoints, envType);
		adsr2(Csound, D[i], sr, T, A, *nPoints, envType);
		adsr2(Dsound, D[i], sr, T, A, *nPoints, envType);
		adsr2(Esound, D[i], sr, T, A, *nPoints, envType);
		adsr2(Fsound, D[i], sr, T, A, *nPoints, envType);
		adsr2(Gsound, D[i], sr, T, A, *nPoints, envType);
		adsr2(Hsound, D[i], sr, T, A, *nPoints, envType);
		printf("Done.\n");
		//panning
		printf("\tRunning panning algorithms...");
		a = 0;
		//double *starthree = 3.0;
		printf("startc = %d\n", startc);
	///*
		for(c = startc; c < startc + (D[i] * sr) ; c++){
			//printf("\tpan[c] = %f\n", pan[c]);
			//printf("%d\n", c);
			//printf("%f\n", ypan[c]);
			
			//xpan[c] *= 3.0;
			ypan[c] *= 3.0;
			
			if(xpan[c] > 1.0){
				xpan[c] = 1.0;
			}
			if(xpan[c] < 0.0){
				xpan[c] = 0.0;
			}
			
			if(ypan[c] > 3.0){
				ypan[c] = 3.0;
			}
			if(xpan[c] < 0.0){
				ypan[c] = 0.0;
			}
		
			//printf("a= %d , c= %d\n",a,c);
			// TODO x and y panning for 8 chan
			//equalpowerpan = 1.0/(sqrt(pow(pan[c],2.0)+pow((1.0 - pan[c]),2.0)));
			//x axis panning	
			Asound[a] *= (1.0 - xpan[c]);
			Fsound[a] *= (1.0 - xpan[c]);
			Gsound[a] *= (1.0 - xpan[c]);
			Hsound[a] *= (1.0 - xpan[c]);
			
			Bsound[a] *= xpan[c];
			Csound[a] *= xpan[c];
			Dsound[a] *= xpan[c];
			Esound[a] *= xpan[c];
		/*	
			if(xpan[c] < 1.0){
			
				Bsound[a] *= 0.0;
				Esound[a] *= 0.0;
				Csound[a] *= 0.0;
				Dsound[a] *= 0.0;
				
				Hsound[a]  *= (1.0 - xpan[c]);
				Gsound[a]  *= (1.0 - xpan[c]);
			
				Asound[a]  *= xpan[c];
				Fsound[a]  *= xpan[c];
			
			}else if(xpan[c] >= 1.0 && xpan[c] <= 2.0){
				
				Hsound[a]  *= 0.0;
				Gsound[a]  *= 0.0;
				Csound[a]  *= 0.0;
				Dsound[a]  *= 0.0;
			
				Asound[a]  *= (1.0 - (xpan[c] - 1.0));
				Fsound[a]  *= (1.0 - (xpan[c] - 1.0));

				Bsound[a]  *= (xpan[c] - 1.0);
				Esound[a]  *= (xpan[c] - 1.0);
			
			} else if(xpan[c] > 2.0){
				
				Hsound[a]  *= 0.0;
				Gsound[a]  *= 0.0;
				Asound[a]  *= 0.0;
				Fsound[a]  *= 0.0;
			
				Bsound[a]  *= (1.0 - (xpan[c] - 2.0));
				Esound[a]  *= (1.0 - (xpan[c] - 2.0));

				Csound[a]  *= (xpan[c] - 2.0);
				Dsound[a]  *= (xpan[c] - 2.0);
			
			}
		*/	
			//y axis panning
			
			if(ypan[c] < 1.0){
			
				Hsound[a]  *= 0.0;
				Csound[a]  *= 0.0;
				Asound[a]  *= 0.0;
				Bsound[a]  *= 0.0;
				
				Fsound[a]  *= (1.0 - ypan[c]);
				Esound[a]  *= (1.0 - ypan[c]);
			
				Gsound[a]  *= ypan[c];
				Dsound[a]  *= ypan[c];
			
			}else if(ypan[c] >= 1.0 && ypan[c] <= 2.0){
				
				Fsound[a]  *= 0.0;
				Esound[a]  *= 0.0;
				Asound[a]  *= 0.0;
				Bsound[a]  *= 0.0;
			
				Gsound[a]  *= (1.0 - (ypan[c] - 1.0));
				Dsound[a]  *= (1.0 - (ypan[c] - 1.0));

				Hsound[a]  *= (ypan[c] - 1.0);
				Csound[a]  *= (ypan[c] - 1.0);
			
			} else if(ypan[c] > 2.0){
				
				Gsound[a]  *= 0.0;
				Dsound[a]  *= 0.0;
				Fsound[a]  *= 0.0;
				Esound[a]  *= 0.0;
			
				Hsound[a]  *= (1.0 - (ypan[c] - 2.0));
				Csound[a]  *= (1.0 - (ypan[c] - 2.0));

				Asound[a]  *= (ypan[c] - 2.0);
				Bsound[a]  *= (ypan[c] - 2.0);
			
			}
			a++;
		}
	//*/
		startc = c;
		c = (c + (D[i] * sr));
		panCount += (D[i] * sr); 
		printf("Done.\n");
		printf("\tRunning mixAdd()...");	
		//mixAdd
		printf("\tstarttime = %f\n", toneStartTime);
		mixAdd(Asound,toneStartTime,D[i],sr,0);
		mixAdd(Bsound,toneStartTime,D[i],sr,1);
		mixAdd(Csound,toneStartTime,D[i],sr,2);
		mixAdd(Dsound,toneStartTime,D[i],sr,3);
		mixAdd(Esound,toneStartTime,D[i],sr,4);
		mixAdd(Fsound,toneStartTime,D[i],sr,5);
		mixAdd(Gsound,toneStartTime,D[i],sr,6);
		mixAdd(Hsound,toneStartTime,D[i],sr,7);
		
		printf("Done.\n");
		toneStartTime = toneStartTime + D[i];
		printf("\tFreeing memory...\n");
		Free(Asound);
		Free(Bsound);
		Free(Csound);
		Free(Dsound);
		Free(Esound);
		Free(Fsound);
		Free(Gsound);
		Free(Hsound);
		printf("\tMemory freed.\n");	
	}
	printf("////////Done.\n");	
}
