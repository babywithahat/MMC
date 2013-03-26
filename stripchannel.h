
extern int filecount(char* file_name);
extern void file2array(char* file_name, double buffer[]);

extern double durOfTone(int file_letter,double startTime);

extern void allocSetup(char file_letter,
		int *mallocCar,
		int *mallocMod,
		int *mallocPD,
		int *mallocD,
		int *mallocAmp,
		int *mallocT,
		int *mallocA,
		int *mallocPan);

extern void loadarray(char file_letter,
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
		double ypan[]);

extern void stripchannel(double db,
		double FMstartTime,
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
		double ypan[]);
