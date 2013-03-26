extern void adsr2(short sound[], double dur, int sr, 
		double T[], double A[], int N, int flag) ;

extern double 
	*linearEnv(double dur, int sr, double X[], double Y[], int N) ;

extern double *exponentialEnv0(
		double dur,			// duration of envelope in seconds
		int sr,				// sampling rate
		double X[],			// array of X values
		double Y[],			// array of Y values
		int N);				// number of X,Y pairs
