#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sf.h>
#include <sfmisc.h>
#include <sfmalloc.h>


double *exponentialEnv0(
		double dur,			// duration of envelope in seconds
		int sr,				// sampling rate
		double X[],			// array of X values
		double Y[],			// array of Y values
		int N)				// number of X,Y pairs
{
	int i,					// counter
		samples ;			// length of envelope in samples
	double *env ;			// array for the envelope
	double x1, x2, y1, y2 ;	// variables for x,y pairs

	double m, m_inc ;		// slope, slope increment
	int x, xrange, xIndx ;	// integerst to count 'x' values

	double logRatio,		// log of the ratio of y1 and y2
		   yrange ;			// range of y values for one segment

	/*
	 *	function that will be used to calculate the
	 *	exponential curves between the points:  
	 *	f(x) = y1 + (y2-y1) * 
	 *			[ (1 - e^(x * log(y2/y1))) / (1 - e^(log(y2/y1))) ]
	 *	where
	 *		x : time, moving from 0 to 1.0
	*			when x == 0, y = y1; when x == 1, y = y2	
	*		y1, y2 : initial and final values for y
	 */

	// calculate the duration of the sound in samples,
	// and allocate memory
	samples = dur * sr ;
	env = (double *)Malloc(samples * sizeof(double)) ;

	x1 = X[0] ;		// initialize x1, y1
	y1 = Y[0] ;
	if ( y1 == 0 )		// make sure y1 is not zero
		y1 = 0.000001 ;	

	xIndx = 0 ;		// initialize xIndx

	for ( i = 1 ; i < N ; i++ ) { // loop through the xy points

		x2 = X[i] ;		// set values for x2, y2
		y2 = Y[i] ;
		if ( y2 == 0 )	// make sure y2 is not zero
			y2 = 0.000001 ;	

		// set xrange and yrange
		xrange = (int)((x2 - x1) * samples + 0.5) ; // roundoff
		yrange = y2 - y1 ;

		if ( xIndx + xrange > samples )	// check for roundoff error
			xrange = samples - xIndx ;

		m = 0 ;							// initialize the slope
		m_inc = 1.0 / ( xrange - 1.0 ) ;	// slope increment

		logRatio = log(y2/y1) ;			// set the log of the ratio
		yrange = yrange / (1.0 - exp(logRatio)) ;	// scale yrange

		// create exponential segment 
		for ( x = xIndx ; x < xIndx + xrange ; x++ ) {
			env[x] = y1 + yrange * ( 1 - exp( m * logRatio ) ) ;
			m = m + m_inc ;
		}

		xIndx = xIndx + xrange ;	// increment the xIndex

		x1 = x2 ;					// x2,y2 -> x1,y1
		y1 = y2 ;

	}

	return env ;
}

double *exponentialEnv(
		double dur,			// duration of envelope in seconds
		int sr,				// sampling rate
		double X[],			// array of X values
		double Y[],			// array of Y values
		int N,				// number of X,Y pairs
		int alpha)			// number of exponential time constants used
{
	int i,					// counter
		samples ;			// length of envelope in samples
	double *env ;			// array for the envelope
	double x1, x2, y1, y2 ;	// variables for x,y pairs

	double m, m_inc ;		// slope, slope increment
	int x, xrange, xIndx ;	// integerst to count 'x' values

	double logRatio,		// log of the ratio of y1 and y2
		   yrange ;			// range of y values for one segment

	/*
	 *	function that will be used to calculate the
	 *	exponential curves between the points:  
	 *	f(x) = y1 + (y2-y1) * 
	 *			[ (1 - e^(x * alpha)) / ( 1 - e^alpha ) ]
	 *	where
	 *		x : time, moving from 0 to 1.0
	*			when x == 0, y = y1; when x == 1, y = y2	
	*		y1, y2 : initial and final values for y
	 */

	// calculate the duration of the sound in samples,
	// and allocate memory
	samples = dur * sr ;
	env = (double *)Malloc(samples * sizeof(double)) ;

	x1 = X[0] ;		// initialize x1, y1
	y1 = Y[0] ;
	if ( y1 == 0 )		// make sure y1 is not zero
		y1 = 0.000001 ;	

	xIndx = 0 ;		// initialize xIndx

	for ( i = 1 ; i < N ; i++ ) { // loop through the xy points

		x2 = X[i] ;		// set values for x2, y2
		y2 = Y[i] ;
		if ( y2 == 0 )	// make sure y2 is not zero
			y2 = 0.000001 ;	

		// set xrange and yrange
		xrange = (int)((x2 - x1) * samples + 0.5) ; // roundoff
		yrange = y2 - y1 ;

		if ( xIndx + xrange > samples )	// check for roundoff error
			xrange = samples - xIndx ;

		m = 0 ;							// initialize the slope
		m_inc = 1.0 / ( xrange - 1.0 ) ;	// slope increment

		// logRatio = log(y2/y1) ;			// set the log of the ratio
		logRatio = alpha ;
		yrange = yrange / (1.0 - exp(logRatio)) ;	// scale yrange

		// create exponential segment 
		for ( x = xIndx ; x < xIndx + xrange ; x++ ) {
			env[x] = y1 + yrange * ( 1 - exp( m * logRatio ) ) ;
			m = m + m_inc ;
		}

		xIndx = xIndx + xrange ;	// increment the xIndex

		x1 = x2 ;					// x2,y2 -> x1,y1
		y1 = y2 ;

	}

	return env ;
}

double *linearEnv(double dur, int sr, double X[], double Y[], int N)
{
	double *env ;
	int i, x, xIndx, xrange ;
	double x1, y1, x2, y2 ;
	double yrange, m, m_inc ;
	int samples ;

	samples = dur * sr ;
    env = (double *)Malloc(samples * sizeof(double)) ;

    x1 = X[0] ;			// initial X and Y
    y1 = Y[0] ;	

    xIndx = 0 ;				// starting value for index

    for ( i = 1 ; i < N ; i++ ) {		// Create the envelope shape

		x2 = X[i] ;						// set x2 & y2
		y2 = Y[i] ;					

        xrange = (x2 - x1) * samples ;		// xrange (in samples)
		yrange = (y2 - y1) ;				// yrange (between 0-1)

		m = 0 ;								// initialize slope
		m_inc = 1.0 / ( xrange - 1.0 ) ;	// slope increment

		// create linear segment 
        for ( x = xIndx ; x < xIndx+xrange ; x++ ) {
            env[x] = y1 + (m * yrange) ;
			m = m + m_inc ;
		}
		
        xIndx = xIndx + xrange ;     // update the x index

        y1 = y2 ;					// previous y2/x2 become next y1/x1
        x1 = x2 ;     
    }

	return env ;
}

void adsr2(short sound[], double dur, int sr, 
		double T[], double A[], int N, int flag)
{
    double *envelope ;		// array to create the envelope
    int t, samples ;

	if ( flag == 0 )
		envelope = linearEnv(dur, sr, T, A, N) ;
	else
		envelope = exponentialEnv0(dur, sr, T, A, N) ;

    samples = dur * sr ;    // duration in samples

	// apply the envelope to the sound
    for ( t = 0 ; t < samples ; t++ )	
        sound[t] = envelope[t] * sound[t] ;

    Free(envelope) ;    // free the memory used for 'envelope'

    return ;
}

void copy(short *output, short *input, double d, int sr)
{
	memcpy(output, input, d * sr * sizeof(short)) ;

	return ;
}

#ifdef TEST

int main(void)
{
	int N = 10 ;
	short *wave, *input, *output ;
	double st = 0 , d = 1, f = 440, db = 96 ;
	int samples, i ;
	int nPoints = 5, sr = 44100 ;
	double T[] = {0.0, 0.1, 0.25, 0.75, 1.0} ;
	double A[] = {0.0, 1.0, 0.75, 0.50, 0.0} ;

	samples = d * sr ;

	wave = (short *)Malloc(samples * sizeof(short)) ;
	input = (short *)Malloc(samples * sizeof(short)) ;
	output = (short *)Malloc(N * samples * sizeof(short)) ;


	Sine(wave, st, d, sr, f, db) ;

	for ( i = 0 ; i < N ; i++ ) {
		copy(input, wave, d, sr) ;
		adsr2(input, d, sr, T, A, nPoints, -((i+1))) ;
		concatenate(output, input, st, d, sr) ;
		st += d ;
	}

	sfsave("atest.wav", output, N*d, sr, 1) ;

	Free(wave) ;
	Free(input) ;
	Free(output) ;

	return 0 ;
}
#endif
