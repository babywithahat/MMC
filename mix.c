#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sf.h>
#include <sfmisc.h>

static double *mixBuffer ;		// global (for this file) *mixBuffer
static int totalFrames ;		// global (for this file) number of frames
static int totalChannels ;		// global (for this file) number of channels

/*
 * mixSetup(totalDuration, sr, nchan) 
 * initializes the output buffer, and saves the variables
 */
void mixSetup(double totalDuration, int sr, int nchan)
{
	extern double *mixBuffer ;
	extern int totalFrames ;
	extern int totalChannels ;

	totalChannels = nchan ;				// save values for nchan & frames
	totalFrames = totalDuration * sr ;

	// allocate memory for 'mixBuffer'
	mixBuffer = (double *)Malloc(totalFrames * totalChannels * sizeof(double)) ;

	return ;
}

/*
 * mixAdd(input, startTime, dur, sr, thisChannel)
 * adds the 'input' sound to the 'mixBuffer' 
 */

void mixAdd(short *input, double startTime, double dur, int sr, int thisChannel)
{
	extern double *mixBuffer ;
	extern int totalChannels ;
	extern int totalFrames ;

	static int functionCounter = 0 ;
	int start, frames, t, n ;

	functionCounter++ ;		// increment the function counter

	// make sure the value for 'thisChannel' is acceptable
	if ( thisChannel >= totalChannels || thisChannel < 0 ) {
		fprintf(stderr,
			"addToMix Errot(%d): Bad value for 'thisChannel'\n",
				functionCounter);
		fprintf(stderr,"\t'totalChannels: %d  'thisChannel': %d\n",
				totalChannels, thisChannel) ;
		exit(EXIT_FAILURE) ;
	}

	start = startTime * sr ;		// set the start time for the new sound
	frames = dur * sr ;			// set the length of the new sound

	// make sure the new sound fits in the mixBuffer
	if ( start + frames > totalFrames ) {
		fprintf(stderr,
			"addToMix Error(%d): 'start' + 'frames' exceeds mixBuffer\n", 
				functionCounter);
		fprintf(stderr,"\t'start' + 'frames' = %d\n", start + frames) ;
		fprintf(stderr,"|t'totalFrames' = %d\n", totalFrames) ;
		exit(EXIT_FAILURE) ;
	}

	t = start * totalChannels + thisChannel ; // start time for the output mix

	for ( n = 0 ; n < frames ; n++ ) { // add the new sound to the output mix
		mixBuffer[t] = mixBuffer[t] + input[n] ;
		t += totalChannels ;
	}

	return ;
}

void mixScale(void)		// scale the output buffer to be within 16 bits
{
	double xmax, xmin, scale ;
	const double ymax = pow(2,15) - 1;		// max for 16-bits: 32767
	const double ymin = -pow(2,15) ;		// min for 16-bits: -32768
	int t ;

	extern double *mixBuffer ;
	extern int totalFrames ;
	extern int totalChannels ;

	// find the current min/max values in the mixBuffer
	xmax = xmin = mixBuffer[0] ;
	for ( t = 0 ; t < totalFrames * totalChannels ; t++ ) {
		if ( mixBuffer[t] > xmax ) 
			xmax = mixBuffer[t] ;
		else if ( mixBuffer[t] < xmin ) 
			xmin = mixBuffer[t] ;
	}

	// set the scale factor
	scale = ( ymax - ymin ) / ( xmax - xmin ) ;

	if ( scale >= 1 )	// no scaling needs: values are within 16-bits
		return ;

	// scale all the values in the mixBuffer ;
	for ( t = 0 ; t < totalFrames * totalChannels ; t++ )
		mixBuffer[t] = scale * ( mixBuffer[t] - xmin ) + ymin ;

	return ;
}

short *mixConvertTo16Bits(void)
{
	extern double *mixBuffer ;
	extern int totalFrames ;
	extern int totalChannels ;
	
	short *output ;
	int t ;

	mixScale() ;	// scale all the values to 16-bits

	// allocate memory for the output buffer
	output = (short *)Malloc(totalFrames * totalChannels * sizeof(short)) ;

	// copy all the values in 'mixBuffer' to 'output'
	for ( t = 0 ; t < totalFrames * totalChannels ; t++ )
		output[t] = mixBuffer[t] ;

	Free(mixBuffer) ;	// release memory

	return output ;
}
