
extern void mixSetup(double totalDuration, int sr, int nchan) ;
extern void mixAdd(short *input, 
		double startTime, double dur, int sr, int thisChannel) ;
extern short *mixConvertTo16Bits(void) ;

