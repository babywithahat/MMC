
void concatenate(short output[], short input[], double st, double dur, int sr)
{
	int start, length, n, i ;

	start = st * sr ;		// start time of tone in samples
	length = dur * sr ;		// duration of tone in samples

	// using 2 indexes (i and n), copy 'input' into 'output'
	for ( n = start, i = 0 ; n < start + length ; n++, i++ )
		output[n] = input[i] ;

	return ;
}

