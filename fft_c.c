/******************************************************************************
   This function computes the FFT and IFFT.
   dir: direction of transformation, 1-FFT, 0-IFFT 
   pow: number of data points in power, total data points = 2^pow
   re: real part of the input data array
   im: imaginary part of the input data array
******************************************************************************/
void fft_c(unsigned int dir,
		   unsigned int pow,
		   double *re,
           double *im)
{
    int num = 1, num_hf; // number of points
	int i = 0, j = 0;
	int a = 0, b = 0, k = 0;
	int l1, l2 = 1;
	
	double temp_x, temp_y; // temporarily saves the input data
	double c1 = -1.0, c2 = 0.0;
	double u1, u2, temp_u1;
	double tmp1, tmp2;
	
	/* number of points of the input data */
	for (i = 0; i < pow; i++) 
    {
		num *= 2;
	}

	/* Bit reverse */
	num_hf = num >> 1;
	b = 0;

	for (a = 0; a < num - 1; a++)
	{
		if (a < b)
		{
			// exchange the data
			temp_x = re[a];
			temp_y = im[a];
			re[a] = re[b];
			im[a] = im[b];
			re[b] = temp_x;
			im[b] = temp_y;
		}

		// find position b
		k = num_hf;
		while (k <= b)
		{
			b -= k;
			k >>= 1;
		}

		b += k;
    }

    /* Calculates FFT */
	for (i = 0; i < pow; i++)
	{
		u1 = 1.0; 
		u2 = 0.0;

		l1 = l2;
		l2 <<= 1;

		for (j = 0; j < l1; j++)
		{
			for (a = j; a < num; a += l2)
			{
				b = a + l1;
				tmp1 = u1 * re[b] - u2 * im[b];
				tmp2 = u1 * im[b] + u2 * re[b];
				re[b] = re[a] - tmp1; 
				im[b] = im[a] - tmp2;
				re[a] += tmp1;
				im[a] += tmp2;
			}

			temp_u1 =  u1 * c1 - u2 * c2;
			u2 = u1 * c2 + u2 * c1;
			u1 = temp_u1;
		}

		c2 = sqrt((1.0 - c1) / 2.0);
		c1 = sqrt((1.0 + c1) / 2.0);

		if (dir == 1) 
		{
			// forward
			c2 = -c2;
		}
	}

	/* Scaling for the inverse transformation */
    if (dir == -1)
	{
		for (i = 0; i < num; i++)
		{
			re[i] /= num;
			im[i] /= num;
		}
    }
   
    return;
}

