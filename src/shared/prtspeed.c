#include <stdio.h>

#define max(a,b)( ((a) > (b)) ? (a) : (b) )

void prtspeed(int m, int l, int n, double time, int ok) {
	double		speed;

	speed = 1.0E-6 * (2 * (double)m * (double)n * (double)l) / (max(time, 1.0E-9));

	printf(" %6d| %6d| %6d| %11.4f| %11.4e| ", m, l, n, time, speed);

	if (ok == 0)
		printf(" T |\n");
	else
		printf(" F |\n");

	fflush(stdout);
}
