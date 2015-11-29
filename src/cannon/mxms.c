
void mxms(int m, int l, int n, double *a, double *b, double *c) {
	int		i         , j, k;
	double		t;

#if defined __OPENMP_OUTER

#pragma omp parallel for default(shared) private(t,k,j) schedule(runtime)
	for (i = 0; i < m; i++) {
		for (j = 0; j < l; j++) {
			t = a[i * l + j];
			for (k = 0; k < n; k++) {
				c[i * n + k] = c[i * n + k] + t * b[j * n + k];
			}
		}
	}

#elif defined __OPENMP_MIDDLE

	for (i = 0; i < m; i++) {
#pragma omp parallel for default(shared) private(t,j) schedule(runtime)
		for (k = 0; k < n; k++) {
			t = c[i * n + k];
			for (j = 0; j < l; j++) {
				t = t + a[i * l + j] * b[j * n + k];
			}
			c[i * n + k] = t;
		}
	}


#elif defined __OPENMP_INNER

	for (i = 0; i < m; i++) {
		for (j = 0; j < l; j++) {
			t = a[i * l + j];
#pragma omp parallel for default(shared) firstprivate(t) schedule(runtime)
			for (k = 0; k < n; k++) {
				c[i * n + k] = c[i * n + k] + t * b[j * n + k];
			}
		}
	}

#elif defined __OPENMP_NESTED

#pragma omp parallel for default(shared) private(t,k,j) num_threads(2) schedule(runtime)
	for (i = 0; i < m; i++) {
#pragma omp parallel for default(shared) private(t,j) num_threads(2) schedule(runtime)
		for (k = 0; k < n; k++) {
			t = c[i * n + k];
			for (j = 0; j < l; j++) {
				t = t + a[i * l + j] * b[j * n + k];
			}
			c[i * n + k] = t;
		}
	}


#else

	for (i = 0; i < m; i++)
		for (k = 0; k < n; k++) {
			t = c[i * n + k];
			for (j = 0; j < l; j++)
				t = t + a[i * l + j] * b[j * n + k];
			c[i * n + k] = t;
		}

#endif
}
