void debug_printf(const char* func_name, int mpi_rank, const char* fmt, ...);

void header(void);

void results(int m, int l, int n, double time, int ok);

void footer(void);
