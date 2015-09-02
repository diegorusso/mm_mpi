#include <time.h>
#include <stddef.h>
#include <stdio.h>

void state(char *progname) {
	time_t		now;

	printf(" Single-CPU benchmark (C-version), program %s\n", progname);
	printf(" -----------------------------------------------------------");
	printf("------------\n");

	char		machin    [] = "CPU type";
	char		memory    [] = "XGB of memory per node";
	char		compil    [] = "compile version";
	char		option    [] = "flags";
	char		os        [] = "-";
	char		runby     [] = "Diego Russo (me@diegor.it)";
	char		comins    [] = "UniPG";

	printf(" Computer                  %s\n", machin);
	printf(" Memory size               %s\n", memory);
	printf(" Compiler version          %s\n", compil);
	printf(" Compiler options          %s\n", option);
	printf(" Operating System version  %s\n", os);
	printf(" Working precision         64 bits\n");
	printf(" Run by                    %s\n", runby);
	printf(" Company/Institute         %s\n", comins);
	printf("\n This program is run at: %s", asctime(localtime(&now)));
	printf(" -----------------------------------------------------------");
	printf("------------\n");

	now = time(NULL);
}
