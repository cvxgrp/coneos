#ifndef RUN_coneOS_H_GUARD
#define RUN_coneOS_H_GUARD

int main(int argc, char **argv);
void freeData(Data * d, Cone * k);
int read_in_data(FILE * fp,Data * d, Cone * k);
int open_file(int argc, char ** argv, int idx, char * default_file, FILE ** fb);
void freeSol(Sol *sol);

#endif
