
#include "gft.h"



int main(int argc, char **argv){
  gft::FileList::FileList *L, *C;
  char *file1, *file2;
  char command[1024];
  char tmp[512];
  float per,s1,s2;
  FILE *fp1;
  int i,method;

  if(argc < 4){
    fprintf(stdout,"usage:\n");
    fprintf(stdout,"runexp <image_list> <config_list> <method>\n");
    fprintf(stdout,"method: 0 ... Riverbed\n");
    fprintf(stdout,"        1 ... LiveWire\n");
    fprintf(stdout,"        2 ... G-wire\n");
    fprintf(stdout,"        3 ... Bandeirantes\n");
    exit(0);
  }
  
  L = gft::FileList::Read(argv[1]);
  C = gft::FileList::Read(argv[2]);
  method = atoi(argv[3]);
  
  s1 = s2 = 0.0;
  for(i = 0; i < L->n; i++){
    file1 = gft::FileList::GetFile(L, i);
    file2 = gft::FileList::GetFile(C, i);
    sprintf(command, "./proj1 %s %s %d > out.txt", file1, file2, method);
    system(command);
    
    fp1 = fopen("out.txt", "r");
    if(fp1 == NULL)
      gft::Error("Error in file out.txt","runexp");

    fscanf(fp1, "%s %f%%", tmp, &per);
    s1 += per;
    
    strcpy(tmp, file1);
    gft::FileList::RemoveFileDirectory(tmp);
    printf("%s; %f; ", tmp, per);
    
    fscanf(fp1, "%s %f%%", tmp, &per);
    s2 += per;
    
    printf("%f\n", per);
    
    fclose(fp1);
  }
  printf("mean; %f; %f\n", s1/L->n, s2/L->n);
  
  gft::FileList::Destroy(&L);
  gft::FileList::Destroy(&C);

  return 0;
}

