#include "gft_heap64f.h"

namespace gft{
  namespace Heap64f{

    void SetRemovalPolicy(Heap64f *H, char policy){
      if(H->removal_policy != policy){
	H->removal_policy = policy;
      }
    }

    void GoUp(Heap64f *H, int i) {
      int j = HEAP_DAD(i);
      
      if(H->removal_policy == MINVALUE){
	
	while ((i > 0) && (H->cost[H->pixel[j]] > H->cost[H->pixel[i]])) {
	  gft::SwapInt(&H->pixel[j], &H->pixel[i]);
	  H->pos[H->pixel[i]] = i;
	  H->pos[H->pixel[j]] = j;
	  i = j;
	  j = HEAP_DAD(i);
	}
      }
      else{ /* removal_policy == MAXVALUE */
	
	while ((i > 0) && (H->cost[H->pixel[j]] < H->cost[H->pixel[i]])) {
	  gft::SwapInt(&H->pixel[j], &H->pixel[i]);
	  H->pos[H->pixel[i]] = i;
	  H->pos[H->pixel[j]] = j;
	  i = j;
	  j = HEAP_DAD(i);
	}
      }
    }

    void GoDown(Heap64f *H, int i) {
      int j, left = HEAP_LEFTSON(i), right = HEAP_RIGHTSON(i);
      
      j = i;
      if(H->removal_policy == MINVALUE){
	
	if ((left <= H->last) && 
	    (H->cost[H->pixel[left]] < H->cost[H->pixel[i]]))
	  j = left;
	if ((right <= H->last) && 
	    (H->cost[H->pixel[right]] < H->cost[H->pixel[j]]))
	  j = right;
      }
      else{ /* removal_policy == MAXVALUE */
	
	if ((left <= H->last) && 
	    (H->cost[H->pixel[left]] > H->cost[H->pixel[i]]))
	  j = left;
	if ((right <= H->last) && 
	    (H->cost[H->pixel[right]] > H->cost[H->pixel[j]]))
	  j = right;
      }
      
      if(j != i) {
	gft::SwapInt(&H->pixel[j], &H->pixel[i]);
	H->pos[H->pixel[i]] = i;
	H->pos[H->pixel[j]] = j;
	GoDown(H, j);
      }
    }

    char IsFull(Heap64f *H) {
      if (H->last == (H->n - 1))
	return 1;
      else
	return 0;
    }
    
    char IsEmpty(Heap64f *H) {
      if (H->last == -1){
	Reset(H); 
	return 1;
      }else
	return 0;
    }
    
    Heap64f *Create(int n, double *cost) {
      Heap64f *H = NULL;
      int i;
      
      if (cost == NULL) {
	fprintf(stdout,"Cannot create heap without cost map in Heap64f::Create");
	return NULL;
      }
      
      H = (Heap64f *) malloc(sizeof(Heap64f));
      if (H != NULL) {
	H->n       = n;
	H->cost    = cost;
	H->color   = (char *) malloc(sizeof(char) * n);
	H->pixel   = (int *) malloc(sizeof(int) * n);
	H->pos     = (int *) malloc(sizeof(int) * n);
	H->last    = -1;
	H->removal_policy = MINVALUE;
	if (H->color == NULL || H->pos == NULL || H->pixel == NULL)
	  gft::Error((char *)MSG1,(char *)"Heap64f::Create");
	for (i = 0; i < H->n; i++) {
	  H->color[i] = WHITE;
	  H->pos[i]   = -1;
	  H->pixel[i] = -1;
	}    
      } 
      else
	gft::Error((char *)MSG1,(char *)"Heap64f::Create");

      return H;
    }

    void Destroy(Heap64f **H) {
      Heap64f *aux = *H;
      if (aux != NULL) {
	if (aux->pixel != NULL) free(aux->pixel);
	if (aux->color != NULL) free(aux->color);
	if (aux->pos != NULL)   free(aux->pos);
	free(aux);
	*H = NULL;
      }
    }
    
    char Insert(Heap64f *H, int pixel) {
      if (!IsFull(H)) {
	H->last++;
	H->pixel[H->last] = pixel;
	H->color[pixel]   = GRAY;
	H->pos[pixel]     = H->last;
	GoUp(H, H->last); 
	return 1;
      } else 
	return 0;
    }

    char Remove(Heap64f *H, int *pixel) {
      if (!IsEmpty(H)) {
	*pixel = H->pixel[0];
	H->pos[*pixel]   = -1;
	H->color[*pixel] = BLACK;
	H->pixel[0]      = H->pixel[H->last];
	H->pos[H->pixel[0]] = 0;
	H->pixel[H->last] = -1;
	H->last--;
	GoDown(H, 0);
	return 1;
      } else 
	return 0;
    }


    void Update(Heap64f *H, int p, double value){
      H->cost[p] = value;
      
      if (H->color[p] == BLACK) printf("ferrou\n");
      
      if(H->color[p] == WHITE)
	Insert(H, p);
      else
	GoUp(H, H->pos[p]);
    }
    
    void Reset(Heap64f *H){
      int i;
      
      for (i=0; i < H->n; i++) {
	H->color[i] = WHITE;
	H->pos[i]   = -1;
	H->pixel[i] = -1;
      }
      H->last = -1;
    }
    

  } /*end Heap64f namespace*/
} /*end gft namespace*/


