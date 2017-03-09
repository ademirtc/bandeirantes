#ifndef _GFT_HEAP64F_H_
#define _GFT_HEAP64F_H_

#include "gft_common.h"
#include "gft_gpqueue_by_Falcao.h"

#include "gft_heap.h"

namespace gft{
  namespace Heap64f{

    typedef struct _heap64f {
      double *cost;
      char *color;
      int *pixel;
      int *pos;
      int last;
      int n;
      char removal_policy; /* 0 is MINVALUE and 1 is MAXVALUE */
    } Heap64f;


    /* Auxiliary Functions */

    //#define HEAP_DAD(i) ((i - 1) / 2)
    //#define HEAP_LEFTSON(i) (2 * i + 1)
    //#define HEAP_RIGHTSON(i) (2 * i + 2)

    void SetRemovalPolicy(Heap64f *H, char policy);
    char IsFull(Heap64f *H);
    char IsEmpty(Heap64f *H);
    Heap64f *Create(int n, double *cost);
    void Destroy(Heap64f **H);
    char Insert(Heap64f *H, int pixel);
    char Remove(Heap64f *H, int *pixel);
    void Update(Heap64f *H, int p, double value);
    void GoUp(Heap64f *H, int i);
    void GoDown(Heap64f *H, int i);
    void Reset(Heap64f *H);

  } //end Heap64f namespace
} //end gft namespace

#endif

