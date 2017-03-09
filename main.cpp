
#include "gft.h"

#include "munkres.h"


#define METHOD_RIVERBED      0
#define METHOD_LIVEWIRE      1
#define METHOD_GWIRE         2
#define METHOD_BANDEIRANTES  3


int get_arc_index(int dx, int dy, gft::AdjRel::AdjRel *A){
  int i;
  for(i = 1; i < A->n; i++){
    if(dx == A->dx[i] && dy == A->dy[i])
      return i;
  }
  return -1;
}


typedef struct _GWireData{
  int n;
  gft::AdjRel::AdjRel *A;
  int *pred;
  double *cost;
} GWireData;


typedef struct _Path{
  int nnodes;
  int *nodes;
  double *cost;
} Path;



double ComputePathCostByBandeirantes(Path *P,
				     gft::Image32::Image32 *img){
  double power = 9.0;
  double alpha = 0.03;
  int ntimes = 16;
  int p1,p2,p1_x,p1_y,p2_x,p2_y,dx,dy,v_x,v_y,nt;
  double ds,edge,Imean,sum,tmp;
  int i,p,q,k,r;
  tmp = 0.0;
  P->cost[P->nnodes-1] = 0.0;
  for(i = P->nnodes-2; i >= 0; i--){
    q = P->nodes[i];
    p = P->nodes[i+1];

    v_x = q%img->ncols;
    v_y = q/img->ncols;
    
    edge = img->data[p] + img->data[q];
    
    nt = ntimes-1;
    k = i+1;
    r = p;
    while(nt > 0 && k < P->nnodes-1){ 
      k++;
      r = P->nodes[k];
      nt--;
    }
    p1 = r;
    
    nt = ntimes;
    sum = img->data[p1];
    r = p1;
    while(nt > 0 && k < P->nnodes-1){ 
      k++;
      r = P->nodes[k];
      sum += img->data[r];
      nt--;
    }
    Imean = sum/(ntimes - nt + 1);
    p2 = r;
    
    p1_x = p1%img->ncols;
    p1_y = p1/img->ncols;
    p2_x = p2%img->ncols;
    p2_y = p2/img->ncols;
    dx = p1_x*2 - p2_x - v_x;
    dy = p1_y*2 - p2_y - v_y;
    ds = sqrt(dx*dx + dy*dy);
    
    edge += fabs(img->data[q] - Imean);
    
    tmp = tmp + pow(alpha*edge, power);
    if(nt == 0) tmp += pow(ds, power);

    P->cost[i] = tmp;  // /(P->nnodes - i);
  }
}



Path *GetIFTPaths(gft::Image32::Image32 *pred,
		  double *cost, int S[]){
  Path *P;
  int j,p,n;
  P = (Path *)calloc(S[0], sizeof(Path));
  for(j = 1; j <= S[0]; j++){
    p = S[j];
    n = 0;
    while(p != NIL){
      n++;
      p = pred->data[p];
    }

    P[j-1].nnodes = n;
    P[j-1].nodes = gft::AllocIntArray(n);
    P[j-1].cost = gft::AllocDoubleArray(n);

    p = S[j];
    n = 0;
    while(p != NIL){
      P[j-1].nodes[n] = p;
      P[j-1].cost[n]  = cost[p];
      n++;
      p = pred->data[p];
    }
  }
  return P;
}


Path *GetGWirePaths(GWireData *data, int S[]){
  gft::AdjRel::AdjRel *A;
  Path *P;
  int j,i,n,p,vp,pmin;
  P = (Path *)calloc(S[0], sizeof(Path));
  A = data->A;
  for(j = 1; j <= S[0]; j++){
    vp = S[j];
    
    pmin = vp*A->n + 1;
    for(i = 1; i < A->n; i++){
      p = vp*A->n + i;
      if(data->cost[p] < data->cost[pmin])
	pmin = p;
    }

    n = 0;
    p = pmin;
    while(p != NIL){
      n++;
      p = data->pred[p];
    }

    P[j-1].nnodes = n;
    P[j-1].nodes = gft::AllocIntArray(n);
    P[j-1].cost = gft::AllocDoubleArray(n);

    n = 0;
    p = pmin;
    while(p != NIL){
      P[j-1].nodes[n] = p / A->n;
      P[j-1].cost[n]  = data->cost[p];
      n++;
      p = data->pred[p];
    }
  }
  return P;
}



GWireData *CreateGWireData(int n, gft::AdjRel::AdjRel *A){
  GWireData *GWD;
  GWD = (GWireData *) calloc(1,sizeof(GWireData));
  if (GWD != NULL){
    GWD->pred = gft::AllocIntArray(n * A->n);
    GWD->cost = gft::AllocDoubleArray(n * A->n);
    GWD->n = n;
    GWD->A = A;
  } else {
    gft::Error((char *)MSG1,(char *)"CreateGWireData");
  }
  return GWD;
}


void DestroyGWireData(GWireData **data){
  if((*data)->pred != NULL)
    gft::FreeIntArray(&((*data)->pred));
  if((*data)->cost != NULL)
    gft::FreeDoubleArray(&((*data)->cost));
  if((*data)->A != NULL)
    gft::AdjRel::Destroy(&((*data)->A));
  *data = NULL;
}


/*G-wire: A Livewire Segmentation Algorithm Based on
  a Generalized Graph Formulation*/
GWireData *GWire(gft::Image32::Image32 *img, int s){
  gft::Heap64f::Heap64f *Q=NULL;
  gft::AdjRel::AdjRel *A;
  GWireData *GWD;
  int *Ainv;
  double *Dist;
  int i,k,p,q,t,n;
  double edge,tmp,dx,dy;
  double alpha=0.01, gamma=0.03;
  gft::Pixel v,w,u;
  int vz, wz, vp, wp, up;
  
  n = img->ncols*img->nrows;
  A = gft::AdjRel::Circular(1.5);
  GWD = CreateGWireData(n, A);

  Q = gft::Heap64f::Create(n*A->n, GWD->cost);
  gft::Heap64f::SetRemovalPolicy(Q, MINVALUE);

  Ainv = gft::AllocIntArray(A->n);
  Dist = gft::AllocDoubleArray(A->n);
  for(i = 1; i < A->n; i++){
    k = get_arc_index(-A->dx[i], -A->dy[i], A);
    Ainv[i] = k;
    Dist[i] = sqrt(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i]);
  }
  
  for(p = 0; p < n*A->n; p++){
    GWD->pred[p] = NIL;
    GWD->cost[p] = DBL_MAX;
  }

  for(i = 1; i < A->n; i++){
    p = s*A->n + i;
    GWD->cost[p] = 0.0;
    gft::Heap64f::Insert(Q, p);
  }
  
  while(!gft::Heap64f::IsEmpty(Q)){
    gft::Heap64f::Remove(Q, &p);
    vp = p / A->n;
    vz = p % A->n;
    v.x = vp%img->ncols;
    v.y = vp/img->ncols;

    u.x = v.x + A->dx[vz];
    u.y = v.y + A->dy[vz];
    up = u.x + u.y*img->ncols;
    if(!gft::Image32::IsValidPixel(img, u.x, u.y))
      continue;
    
    for(i = 1; i < A->n; i++){
      w.x = v.x + A->dx[i];
      w.y = v.y + A->dy[i];

      if(gft::Image32::IsValidPixel(img, w.x, w.y)){
	wp = w.x + img->ncols*w.y;
	k = Ainv[i];
	q = wp*A->n + k;
	
	if(Q->color[q] != BLACK){
	  
	  edge = img->data[vp] + img->data[wp];
	  //edge += abs(img->data[wp] - img->data[up]);
	  
	  tmp = GWD->cost[p] + pow(gamma*edge, 9.0);
	  tmp += alpha*Dist[i];
	  dx = w.x - 2.0*v.x + u.x;
	  dy = w.y - 2.0*v.y + u.y;
	  tmp += pow(sqrt(dx*dx + dy*dy), 9.0);
	  //---------------------------------------
	  
	  if(tmp < GWD->cost[q]){
	    gft::Heap64f::Update(Q, q, tmp);
	    GWD->pred[q] = p;
	  }
	}
      }
    }
  }
  gft::FreeIntArray(&Ainv);
  gft::FreeDoubleArray(&Dist);
  gft::Heap64f::Destroy(&Q);
  return GWD;
}



gft::CImage::CImage *ViewPaths(gft::Image32::Image32 *img,
			       Path *P,
			       int npaths){
  gft::CImage::CImage *paths;
  int i,k,p;
  paths = gft::CImage::Clone(img);

  for(i = 1; i <= npaths; i++){
    for(k = 0; k < P[i-1].nnodes; k++){
      p = P[i-1].nodes[k];
      (paths->C[0])->data[p] = 255;
      (paths->C[1])->data[p] = 255;
      (paths->C[2])->data[p] = 0;
    }
  }
  return paths;
}


int BackPred(gft::Image32::Image32 *pred, int p, int *ntimes){
  while((*ntimes) > 0 && pred->data[p] != NIL){
    p = pred->data[p];
    (*ntimes)--;
  }
  return p;
}


int BackPredMean(gft::Image32::Image32 *img,
		 gft::Image32::Image32 *pred, int p, int *ntimes,
		 double *Imean){
  double sum = img->data[p];
  int n = 1;
  while((*ntimes) > 0 && pred->data[p] != NIL){
    p = pred->data[p];
    sum += img->data[p];
    n++;
    (*ntimes)--;
  }
  *Imean = sum/n;
  return p;
}



gft::Image32::Image32 *Bandeirantes(gft::Image32::Image32 *img,
				    gft::AdjRel::AdjRel *A,
				    int s,
				    double *cost){
  gft::Heap64f::Heap64f *Q=NULL;
  gft::Image32::Image32 *pred;
  int i,p,q,n,p1,p2;
  double edge,tmp,ds,Imean;
  double power = 9.0;
  double alpha = 0.03;
  int ntimes = 16;
  int u_x,u_y,v_x,v_y, nt;
  int p1_x,p1_y,p2_x,p2_y,dx,dy;
  //int nties = 0,ntotal = 0;
  double *Dist;
  
  n    = img->ncols*img->nrows;
  pred = gft::Image32::Create(img);
  Q = gft::Heap64f::Create(n, cost);
  gft::Heap64f::SetRemovalPolicy(Q, MINVALUE);

  Dist = gft::AllocDoubleArray(A->n);
  for(i = 1; i < A->n; i++){
    Dist[i] = sqrt(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i]);
  }
  
  
  gft::Image32::Set(pred, NIL);
  for(p=0; p<n; p++){
    cost[p] = DBL_MAX;
  }
  
  cost[s] = 0.0;
  gft::Heap64f::Insert(Q, s);
  
  while(!gft::Heap64f::IsEmpty(Q)){
    gft::Heap64f::Remove(Q, &p);
    u_x = p%img->ncols; 
    u_y = p/img->ncols; 
    
    for(i=1; i<A->n; i++){
      v_x = u_x + A->dx[i];
      v_y = u_y + A->dy[i];
      if(v_x >= 0 && v_x < img->ncols &&
	 v_y >= 0 && v_y < img->nrows){
	q = v_x + img->ncols*v_y;

	//ntotal++;
	
	if(Q->color[q] != BLACK){
	  
	  edge = img->data[p] + img->data[q];
	  
	  nt = ntimes-1;
	  p1 = BackPred(pred, p,  &nt);
	  nt = ntimes;
	  p2 = BackPredMean(img, pred, p1, &nt, &Imean);
	  p1_x = p1%img->ncols;
	  p1_y = p1/img->ncols;
	  p2_x = p2%img->ncols;
	  p2_y = p2/img->ncols;
	  dx = p1_x*2 - p2_x - v_x;
	  dy = p1_y*2 - p2_y - v_y;
	  ds = sqrt(dx*dx + dy*dy);

	  edge += fabs(img->data[q] - Imean);
	  
	  tmp = cost[p] + pow(alpha*edge, power);
	  if(nt == 0) tmp += pow(ds, power);
	  //tmp += Dist[i];

	  /*
	  if(tmp == cost[q]){
	    nties++;
	  }
	  */

	  if(tmp < cost[q]){
	    gft::Heap64f::Update(Q, q, tmp);
	    pred->data[q] = p;
	  }
	}
      }
    }
  }
  //printf("nties: %f%%\n",((float)nties/(float)ntotal)*100.0);
  gft::Heap64f::Destroy(&Q);
  gft::FreeDoubleArray(&Dist);
  return pred;
}


gft::Image32::Image32 *LiveWire(gft::Image32::Image32 *img,
				gft::AdjRel::AdjRel *A,
				int s,
				double *cost){
  gft::Heap64f::Heap64f *Q=NULL;
  gft::Image32::Image32 *pred;
  int i,p,q,n;
  double edge,tmp;
  double power = 9.0;
  double alpha = 0.03;
  int u_x,u_y,v_x,v_y;
  double *Dist;
  
  n    = img->ncols*img->nrows;
  pred = gft::Image32::Create(img);
  Q = gft::Heap64f::Create(n, cost);
  gft::Heap64f::SetRemovalPolicy(Q, MINVALUE);

  Dist = gft::AllocDoubleArray(A->n);
  for(i = 1; i < A->n; i++){
    Dist[i] = sqrt(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i]);
  }
  
  gft::Image32::Set(pred, NIL);
  for(p=0; p<n; p++){
    cost[p] = DBL_MAX;
  }
  
  cost[s] = 0.0;
  gft::Heap64f::Insert(Q, s);
  
  while(!gft::Heap64f::IsEmpty(Q)){
    gft::Heap64f::Remove(Q, &p);
    u_x = p%img->ncols; 
    u_y = p/img->ncols; 
    
    for(i=1; i<A->n; i++){
      v_x = u_x + A->dx[i];
      v_y = u_y + A->dy[i];
      if(v_x >= 0 && v_x < img->ncols &&
	 v_y >= 0 && v_y < img->nrows){
	q = v_x + img->ncols*v_y;
	if(Q->color[q] != BLACK){
	  
	  edge = img->data[p] + img->data[q];
	  
	  tmp = cost[p] + pow(alpha*edge, power) + Dist[i];
	  
	  if(tmp < cost[q]){
	    gft::Heap64f::Update(Q, q, tmp);
	    pred->data[q] = p;
	  }
	}
      }
    }
  }
  gft::Heap64f::Destroy(&Q);
  gft::FreeDoubleArray(&Dist);
  return pred;
}



gft::Image32::Image32 *Riverbed(gft::Image32::Image32 *img,
				gft::AdjRel::AdjRel *A,
				int s,
				double *cost){
  gft::Heap64f::Heap64f *Q=NULL;
  gft::Image32::Image32 *pred;
  int i,p,q,n,p1,p2;
  double edge,tmp;
  int u_x,u_y,v_x,v_y;
  
  n    = img->ncols*img->nrows;
  pred = gft::Image32::Create(img);
  Q = gft::Heap64f::Create(n, cost);
  gft::Heap64f::SetRemovalPolicy(Q, MINVALUE);

  gft::Image32::Set(pred, NIL);
  for(p=0; p<n; p++){
    cost[p] = DBL_MAX;
  }
  
  cost[s] = 0.0;
  gft::Heap64f::Insert(Q, s);
  
  while(!gft::Heap64f::IsEmpty(Q)){
    gft::Heap64f::Remove(Q, &p);
    u_x = p%img->ncols; 
    u_y = p/img->ncols; 
    
    for(i=1; i<A->n; i++){
      v_x = u_x + A->dx[i];
      v_y = u_y + A->dy[i];
      if(v_x >= 0 && v_x < img->ncols &&
	 v_y >= 0 && v_y < img->nrows){
	q = v_x + img->ncols*v_y;
	if(Q->color[q] != BLACK){
	  
	  edge = img->data[p] + img->data[q];
	  
	  tmp = edge;
	  
	  if(tmp < cost[q]){
	    gft::Heap64f::Update(Q, q, tmp);
	    pred->data[q] = p;
	  }
	}
      }
    }
  }
  gft::Heap64f::Destroy(&Q);
  return pred;
}



void Preproc(gft::Image32::Image32 *img){
  gft::Image32::Image32 *dil,*tmp;
  int *hist;
  int Imin,Imax,p,sum,lmin,lmax,lmax_p,val;
  gft::AdjRel::AdjRel *A;
  Imin = gft::Image32::GetMinVal(img);
  Imax = gft::Image32::GetMaxVal(img);
  //printf("Imin: %d\n", Imin);
  //printf("Imax: %d\n", Imax);

  A = gft::AdjRel::Circular(5.0);
  gft::AdjRel::Mult(A, 3);
  dil = gft::Image32::Dilate(img, A);
  tmp = gft::Image32::GaussianBlur(dil);
  gft::Image32::Destroy(&dil);
  dil = tmp;
  gft::Image32::Write(dil, (char *)"paper.pgm");
  
  hist = (int *)calloc((Imax-Imin+1), sizeof(int));
  for(p = 0; p < img->n; p++){
    hist[img->data[p]-Imin]++;
  }
  sum = 0;
  lmin = Imin;
  while(sum < 0.005*img->n){
    sum += hist[lmin-Imin];
    lmin++;
  }
  lmin = MAX(lmin-1, Imin);

  sum = 0;
  lmax = Imax;
  while(sum < 0.005*img->n){
    sum += hist[lmax-Imin];
    lmax--;
  }
  lmax = MIN(lmax+1, Imax);

  //printf("lmin: %d\n", lmin);
  //printf("lmax: %d\n", lmax);
  
  for(p = 0; p < img->n; p++){
    lmax_p = MIN(lmax, dil->data[p]);
    val = img->data[p];
    if(val <= lmin)        img->data[p] = 0;
    else if(val >= lmax_p) img->data[p] = 255;
    else
      img->data[p] = ROUND(((float)(img->data[p] - lmin)/(float)(lmax_p-lmin))*255.0);
  }
  free(hist);
  gft::Image32::Destroy(&dil);
  gft::AdjRel::Destroy(&A);
}


void ReadConfigFile(char *filename, gft::Image32::Image32 *img,
		    int S1[], int S2[], int GT[]){
  FILE *fp;
  int i,x,y,gt,p;
  fp = fopen(filename, "r");
  GT[0] = S1[0] = S2[0] = 0;
  if(fp != NULL){
    fscanf(fp, "%d", &S1[0]);
    GT[0] = S1[0];
    for(i = 1; i <= S1[0]; i++){
      fscanf(fp, "%d", &x);
      fscanf(fp, "%d", &y);
      fscanf(fp, "%d", &gt);
      if(!gft::Image32::IsValidPixel(img, x, y))
	gft::Error((char *)"Invalid file", (char *)"ReadConfigFile");
      p = x + y*img->ncols;
      S1[i] = p;
      GT[i] = gt;
    }

    fscanf(fp, "%d", &S2[0]);
    for(i = 1; i <= S2[0]; i++){
      fscanf(fp, "%d", &x);
      fscanf(fp, "%d", &y);
      if(!gft::Image32::IsValidPixel(img, x, y))
	gft::Error((char *)"Invalid file", (char *)"ReadConfigFile");
      p = x + y*img->ncols;
      S2[i] = p;
    }
    
    fclose(fp);
  }
}


void PrintPath(gft::CImage::CImage *paths,
	       int linecolor, int linespacing, float lineradius,
	       Path *P){
  gft::AdjRel::AdjRel *A = gft::AdjRel::Circular(lineradius);
  int i,q;
  for(i = 0; i < P->nnodes; i+= linespacing){
    q = P->nodes[i];
    gft::CImage::DrawAdjRel(paths, A, q, linecolor);
  }
  gft::AdjRel::Destroy(&A);
}



float ComputeOverlapRatio(gft::Image32::Image32 *img,
			  Path *P1, Path *P2){
  gft::Image32::Image32 *tmp;
  gft::AdjRel::AdjRel *A;
  A = gft::AdjRel::Circular(1.0);
  tmp = gft::Image32::Create(img);
  int i,p,overlap;
  overlap = 0;

  for(i = 0; i < P1->nnodes; i++){
    p = P1->nodes[i];
    gft::Image32::DrawAdjRel(tmp, A, p, 1);
  }

  for(i = 0; i < P2->nnodes; i++){
    p = P2->nodes[i];
    if(tmp->data[p] > 0)
      overlap++;
  }

  gft::AdjRel::Destroy(&A);
  gft::Image32::Destroy(&tmp);

  return (float)overlap/((P1->nnodes+P2->nnodes)/2.0);
}



int main(int argc, char **argv){
  gft::Image32::Image32 *img, *pred;
  gft::CImage::CImage *paths;
  gft::AdjRel::AdjRel *A;
  double *cost;
  char file[512];
  int S1[1024];
  int S2[1024];
  int GT[1024];
  int OUT[1024];
  int S[2];
  int i,j,k,n,correct,nt;
  int   LineColor[]   = {9, 0xFF0000, 0x00FF00, 0x0000FF, 0xFFFF00, 0x00FFFF, 0xFF00FF, 0x008888, 0x880088, 0x888800};
  int   LineSpacing[] = {9, 1, 1, 1, 1, 1, 1, 10, 10, 10};
  float LineRadius[]  = {9, 1, 1, 1, 1, 1, 1, 5, 5, 5};
  float ratio;
  GWireData *GWD;
  Path *P[100];
  int s,method;
  
  if(argc < 4){
    fprintf(stdout,"usage:\n");
    fprintf(stdout,"proj1 <image_file> <config_file> <method>\n");
    fprintf(stdout,"method: 0 ... Riverbed\n");
    fprintf(stdout,"        1 ... LiveWire\n");
    fprintf(stdout,"        2 ... G-wire\n");
    fprintf(stdout,"        3 ... Bandeirantes\n");
    exit(0);
  }

  A   = gft::AdjRel::Circular(1.5);
  img = gft::Image32::Read(argv[1]);
  ReadConfigFile(argv[2], img, S1, S2, GT);
  Preproc(img);

  method = atoi(argv[3]);
  
  n = img->ncols*img->nrows;
 
  Matrix<double> matrix(S1[0], S2[0]);
  Matrix<double> matrix2(S1[0], S2[0]);

  
  for(i = 1; i <= S1[0]; i++){
    s = S1[i];

    if(method != METHOD_GWIRE){
      cost = gft::AllocDoubleArray(n);

      if(method == METHOD_RIVERBED)
	pred = Riverbed(img, A, s, cost);
      else if(method == METHOD_LIVEWIRE)
	pred = LiveWire(img, A, s, cost);
      else if(method == METHOD_BANDEIRANTES)
	pred = Bandeirantes(img, A, s, cost);
      
      P[i-1] = GetIFTPaths(pred, cost, S2);
      gft::FreeDoubleArray(&cost);
      gft::Image32::Destroy(&pred);
    }
    else{
      GWD = GWire(img, s);
      P[i-1] = GetGWirePaths(GWD, S2);
      DestroyGWireData(&GWD);
    }
    
    for(j = 1; j <= S2[0]; j++){
      nt = MIN(15, P[i-1][j-1].nnodes-1);
      ComputePathCostByBandeirantes(&P[i-1][j-1], img);
      matrix(i-1, j-1) = P[i-1][j-1].cost[ nt ];
      matrix2(i-1, j-1) = matrix(i-1, j-1);
    }
    paths = ViewPaths(img, P[i-1], S2[0]);
    sprintf(file, "paths%02d.ppm", i);
    gft::CImage::Write(paths, file);

    gft::CImage::Destroy(&paths);
  }
  
  
  Munkres<double> m;
  
  m.solve(matrix);

  paths = gft::CImage::Clone(img);

  correct = 0;
  for(i = 1; i <= S1[0]; i++){
    for(j = 1; j <= S2[0]; j++){
      if(matrix(i-1, j-1) == 0){
	PrintPath(paths, LineColor[i], LineSpacing[i], LineRadius[i],
		  &P[i-1][j-1]);
	OUT[i] = j;
	if(GT[i] == j)
	  correct++;
      }
    }
  }

  printf("Score1: %.2f%%\n", 100.0*((float)correct/(float)S1[0]));

  gft::CImage::Write(paths, (char *)"best1.ppm");
  gft::CImage::Destroy(&paths);
  
  //printf("Ratio:\n");
  for(i = 1; i <= S1[0]; i++){
    for(j = i+1; j <= S1[0]; j++){
      ratio = ComputeOverlapRatio(img,
				  &P[i-1][OUT[i]-1],
				  &P[j-1][OUT[j]-1]);
      if(ratio > 0.10){
	//printf("Overlap detected\n");
	if(matrix2(i-1, OUT[i]-1) > matrix2(j-1, OUT[j]-1)){
	  for(k = 1; k <= S2[0]; k++){
	    ratio = ComputeOverlapRatio(img,
					&P[i-1][k-1],
					&P[j-1][OUT[j]-1]);
	    if(ratio > 0.10)
	      matrix2(i-1, k-1) = DBL_MAX;
	  }
	}
	else{
	  for(k = 1; k <= S2[0]; k++){
	    ratio = ComputeOverlapRatio(img,
					&P[j-1][k-1],
					&P[i-1][OUT[i]-1]);
	    if(ratio > 0.10)
	      matrix2(j-1, k-1) = DBL_MAX;
	  }
	}
      }
      //printf("%f ", ratio);
    }
    //printf("\n");
  }
  
  Munkres<double> m2;
  
  m2.solve(matrix2);

  paths = gft::CImage::Clone(img);

  correct = 0;
  for(i = 1; i <= S1[0]; i++){
    for(j = 1; j <= S2[0]; j++){
      if(matrix2(i-1, j-1) == 0){
	PrintPath(paths, LineColor[i], LineSpacing[i], LineRadius[i],
		  &P[i-1][j-1]);
	OUT[i] = j;
	if(GT[i] == j)
	  correct++;
      }
    }
  }
  
  printf("Score2: %.2f%%\n", 100.0*((float)correct/(float)S1[0]));

  gft::CImage::Write(paths, (char *)"best2.ppm");
  gft::CImage::Destroy(&paths);
  
  for(i = 1; i <= S1[0]; i++){
    for(j = 1; j <= S2[0]; j++){
      gft::FreeIntArray(&(P[i-1][j-1].nodes));
      gft::FreeDoubleArray(&(P[i-1][j-1].cost));
    }
    free(P[i-1]);
  }
  
  gft::Image32::Destroy(&img);
  gft::AdjRel::Destroy(&A);
  
  return 0;
}

