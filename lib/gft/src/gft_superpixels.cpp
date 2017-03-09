
#include "gft_superpixels.h"


namespace gft{
  namespace Superpixels{

    float ColorDistance(CImage::CImage *cimg_lab, 
			SPixelSeed seed, int q){
      float dc;
      dc = sqrtf(SQUARE(cimg_lab->C[0]->data[q] - seed.actual.l) + 
		 SQUARE(cimg_lab->C[1]->data[q] - seed.actual.a) +
		 SQUARE(cimg_lab->C[2]->data[q] - seed.actual.b));
      return dc;
    }

    float LastColorDistance(SPixelSeed seed){
      float dc;
      dc = sqrtf(SQUARE(seed.last_computed.l - seed.actual.l) + 
		 SQUARE(seed.last_computed.a - seed.actual.a) +
		 SQUARE(seed.last_computed.b - seed.actual.b));
      return dc;
    }

    float LastSpatialDistance(SPixelSeed seed){
      float ds;
      ds = sqrtf(SQUARE(seed.last_computed.i - seed.actual.i) + 
		 SQUARE(seed.last_computed.j - seed.actual.j));
      return ds;
    }
    
    void RunSPixelsByIFT(Image32::Image32 *grad,
			 CImage::CImage *cimg_lab,
			 Image32::Image32 *label,
			 Heap::Heap *Q,
			 Image32::Image32 *pred,
			 float *cost,
			 float alpha,
			 float beta,
			 SPixelSeed *seeds, int nseeds){
      //Image32::Image32 *root;
      AdjRel::AdjRel *A;
      float *dpq = NULL;
      int i, p,q, s;
      int px,py,qx,qy;
      float edge,wl,wg,w,tmp;
      float Cmax, Gmax;
      
      Cmax = sqrtf(100*100 + 256*256 + 256*256);
      Gmax = Image32::GetMaxVal(grad);
      
      A = AdjRel::Neighborhood_4();
      dpq = (float *)AllocFloatArray(A->n);
      for(i = 1; i < A->n; i++){
	dpq[i] = sqrtf(A->dx[i]*A->dx[i] + 
		       A->dy[i]*A->dy[i]);
      }
      
      //root = Image32::Create(grad->ncols, grad->nrows);
      
      for(s = 0; s < nseeds; s++){
	if(seeds[s].recompute){
	  p = seeds[s].actual.j + grad->ncols*seeds[s].actual.i;
	  label->data[p] = s;
	  pred->data[p] = NIL;
	  //Heap::Insert(Q, p);
	  cost[p] = 0.0;
	  Heap::Update(Q, p, 0.0);
	}
      }
      
      while(!Heap::IsEmpty(Q)){
	Heap::Remove(Q, &p);
	px = p%grad->ncols;
	py = p/grad->ncols;
	s = label->data[p];
	
	for(i = 1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  if(Image32::IsValidPixel(grad, qx, qy)){
	    q = qx + qy*grad->ncols;
	    if(Q->color[q] != BLACK){

	      //wl = ColorDistance(cimg_lab, root->data[p], q);
	      wl = ColorDistance(cimg_lab, seeds[s], q);
	      //wg = ((grad->data[p] + grad->data[q])/2.0)/(Gmax/2.0);
	      edge = alpha*wl; //alpha*wg;
	      w = powf(edge, beta) + dpq[i];
	      tmp = cost[p] + w;
	      if(tmp < cost[q] || pred->data[q] == p){
		Heap::Update(Q, q, tmp);
		//Q->cost[q] = tmp;
		//if(Q->color[q] == GRAY)
		//  GoUp(Q, Q->pos[q]);
		//else
		//  Insert(Q, q);

		//root->data[q] = root->data[p];
		pred->data[q] = p;
		label->data[q] = label->data[p];
	      }
	    }
	  }
	}
      }

      for(s = 0; s < nseeds; s++){
	if(seeds[s].recompute){
	  seeds[s].last_computed = seeds[s].actual;
	  seeds[s].recompute = false;
	}
      }
      //Image32::Destroy(&root);
      FreeFloatArray(&dpq);
      AdjRel::Destroy(&A);
    }
    
    
    float WeightedDistanceMeasure(CImage::CImage *cimg_lab, 
				  SPixelSeed seed, int i, int j,
				  float m, float S){
      Color::ColorLab lab;
      float dc,ds,D;
      lab.l = cimg_lab->C[0]->array[i][j];
      lab.a = cimg_lab->C[1]->array[i][j];
      lab.b = cimg_lab->C[2]->array[i][j];
      dc = sqrtf(SQUARE(seed.actual.l - lab.l) + 
		 SQUARE(seed.actual.a - lab.a) +
		 SQUARE(seed.actual.b - lab.b));
      ds = sqrtf(SQUARE(seed.actual.j - j) +
		 SQUARE(seed.actual.i - i));
      D = sqrtf(dc*dc + (ds/S)*(ds/S)*m*m);  
      //Na verdade, nao precisa do sqrtf pois nao vai mudar a ordem entre diferentes distancias.
      return D;
    }
    

    void Move2LowerGradient(Image32::Image32 *grad, 
			    int *j, int *i){
      AdjRel::AdjRel *A;
      int u_x, u_y;
      int v_x, v_y;
      int k, g, Gmin;
      u_x = *j;
      u_y = *i;
      Gmin = grad->array[u_y][u_x];
      A = AdjRel::Circular(1.5);
      for(k = 0; k < A->n; k++){
	v_x = u_x + A->dx[k];
	v_y = u_y + A->dy[k];
	if(Image32::IsValidPixel(grad, v_x, v_y)){
	  g = grad->array[v_y][v_x];
	  if(g < Gmin){
	    Gmin = g;
	    *j = v_x;
	    *i = v_y;
	  }
	}
      }
      AdjRel::Destroy(&A);
    }
    
    
    void Postprocessing_CC(Image32::Image32 *label,
			   SPixelSeed *seeds, int nseeds){
      Image32::Image32 *CC;
      AdjRel::AdjRel *A;
      int p,pp,q, n, i, lb = 0;
      int px,py,ppx,ppy,qx,qy;
      int adjlabel = 0;
      int SUPSZ;
      int *Q;
      int inic,fim;
      
      n = label->n;
      SUPSZ = n/nseeds;
      A = AdjRel::Neighborhood_4();
      A->dx[1] = -1; A->dy[1] = 0;  /* left */
      A->dx[2] = 0;  A->dy[2] = -1; /* top */
      A->dx[3] = 1;  A->dy[3] = 0;  /* right */
      A->dx[4] = 0;  A->dy[4] = 1;  /* bottom */
      CC = Image32::Create(label->ncols, label->nrows);
      Q = (int *)malloc(sizeof(int)*n);
      Image32::Set(CC, NIL);
      for(p = 0; p < n; p++){
	if(CC->data[p] != NIL) continue;
	
	CC->data[p] = lb;
	px = p%label->ncols;
	py = p/label->ncols;
	for(i = 1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  if(IsValidPixel(label, qx, qy)){
	    q = qx + qy*label->ncols;
	    if(CC->data[q] != NIL)
	      adjlabel = CC->data[q];
	  }
	}
	
	inic = fim = 0;
	Q[fim] = p;
	fim++;
	
	while(inic < fim){
	  pp = Q[inic];
	  inic++;
	  ppx = pp%label->ncols;
	  ppy = pp/label->ncols;
	  for(i = 1; i < A->n; i++){
	    qx = ppx + A->dx[i];
	    qy = ppy + A->dy[i];
	    if(IsValidPixel(label, qx, qy)){
	      q = qx + qy*label->ncols;
	      if(CC->data[q] == NIL &&  
		 label->data[pp] == label->data[q]){
		Q[fim] = q;
		fim++;
		CC->data[q] = lb;
	      }
	    }
	  }
	}
	//----------------------------------------------------
	// If segment size is less then a limit, assign an
	// adjacent label found before, and decrement label count.
	//----------------------------------------------------
	if(fim <= SUPSZ/4){
	  for(inic = 0; inic < fim; inic++){
	    pp = Q[inic];
	    CC->data[pp] = adjlabel;
	  }
	  lb--;
	}
	lb++;
      }
      
      for(p = 0; p < n; p++){
	label->data[p] = CC->data[p];
      }
      
      Image32::Destroy(&CC);  
      AdjRel::Destroy(&A);
      free(Q);
    }
    

    void Postprocessing(Image32::Image32 *label,
			SPixelSeed *seeds, int nseeds){
      Image32::Image32 *tmp, *orphans;
      Queue::Queue *Q;
      AdjRel::AdjRel *A;
      int s,p,q,k;
      int p_x,p_y,q_x,q_y;
      int orphaned = 0;
      A = AdjRel::Neighborhood_4();
      Q = Queue::Create(label->n);
      tmp = Image32::Create(label->ncols, label->nrows);
      orphans = Image32::Create(label->ncols, label->nrows);
      Image32::Set(tmp, NIL);
      
      /*
	for(p = 0; p < label->n; p++){
	if(label->data[p] == NIL)
	printf("NIL\n");
	}
      */
      
      for(s = 0; s < nseeds; s++){
	p = seeds[s].last_computed.j + label->ncols*seeds[s].last_computed.i;
	if(label->data[p] != s)
	  printf("Seed %d=(%d, %d) outside its superpixel\n", 
		 s, seeds[s].last_computed.j, seeds[s].last_computed.i);
	tmp->data[p] = s;
	Queue::Push(Q, p);
      }
      
      while(!Queue::IsEmpty(Q)){
	p = Queue::Pop(Q);
	p_x = p%label->ncols;
	p_y = p/label->ncols;
	for(k = 1; k < A->n; k++){
	  q_x = p_x + A->dx[k];
	  q_y = p_y + A->dy[k];
	  if(Image32::IsValidPixel(label, q_x, q_y)){
	    q = q_x + q_y*label->ncols;
	    if(label->data[p] == label->data[q] && tmp->data[q] == NIL){
	      tmp->data[q] = label->data[p];
	      Queue::Push(Q, q);
	    }
	  }
	}
      }
      
      for(p = 0; p < label->n; p++){
	if(tmp->data[p] == NIL){
	  orphans->data[p] = 255;
	  orphaned++;
	}
      }
      
      Image32::Write(orphans, (char *)"orphans.pgm");
      printf("orphaned: %d\n", orphaned);
      
      gft::Queue::Reset(Q);
      
      for(p = 0; p < label->n; p++){
	if(tmp->data[p] == NIL) continue;
	p_x = p%label->ncols;
	p_y = p/label->ncols;
	for(k = 1; k < A->n; k++){
	  q_x = p_x + A->dx[k];
	  q_y = p_y + A->dy[k];
	  if(Image32::IsValidPixel(label, q_x, q_y)){
	    if(tmp->array[q_y][q_x] == NIL){
	      Queue::Push(Q, p);
	      break;
	    }
	  }
	}
      }
      
      while(!Queue::IsEmpty(Q)){
	p = Queue::Pop(Q);
	p_x = p%label->ncols;
	p_y = p/label->ncols;
	for(k = 1; k < A->n; k++){
	  q_x = p_x + A->dx[k];
	  q_y = p_y + A->dy[k];
	  if(Image32::IsValidPixel(label, q_x, q_y)){
	    q = q_x + q_y*label->ncols;
	    if(tmp->data[q] == NIL){
	      tmp->data[q] = tmp->data[p];
	      Queue::Push(Q, q);
	    }
	  }
	}
      }
      
      for(p = 0; p < label->n; p++){
	label->data[p] = tmp->data[p];
      }
      
      AdjRel::Destroy(&A);
      Image32::Destroy(&tmp);
      Image32::Destroy(&orphans);
      Queue::Destroy(&Q);
    }
    
    
    SPixelSeed *GetSPixelsSeeds(CImage::CImage *cimg_lab,
				Image32::Image32 *grad,
				int k,
				int *nseeds){
      Image32::Image32 *label;
      int N,ncols,nrows,i,j;
      int S;
      SPixelSeed *seeds;
      int s;
      int xstrips, ystrips;
      int xerr, yerr;
      double xerrperstrip, yerrperstrip;
      int xoff, yoff;
      int ye, xe;
      int perturbseeds = false;
      
      *nseeds = 0;

      ncols = cimg_lab->C[0]->ncols;
      nrows = cimg_lab->C[0]->nrows;
      N = cimg_lab->C[0]->n;
      S = ROUND( sqrtf((float)N/(float)k) );
      
      xstrips = (0.5 + double(ncols)/double(S));
      ystrips = (0.5 + double(nrows)/double(S));
      
      xerr = ncols - S*xstrips; 
      if(xerr < 0){ xstrips--; xerr = ncols - S*xstrips;}
      yerr = nrows - S*ystrips;
      if(yerr < 0){ ystrips--; yerr = nrows - S*ystrips;}
      
      xerrperstrip = double(xerr)/double(xstrips);
      yerrperstrip = double(yerr)/double(ystrips);
      
      label = Image32::Create(ncols, nrows);
      
      xoff = S/2;
      yoff = S/2;
      
      *nseeds = xstrips*ystrips;
      
      seeds = (SPixelSeed *)malloc((*nseeds)*sizeof(SPixelSeed));
      
      s = 0;
      for(i = 0; i < ystrips; i++){
	ye = i*yerrperstrip;
	for(j = 0; j < xstrips; j++){
	  xe = j*xerrperstrip;
	  seeds[s].actual.j = (j*S + xoff+xe);
	  seeds[s].actual.i = (i*S + yoff+ye);
	  
	  //printf("Seeds: (%d,%d)\n", seeds[s].actual.j, seeds[s].actual.i);
	  
	  if(Image32::IsValidPixel(label, seeds[s].actual.j, seeds[s].actual.i)){
	    label->array[seeds[s].actual.i][seeds[s].actual.j] = 128;
	    Image32::DrawCircle(label, seeds[s].actual.j, seeds[s].actual.i, 4.5, 128);
	  }
	  s++;
	}
      }
      
      for(s = 0; s < *nseeds; s++){
	if(perturbseeds)
	  Move2LowerGradient(grad, &(seeds[s].actual.j), &(seeds[s].actual.i));
	i = seeds[s].actual.i;
	j = seeds[s].actual.j;
	//label->array[i][j] = 255;
	seeds[s].actual.l = cimg_lab->C[0]->array[i][j];
	seeds[s].actual.a = cimg_lab->C[1]->array[i][j];
	seeds[s].actual.b = cimg_lab->C[2]->array[i][j];
	seeds[s].recompute = true;
      }
      
      Image32::Write(label, (char *)"seeds.pgm");
      Image32::Destroy(&label);
      
      printf("Initial number of seeds: %d\n", *nseeds);
      return seeds;
    }


    Image32::Image32 *IFT_SLIC(CImage::CImage *cimg, bool colored,
			       int k, float alpha, float beta,
			       float err_dc, float err_ds, int itMax){
      CImage::CImage *cimg_lab;
      Image32::Image32 *label, *img, *grad;
      Heap::Heap *Q;
      Image32::Image32 *pred;
      float *cost;
      int N,ncols,nrows,S,p,q,i;
      SPixelSeed *seeds;
      int it, s, nseeds = 0;
      int p_x,p_y,q_x,q_y;
      AdjRel::AdjRel *A;
      //---------------------
      char filename[512];
      
      A = AdjRel::Neighborhood_4();      
      if(colored){
	cimg_lab = CImage::RGB2Lab(cimg);
	img = CImage::Luminosity(cimg);
      }
      else{
	cimg_lab = CImage::Clone(cimg);
        img = Image32::Clone(cimg->C[0]);
      }
      //Image32::Write(img, (char *)"luminosity.pgm");
      grad = Image32::SobelFilter(img);
      //Image32::Write(grad, (char *)"grad.pgm");
      
      ncols = cimg->C[0]->ncols;
      nrows = cimg->C[0]->nrows;
      N = cimg->C[0]->n;
      S = ROUND( sqrtf((float)N/(float)k) );
      
      seeds = GetSPixelsSeeds(cimg_lab, grad, k, &nseeds);
      
      label = Image32::Create(ncols, nrows);

      pred = Image32::Create(ncols, nrows);
      cost = AllocFloatArray(N);
      for(p = 0; p < N; p++){
	cost[p] = FLT_MAX;
	//root->data[p] = p;
      }
      Image32::Set(pred, NIL);
      Image32::Set(label, NIL);
      Q = Heap::Create(N, cost);
      Heap::SetRemovalPolicy(Q, MINVALUE);

      for(it = 0; it < itMax; it++){
	RunSPixelsByIFT(grad, cimg_lab, label, 
			Q, pred, cost,
			alpha, beta,
			seeds, nseeds);

	//---------------------
	//sprintf(filename, "label%02d.pgm", it);
	//Image32::Write(label, filename);
	//---------------------
	
	UpdateSPixelsSeeds(cimg_lab, label,
			   &seeds, nseeds, true); //false);
	if(it < itMax - 1){
	  for(s = 0; s < nseeds; s++){
	    if( LastColorDistance(seeds[s]) > err_dc ||
		LastSpatialDistance(seeds[s]) > err_ds ||
		seeds[s].n < (N/nseeds)/5.0 ){
	      seeds[s].recompute = true;
	    }
	  }
	  for(p = 0; p < N; p++){
	    Q->color[p] = WHITE;
	    s = label->data[p];
	    if(s == NIL) continue;
	    if(seeds[s].recompute){
	      cost[p] = FLT_MAX;
	      pred->data[p] = NIL;
	      label->data[p] = NIL;
	    }
	  }

	  //Percorre pixels:
	  //Para pixels com NIL, verifica vizinhos,
	  //Se tem vizinho != NIL, entao insere na fila esse vizinho
          //com o custo que ele tinha.
	  for(p = 0; p < N; p++){
	    if(label->data[p] != NIL) continue;
	    p_x = p%label->ncols;
	    p_y = p/label->ncols;
	    for(i = 1; i < A->n; i++){
	      q_x = p_x + A->dx[i];
	      q_y = p_y + A->dy[i];
	      if(Image32::IsValidPixel(label, q_x, q_y)){
		q = q_x + q_y*label->ncols;
		if(label->data[q] != NIL)
		  Heap::Update(Q, q, cost[q]);
	      }
	    }
	  }

	}
      }
      //gft::Image32::Write(label, (char *)"label_ift_before.pgm");
      //Postprocessing_CC(label, seeds, nseeds);

      Heap::Destroy(&Q);
      Image32::Destroy(&pred);
      FreeFloatArray(&cost);
      
      Image32::Destroy(&img);
      Image32::Destroy(&grad);
      CImage::Destroy(&cimg_lab);
      AdjRel::Destroy(&A);
      free(seeds);
      return label;
    }


    Image32::Image32 *IFT_SLIC(CImage::CImage *cimg,
			       int k, float alpha, float beta,
			       float err_dc, float err_ds, int itMax){
      return IFT_SLIC(cimg, true,
		      k, alpha, beta,
		      err_dc, err_ds, itMax);
    }


    Image32::Image32 *IFT_SLIC(Image32::Image32 *img,
			       int k, float alpha, float beta,
			       float err_dc, float err_ds, int itMax){
      CImage::CImage *cimg;
      Image32::Image32 *label;
      cimg = CImage::Clone(img);
      label = IFT_SLIC(cimg, false,
		       k, alpha/sqrtf(3.0), beta,
		       err_dc, err_ds, itMax);
      CImage::Destroy(&cimg);
      return label;
    }


    void UpdateSPixelsSeeds(CImage::CImage *cimg_lab,
			    Image32::Image32 *label,
			    SPixelSeed **seeds, int nseeds,
			    int inside){
      int s,p,N,i,j;
      int ncols,nrows;
      SPixelSeed *newseeds;
      //------------------------
      Image32::Image32 *tmp;
      static int t = 0;
      char filename[512];
      
      ncols = cimg_lab->C[0]->ncols;
      nrows = cimg_lab->C[0]->nrows;
      N = cimg_lab->C[0]->n;
      newseeds = (SPixelSeed *)malloc(nseeds*sizeof(SPixelSeed));
      for(s = 0; s < nseeds; s++){
	newseeds[s].n = 0;
	newseeds[s].actual.i = 0;
	newseeds[s].actual.j = 0;
	newseeds[s].actual.l = 0.0;
	newseeds[s].actual.a = 0.0;
	newseeds[s].actual.b = 0.0;
      }
	
      for(p = 0; p < N; p++){
	s = label->data[p];
	if(s == NIL) continue;
	i = p/ncols;
	j = p%ncols;
	
	newseeds[s].n++;
	newseeds[s].actual.i += i;
	newseeds[s].actual.j += j;
	newseeds[s].actual.l += cimg_lab->C[0]->array[i][j];
	newseeds[s].actual.a += cimg_lab->C[1]->array[i][j];
	newseeds[s].actual.b += cimg_lab->C[2]->array[i][j];
      }
      
      for(s = 0; s < nseeds; s++){
	if(newseeds[s].n > 0){
	  newseeds[s].actual.i /= newseeds[s].n;
	  newseeds[s].actual.j /= newseeds[s].n;
	  newseeds[s].actual.l /= newseeds[s].n;
	  newseeds[s].actual.a /= newseeds[s].n;
	  newseeds[s].actual.b /= newseeds[s].n;
	}
	else
	  printf("s: %d at (%d, %d)\n", 
		 s, (*seeds)[s].actual.j, (*seeds)[s].actual.i);
      }
	
      //------------------------------------------
      //To enforce the placement of seeds inside superpixels:
      if(inside){
	for(p = 0; p < N; p++){
	  s = label->data[p];
	  if(s == NIL) continue;
	  i = p/ncols;
	  j = p%ncols;
	  
	  if( SQUARE(j-newseeds[s].actual.j) + SQUARE(i-newseeds[s].actual.i)  < 
	      SQUARE((*seeds)[s].actual.j-newseeds[s].actual.j) + SQUARE((*seeds)[s].actual.i-newseeds[s].actual.i) ){
	    (*seeds)[s].actual.j = j;
	    (*seeds)[s].actual.i = i;
	  }
	}
	for(s = 0; s < nseeds; s++){
	  newseeds[s].actual.j = (*seeds)[s].actual.j;
	  newseeds[s].actual.i = (*seeds)[s].actual.i;
	}
      }

      for(s = 0; s < nseeds; s++){
	newseeds[s].last_computed = (*seeds)[s].last_computed;
	newseeds[s].recompute = (*seeds)[s].recompute;
      }
      free(*seeds);
      *seeds = newseeds;

      //--------------------
      tmp = Image32::Create(ncols, nrows);
      for(s = 0; s < nseeds; s++){
	j = newseeds[s].actual.j;
	i = newseeds[s].actual.i;
	p = j + i*ncols;
	tmp->data[p] = 128;
	Image32::DrawCircle(tmp, j, i, 4.5, 128);
      }
      t++;
      //sprintf(filename,"seeds%02d.pgm", t);
      //Image32::Write(tmp, filename);
      Image32::Destroy(&tmp);
    }


    Image32::Image32 *mySLIC(CImage::CImage *cimg, bool colored,
			     int k, float m){
      CImage::CImage *cimg_lab;
      Image32::Image32 *label, *img, *grad;
      float *dist, D;
      int N,ncols,nrows,p,i,j,si,sj;
      int S;
      SPixelSeed *seeds;
      int it, s, nseeds = 0;

      if(colored){      
	cimg_lab = CImage::RGB2Lab(cimg);
	img = CImage::Luminosity(cimg);
      }
      else{
	cimg_lab = CImage::Clone(cimg);
        img = Image32::Clone(cimg->C[0]);
      }
      //Image32::Write(img, (char *)"luminosity.pgm");
      grad = Image32::SobelFilter(img);
      //Image32::Write(grad, (char *)"grad.pgm");
      
      ncols = cimg->C[0]->ncols;
      nrows = cimg->C[0]->nrows;
      N = cimg->C[0]->n;
      S = ROUND( sqrtf((float)N/(float)k) );
      
      seeds = GetSPixelsSeeds(cimg_lab, grad, k, &nseeds);
      
      Image32::Destroy(&img);
      Image32::Destroy(&grad);
      
      label = Image32::Create(ncols, nrows);
      dist = AllocFloatArray(N);
      
      for(it = 0; it < 10; it++){
	
	Image32::Set(label, NIL);
	for(p = 0; p < N; p++)
	  dist[p] = FLT_MAX;
	
	for(s = 0; s < nseeds; s++){
	  si = seeds[s].actual.i;
	  sj = seeds[s].actual.j;
	  p = sj + si*ncols;
	  dist[p] = 0;
	  label->data[p] = s;
	  
	  for(i = ROUND(si - S); i <= ROUND(si + S); i++){
	    for(j = ROUND(sj - S); j <= ROUND(sj + S); j++){
	      if(Image32::IsValidPixel(label, j, i)){
		D = WeightedDistanceMeasure(cimg_lab, seeds[s], i, j, m, S);
		p = j + i*ncols;
		if(D < dist[p]){
		  dist[p] = D;
		  label->data[p] = s;
		}
	      }
	    }
	  }
	}

	UpdateSPixelsSeeds(cimg_lab, label,
			   &seeds, nseeds, false);

      }
      
      //Image32::Write(label, (char *)"label_before.pgm");
      Postprocessing_CC(label, seeds, nseeds);
      
      CImage::Destroy(&cimg_lab);
      FreeFloatArray(&dist);
      free(seeds);
      return label;
    }

    
    Image32::Image32 *mySLIC(CImage::CImage *cimg, 
			     int k, float m){
      return mySLIC(cimg, true, k, m);
    }
    
    
    Image32::Image32 *mySLIC(Image32::Image32 *img, 
			     int k, float m){
      CImage::CImage *cimg;
      Image32::Image32 *label;
      cimg = CImage::Clone(img);      
      label = mySLIC(cimg, false, k, m*sqrtf(3.0));
      CImage::Destroy(&cimg);
      return label;
    }
    
    
    int GetNumberOfSuperPixels(Image32::Image32 *label){
      int p;
      int Lmin,Lmax;
      Lmin = INT_MAX;
      Lmax = INT_MIN;
      for(p = 0; p < label->n; p++){
	if(label->data[p] < Lmin)
	  Lmin = label->data[p];
	if(label->data[p] > Lmax)
	  Lmax = label->data[p];
      }
      //printf("Lmin: %d, Lmax: %d\n", Lmin, Lmax);
      return Lmax-Lmin+1;
    }

    //-----------------------------------------------

    float FeatureDistance(Scene32::Scene32 *scn, 
			  SVoxelSeed seed, int q){
      float df;
      df = fabsf(scn->data[q] - seed.actual.l);
      return df;
    }

    float LastFeatureDistance(SVoxelSeed seed){
      float dc;
      dc = fabsf(seed.last_computed.l - seed.actual.l);
      return dc;
    }

    float LastSpatialDistance(SVoxelSeed seed){
      float ds;
      ds = sqrtf(SQUARE(seed.last_computed.i - seed.actual.i) + 
		 SQUARE(seed.last_computed.j - seed.actual.j) +
		 SQUARE(seed.last_computed.k - seed.actual.k));
      return ds;
    }


    void RunSVoxelsByIFT(Scene32::Scene32 *scn,
			 Scene32::Scene32 *label,
			 Heap::Heap *Q,
			 Scene32::Scene32 *pred,
			 float *cost,
			 float alpha,
			 float beta,
			 SVoxelSeed *seeds, int nseeds){
      AdjRel3::AdjRel3 *A;
      float *dpq = NULL;
      int i, p,q, s;
      gft::Voxel u,v;
      float edge,wl,w,tmp;
      
      A = AdjRel3::Spheric(1.0);
      dpq = (float *)AllocFloatArray(A->n);
      for(i = 1; i < A->n; i++){
	dpq[i] = sqrtf(A->d[i].axis.x*A->d[i].axis.x + 
		       A->d[i].axis.y*A->d[i].axis.y +
		       A->d[i].axis.z*A->d[i].axis.z);
      }
      
      for(s = 0; s < nseeds; s++){
	if(seeds[s].recompute){
	  p = Scene32::GetVoxelAddress(scn,
				       seeds[s].actual.j,
				       seeds[s].actual.i,
				       seeds[s].actual.k);
	  label->data[p] = s;
	  pred->data[p] = NIL;
	  cost[p] = 0.0;
	  Heap::Update(Q, p, 0.0);
	}
      }
      
      while(!Heap::IsEmpty(Q)){
	Heap::Remove(Q, &p);
	u.c.x = Scene32::GetAddressX(scn, p);
	u.c.y = Scene32::GetAddressY(scn, p);
	u.c.z = Scene32::GetAddressZ(scn, p);
	s = label->data[p];
	
	for(i = 1; i < A->n; i++){
	  //v.c.x = u.c.x + A->d[i].axis.x;
	  //v.c.y = u.c.y + A->d[i].axis.y;
	  //v.c.z = u.c.z + A->d[i].axis.z;
	  v.v = u.v + A->d[i].v;
	  if(Scene32::IsValidVoxel(scn, v)){
	    q = Scene32::GetVoxelAddress(scn, v);

	    if(Q->color[q] != BLACK){
	      wl = FeatureDistance(scn, seeds[s], q);
	      edge = alpha*wl; 
	      w = powf(edge, beta) + dpq[i];
	      tmp = cost[p] + w;
	      if(tmp < cost[q] || pred->data[q] == p){
		Heap::Update(Q, q, tmp);
		pred->data[q] = p;
		label->data[q] = label->data[p];
	      }
	    }
	  }
	}
      }
      
      for(s = 0; s < nseeds; s++){
	if(seeds[s].recompute){
	  seeds[s].last_computed = seeds[s].actual;
	  seeds[s].recompute = false;
	}
      }
      FreeFloatArray(&dpq);
      AdjRel3::Destroy(&A);
    }
    

    Scene32::Scene32 *IFT_SLIC(Scene32::Scene32 *scn,
			       int k, float alpha, float beta,
			       float err_df, float err_ds, int itMax){
      Scene32::Scene32 *label, *pred;
      Heap::Heap *Q;
      float *cost;
      int N,S,p,q,i;
      SVoxelSeed *seeds;
      int it, s, nseeds = 0;
      Voxel u,v;
      AdjRel3::AdjRel3 *A;
      //------------------
      FILE *fp;
      
      A = AdjRel3::Spheric(1.0);
      N = scn->n;
      S = ROUND( sqrtf((float)N/(float)k) );

      seeds = GetSVoxelsSeeds(scn, k, &nseeds);
      
      label = Scene32::Create(scn);
      pred  = Scene32::Create(scn);
      cost = AllocFloatArray(N);
      for(p = 0; p < N; p++){
	cost[p] = FLT_MAX;
      }
      Scene32::Fill(pred,  NIL);
      Scene32::Fill(label, NIL);
      Q = Heap::Create(N, cost);
      Heap::SetRemovalPolicy(Q, MINVALUE);

      for(it = 0; it < itMax; it++){
	RunSVoxelsByIFT(scn, label,
			Q, pred, cost,
			alpha, beta,
			seeds, nseeds);

	UpdateSVoxelsSeeds(scn, label,
			   &seeds, nseeds, true); //false);
	
	if(it < itMax - 1){
	  for(s = 0; s < nseeds; s++){
	    if( LastFeatureDistance(seeds[s]) > err_df ||
		LastSpatialDistance(seeds[s]) > err_ds ||
		seeds[s].n < (N/nseeds)/5.0 ){
	      seeds[s].recompute = true;
	    }
	  }
	  for(p = 0; p < N; p++){
	    Q->color[p] = WHITE;
	    s = label->data[p];
	    if(s == NIL) continue;
	    if(seeds[s].recompute){
	      cost[p] = FLT_MAX;
	      pred->data[p] = NIL;
	      label->data[p] = NIL;
	    }
	  }

	  //Percorre pixels:
	  //Para pixels com NIL, verifica vizinhos,
	  //Se tem vizinho != NIL, entao insere na fila esse vizinho
          //com o custo que ele tinha.
	  for(p = 0; p < N; p++){
	    if(label->data[p] != NIL) continue;
	    u.c.x = Scene32::GetAddressX(scn, p);
	    u.c.y = Scene32::GetAddressY(scn, p);
	    u.c.z = Scene32::GetAddressZ(scn, p);

	    for(i = 1; i < A->n; i++){
	      v.c.x = u.c.x + A->d[i].axis.x;
	      v.c.y = u.c.y + A->d[i].axis.y;
	      v.c.z = u.c.z + A->d[i].axis.z;

	      if(Scene32::IsValidVoxel(scn, v)){
		q = Scene32::GetVoxelAddress(scn, v);

		if(label->data[q] != NIL)
		  Heap::Update(Q, q, cost[q]);
	      }
	    }
	  }
	}
      }

      Heap::Destroy(&Q);
      Scene32::Destroy(&pred);
      FreeFloatArray(&cost);
      AdjRel3::Destroy(&A);
      //--------------------      
      fp = fopen("SP_seeds_report.txt", "w");
      if(fp != NULL){
	fprintf(fp, "%d\n", nseeds);
	for(i = 0; i < nseeds; i++){
	  fprintf(fp, "%d %d %d %d %f\n",
		  seeds[i].actual.j,
		  seeds[i].actual.i,
		  seeds[i].actual.k,
		  i,
		  seeds[i].actual.l);
	}
	fclose(fp);
      }
      //--------------------      
      free(seeds);
      return label;
    }


    
    SVoxelSeed *GetSVoxelsSeeds(Scene32::Scene32 *scn,
				int k,
				int *nseeds){
      Scene32::Scene32 *label;
      int N,i,j,l,S;
      SVoxelSeed *seeds;
      int s;
      int xstrips, ystrips, zstrips;
      int xerr, yerr, zerr;
      double xerrperstrip, yerrperstrip, zerrperstrip;
      int xoff, yoff, zoff;
      int ye, xe, ze;
      
      *nseeds = 0;

      N = scn->n;
      if(scn->zsize > 1)
	S = ROUND( pow((double)N/(double)k, 1.0/3.0) );
      else
	S = ROUND( pow((double)N/(double)k, 1.0/2.0) );
      
      xstrips = (0.5 + double(scn->xsize)/double(S));
      ystrips = (0.5 + double(scn->ysize)/double(S));
      zstrips = (0.5 + double(scn->zsize)/double(S));
      
      xerr = scn->xsize - S*xstrips; 
      if(xerr < 0){ xstrips--; xerr = scn->xsize - S*xstrips;}

      yerr = scn->ysize - S*ystrips;
      if(yerr < 0){ ystrips--; yerr = scn->ysize - S*ystrips;}

      zerr = scn->zsize - S*zstrips;
      if(zerr < 0){ zstrips--; zerr = scn->zsize - S*zstrips;}
      
      xerrperstrip = double(xerr)/double(xstrips);
      yerrperstrip = double(yerr)/double(ystrips);
      if(zstrips != 0.0)
	zerrperstrip = double(zerr)/double(zstrips);
      else
	zerrperstrip = zerr;
      
      label = Scene32::Create(scn);
      
      xoff = S/2;
      yoff = S/2;
      if(scn->zsize > 1)
	zoff = S/2;
      else
	zoff = 0;

      if(scn->zsize > 1)
	*nseeds = xstrips*ystrips*zstrips;
      else{
	zstrips = 1;
	*nseeds = xstrips*ystrips*zstrips;
      }
	
      seeds = (SVoxelSeed *)malloc((*nseeds)*sizeof(SVoxelSeed));
      
      s = 0;
      for(l = 0; l < zstrips; l++){
	ze = l*zerrperstrip;
	for(i = 0; i < ystrips; i++){
	  ye = i*yerrperstrip;
	  for(j = 0; j < xstrips; j++){
	    xe = j*xerrperstrip;
	  
	    seeds[s].actual.j = (j*S + xoff+xe);
	    seeds[s].actual.i = (i*S + yoff+ye);
	    seeds[s].actual.k = (l*S + zoff+ze);
	  
	    if(Scene32::IsValidVoxel(label,
				     seeds[s].actual.j,
				     seeds[s].actual.i,
				     seeds[s].actual.k)){
	      label->array[seeds[s].actual.k][seeds[s].actual.i][seeds[s].actual.j] = 128;
	      //Image32::DrawCircle(label, seeds[s].actual.j, seeds[s].actual.i, 4.5, 128);
	    }
	    s++;
	  }
	}
      }
      
      for(s = 0; s < *nseeds; s++){
	i = seeds[s].actual.i;
	j = seeds[s].actual.j;
	l = seeds[s].actual.k;
	seeds[s].actual.l = scn->array[l][i][j];
	seeds[s].recompute = true;
      }
      
      Scene32::Write(label, (char *)"seeds.scn");
      Scene32::Destroy(&label);
      
      printf("Initial number of seeds: %d\n", *nseeds);
      return seeds;
    }


    void UpdateSVoxelsSeeds(Scene32::Scene32 *scn,
			    Scene32::Scene32 *label,
			    SVoxelSeed **seeds, int nseeds,
			    int inside){
      int s,p,N,i,j,k;
      SVoxelSeed *newseeds;
      
      N = scn->n;
      newseeds = (SVoxelSeed *)malloc(nseeds*sizeof(SVoxelSeed));
      for(s = 0; s < nseeds; s++){
	newseeds[s].n = 0;
	newseeds[s].actual.i = 0;
	newseeds[s].actual.j = 0;
	newseeds[s].actual.k = 0;
	newseeds[s].actual.l = 0.0;
      }
	
      for(p = 0; p < N; p++){
	s = label->data[p];
	if(s == NIL) continue;
	j = Scene32::GetAddressX(scn, p);
	i = Scene32::GetAddressY(scn, p);
	k = Scene32::GetAddressZ(scn, p);
	
	newseeds[s].n++;
	newseeds[s].actual.i += i;
	newseeds[s].actual.j += j;
	newseeds[s].actual.k += k;
	newseeds[s].actual.l += scn->array[k][i][j];
      }
      
      for(s = 0; s < nseeds; s++){
	if(newseeds[s].n > 0){
	  newseeds[s].actual.i /= newseeds[s].n;
	  newseeds[s].actual.j /= newseeds[s].n;
	  newseeds[s].actual.k /= newseeds[s].n;
	  newseeds[s].actual.l /= newseeds[s].n;
	}
	else
	  printf("s: %d at (%d, %d, %d)\n", 
		 s, (*seeds)[s].actual.j, (*seeds)[s].actual.i, (*seeds)[s].actual.k);
      }
      
      //------------------------------------------
      //To enforce the placement of seeds inside superpixels:
      if(inside){
	for(p = 0; p < N; p++){
	  s = label->data[p];
	  if(s == NIL) continue;

	  j = Scene32::GetAddressX(scn, p);
	  i = Scene32::GetAddressY(scn, p);
	  k = Scene32::GetAddressZ(scn, p);
	  
	  if( SQUARE(j-newseeds[s].actual.j) + SQUARE(i-newseeds[s].actual.i) + SQUARE(k-newseeds[s].actual.k)  < 
	      SQUARE((*seeds)[s].actual.j-newseeds[s].actual.j) +
	      SQUARE((*seeds)[s].actual.i-newseeds[s].actual.i) +
	      SQUARE((*seeds)[s].actual.k-newseeds[s].actual.k)
	      ){
	    (*seeds)[s].actual.j = j;
	    (*seeds)[s].actual.i = i;
	    (*seeds)[s].actual.k = k;
	  }
	}
	for(s = 0; s < nseeds; s++){
	  newseeds[s].actual.j = (*seeds)[s].actual.j;
	  newseeds[s].actual.i = (*seeds)[s].actual.i;
	  newseeds[s].actual.k = (*seeds)[s].actual.k;
	}
      }
      
      for(s = 0; s < nseeds; s++){
	newseeds[s].last_computed = (*seeds)[s].last_computed;
	newseeds[s].recompute = (*seeds)[s].recompute;
      }
      free(*seeds);
      *seeds = newseeds;
    }
    

  } /*end Superpixels namespace*/
} /*end gft namespace*/

