
#ifndef _GFT_SUPERPIXELS_H_
#define _GFT_SUPERPIXELS_H_

#include "gft_common.h"
#include "gft_image32.h"
#include "gft_cimage.h"
#include "gft_color.h"
#include "gft_adjrel.h"
#include "gft_filtering.h"
#include "gft_queue.h"
#include "gft_heap.h"
#include "gft_scene32.h"
#include "gft_adjrel3.h"


namespace gft{
  namespace Superpixels{

    typedef struct _feature{
      int j;
      int i;
      float l;
      float a;
      float b;
    } Feature;

    typedef struct _spixelseed{
      Feature actual;
      Feature last_computed;
      bool recompute;
      int n;
    } SPixelSeed;

    //--------------------------

    float ColorDistance(CImage::CImage *cimg_lab, 
			SPixelSeed seed, int q);

    float LastColorDistance(SPixelSeed seed);
    float LastSpatialDistance(SPixelSeed seed);

    void RunSPixelsByIFT(Image32::Image32 *grad,
			 CImage::CImage *cimg_lab,
			 Image32::Image32 *label,
			 Heap::Heap *Q,
			 Image32::Image32 *pred,
			 float *cost,
			 float alpha,
			 float beta,
			 SPixelSeed *seeds, int nseeds);
    
    Image32::Image32 *IFT_SLIC(CImage::CImage *cimg,
			       int k, float alpha, float beta,
			       float err_dc, float err_ds, int itMax);
    Image32::Image32 *IFT_SLIC(Image32::Image32 *img,
			       int k, float alpha, float beta,
			       float err_dc, float err_ds, int itMax);
    //--------------------------
    float WeightedDistanceMeasure(CImage::CImage *cimg_lab, 
				  SPixelSeed seed, int i, int j,
				  float m, float S);

    Image32::Image32 *mySLIC(CImage::CImage *cimg, 
			     int k, float m);
    Image32::Image32 *mySLIC(Image32::Image32 *img, 
			     int k, float m);
    //--------------------------
    SPixelSeed *GetSPixelsSeeds(CImage::CImage *cimg_lab,
				Image32::Image32 *grad,
				int k,
				int *nseeds);
    void Move2LowerGradient(Image32::Image32 *grad, 
			    int *j, int *i);

    void UpdateSPixelsSeeds(CImage::CImage *cimg_lab,
			    Image32::Image32 *label,
			    SPixelSeed **seeds, int nseeds,
			    int inside);

    int GetNumberOfSuperPixels(Image32::Image32 *label);


    void Postprocessing_CC(Image32::Image32 *label,
			   SPixelSeed *seeds, int nseeds);
    void Postprocessing(Image32::Image32 *label,
			SPixelSeed *seeds, int nseeds);


    //--------------------3D--------------------------------
    typedef struct _feature3D{
      int j;
      int i;
      int k;
      float l;
    } Feature3D;

    typedef struct _svoxelseed{
      Feature3D actual;
      Feature3D last_computed;
      bool recompute;
      int n;
    } SVoxelSeed;

    
    float FeatureDistance(Scene32::Scene32 *scn, 
			  SVoxelSeed seed, int q);
    
    float LastFeatureDistance(SVoxelSeed seed);
    float LastSpatialDistance(SVoxelSeed seed);
    
    void RunSVoxelsByIFT(Scene32::Scene32 *scn,
			 Scene32::Scene32 *label,
			 Heap::Heap *Q,
			 Scene32::Scene32 *pred,
			 float *cost,
			 float alpha,
			 float beta,
			 SVoxelSeed *seeds, int nseeds);
    
    Scene32::Scene32 *IFT_SLIC(Scene32::Scene32 *scn,
			       int k, float alpha, float beta,
			       float err_df, float err_ds, int itMax);

    SVoxelSeed *GetSVoxelsSeeds(Scene32::Scene32 *scn,
				int k,
				int *nseeds);
    
    void UpdateSVoxelsSeeds(Scene32::Scene32 *scn,
			    Scene32::Scene32 *label,
			    SVoxelSeed **seeds, int nseeds,
			    int inside);
    
  } //end Superpixels namespace
} //end gft namespace


#endif

