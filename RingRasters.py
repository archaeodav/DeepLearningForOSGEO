#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 14:12:26 2020

@author: dav
"""
import os
import csv

from osgeo import gdal
from osgeo import osr

import numpy as np

from scipy.ndimage import gaussian_filter
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.filters import minimum_filter
from scipy.ndimage.filters import uniform_filter
from scipy.ndimage.filters import generic_filter

#from skimage.morphology import disk
#from skimage.filters.rank import gradient
#from skimage.draw import disk
#from skimage.draw import circle_perimeter
#from skimage.morphology import binary_erosion
#from skimage.morphology import remove_small_objects


gdal.UseExceptions()

class DtmImage():
    '''Class to process dtm images'''
    def __init__(self,
                 image,
                 overwrite_nd = True):  
        
        #print ('Loading: %s' %(image))
        
        #image
        self.image = None
        
        # window
        self.window = None
        
        # load image
        if self.image is None:
            self.load_image(image)
            
            if overwrite_nd is True:
                self.overwrite_nodata()
        
        # until we say otherwise the window is the image
        if self.window is None:
            self.window = self.image
    
    def overwrite_nodata(self,
                         minval=-100,
                         newval=0):
        
        mask = self.image['dtm']<minval
        
        #print (mask.shape)
        
        self.image['dtm'][mask]=newval
        
    
    def subset_from_point(self,
                          point,
                          dist_xm=150,
                          dist_ym=150):
        ''' method sets current window'''
        # get xy coordinates for window
        # TODO test
        minx = point[0]-dist_xm
        miny = point[1]-dist_ym
        maxx = point[0]+dist_xm
        maxy = point[1]+dist_ym
        
        # slice image using indices of xy coords
        self.window['dtm']=self.image['dtm'][self.get_y_idx(miny):self.get_y_idx(maxy),
                                                         self.get_x_idx(minx):self.get_x_idx(maxx)]
        
        # set spatial transform
        self.window['spat'][0]=minx
        self.window['spat'][3]=maxy
        
    def get_x_idx(self,
                  x):
        #method gets column from x coord
        c = (x-self.image['spat'][0])/self.image['spat'][1]
        return int(c)
        
    def get_y_idx(self,
                  y):
        # method gets row from y coord
        r = (y-self.image['spat'][3])/self.image['spat'][5]
        return int(r)
    
       
    def load_image(self,
                   image):
        ''' method loads image''' 
        
        
        dtm = gdal.Open(image)
    
        dtm_array = dtm.ReadAsArray()
        
            
        spatial = dtm.GetGeoTransform()
        
        proj = dtm.GetProjection()
        
        self.image={'proj':proj,
                    'spat':spatial,
                    'dtm':dtm_array,
                    'name':os.path.split(image)[-1].split('.')[0]}
        
        #print (self.image)
    

    def write_image(self,
                    proj,
                    spat,
                    dtm_array,
                    name,
                    byte_dt=False):        
        '''method writes image'''
        
        #print('ROWS,COLS',dtm_array.shape)
        
       
        #load the driver for the format of choice
        driver = gdal.GetDriverByName("Gtiff") 
        #create an empty output file
        #get the number of bands we'll need:
            
        if len(dtm_array.shape)==2:
            bands = 1    
        else:
            bands = dtm_array.shape[2]
        #print('BANDS OUT', bands)
     
        if byte_dt is True:
            out = driver.Create(name, 
                                dtm_array.shape[1], 
                                dtm_array.shape[0], 
                                bands, 
                                gdal.GDT_Byte,
                                options=['COMPRESS=LZW'])
     
        else:
            out = driver.Create(name, 
                                dtm_array.shape[1], 
                                dtm_array.shape[0], 
                                bands, 
                                gdal.GDT_Float32,
                                options=['COMPRESS=LZW'])
        #define the location using coords of top-left corner
        out.SetGeoTransform(spat)
 
        srs = osr.SpatialReference()
        #get the coodrinate system using the ESPG code
        srs.SetWellKnownGeogCS(proj)
        #set projection of output file 
        out.SetProjection(srs.ExportToWkt())
     
        for band in range(1,bands+1):
            while (band<=bands):
                if bands == 1:
                    data = dtm_array[:,:]
                else:
                    data = dtm_array[:,:,band-1]
                #write values to empty array
                out.GetRasterBand(band).WriteArray(data)    
                #set the no data value
                out.GetRasterBand(band).SetNoDataValue(-9999)
                #apend the statistics to dataset
                out.GetRasterBand(band).GetStatistics(0,1)  
                #print('Saving %s/%s' % (band,bands))
                band = band+1 
        out = None    
     
        #print('Processing of %s complete' % (name))       
     
        return name
    
    
    def rr(self,
           trend_diam = 3,
           pre_diam = 1.2):
        ''' method computes residual relief using difference of gaussians'''
                    
        trend = gaussian_filter(self.window['dtm'],
                                sigma=trend_diam)
                                        
        pre = gaussian_filter(self.window['dtm'], 
                              sigma=pre_diam)
    
        rr = pre-trend
    
        return rr
    
    def multi_res_rr(self,
                     sigmas=[(0.2,1.2),(1.2,3),(3,5)],
                     normalise=True,
                     min_max=[(-0.04,0.04),(-0.08,0.08),(-0.16,0.16)],
                     eight_bit_out=True):
        
        out = None
        
        i = 0
        for s in sigmas:
            rr = self.rr(trend_diam=s[1],pre_diam=s[0])
            #stack.append(rr)
            
            if normalise is True:
                rr = self.normalised_rr(rr=rr,
                                        min_rr=min_max[i][0],
                                        max_rr=min_max[i][1])
            
            
            if out is None:
                out = rr
            else:
                out = np.dstack((out,rr))
                
            i +=1
        
        if eight_bit_out is True and normalise is True:
           out = (out*255).astype(np.uint8)
        
        return out
    
    
    
        
    
    def lrm(self,
            min_thresh = -0.02,
            max_thresh = 0.02):
        
        pass
            
   
    def nhood(self,
              nhood_m):
        '''method calculates neighbourhood in image pixels using distance in m'''
        nhood = int(nhood_m/self.window['spat'][1])
        
        return nhood
    
    def disk_elem(self,
                  rad):
        
        arr = np.zeros((rad*2,rad*2),dtype=np.uint8)
        
        rr,cc = disk((rad,rad),rad)
        
        arr[rr,cc]=1
        
        return arr
    
    def annular_elem(self,
                      rad,
                      width=10):
        
        arr = np.zeros((rad*2,rad*2),dtype=np.uint8)
        
        rr,cc = disk((rad,rad),rad)
        
        brr,bcc =disk((rad,rad),rad-width)
                
        arr[rr,cc]=1
        
        arr[brr,bcc]=0
        
        return arr
    
    def avg(self,
            inarr):
        
        #print (inarr)
        
        return np.average(inarr)
   
    def tpi(self,
            pre_filter=True,
            pre_filter_rad=5,
            nhood_m=150):
        '''method calculates topographic position index'''
        dtm = None
        
        nhood = self.nhood(nhood_m)
        
        fp = self.annular_elem(nhood)
        
        if pre_filter is True:
            dtm = uniform_filter(self.window['dtm'],size=pre_filter_rad)
        else:
            dtm = self.window['dtm']
            
        tpi = dtm-generic_filter(dtm,self.avg,footprint=fp)
            
        #locmin = minimum_filter(self.window['dtm'],size=nhood)
    
        #locmax = maximum_filter(self.window['dtm'],size=nhood)
    
        
        
        return tpi
    
    def grad(self,
             nhood_m=20):
        '''method calculates local gradient '''
        nhood = self.nhood(nhood_m)
        
        g = gradient(self.window['dtm'],self.disk_elem(nhood))
        
        return g
        
    def stack(self):
        '''method calculates stacked raster for ML input- both training data 
        and the candidate features'''
        
        rr = self.rr()
        tpi_100 = self.tpi(nhood_m=50)
        tpi_300 = self.tpi(nhood_m=125)
        
        return np.dstack((rr,tpi_100,tpi_300))
    
    def binary_rr(self,
                  rr=None,
                  threshold = 0.03,
                  min_size = 128,
                  connectivity = 2):
        ''' Class generates binary residual relief for CV operations'''
        
        if rr is None:
            rr = self.rr()
            
        bool_raster = rr>threshold
        
        cleaned = remove_small_objects(bool_raster,min_size, connectivity)
        
        #edges = np.bitwise_xor(cleaned, binary_erosion(cleaned))
        
        return cleaned
    
    def normalise(self):
        '''returns normalised dtm (dtm-dtmmin)'''
        #print (self.window['dtm'])
        n = self.window['dtm']-np.min(self.window['dtm'])
        #print(n)
        
        return n
    
    def normalised_rr(self,
                      rr = None,
                      min_rr=-1,
                      max_rr=1):
        
        '''returns normalised and clipped residual relief'''
        
        if rr is None:
            rr = self.rr()
            
        offset = 0-min_rr
            
        clip_rr = np.clip(rr,min_rr,max_rr)+offset
        
        norm_rr =(clip_rr-(min_rr+offset))/((max_rr+offset)-(min_rr+offset))
        
        return norm_rr
    
    
    
    
class IOhandler():
    def __init__(self,
                 infile,
                 minxcol=5,
                 maxxcol=6,
                 minycol=7,
                 maxycol=8,
                 fcol=0):
    
        self.extents = {}
    
        self.import_extents(infile,
                            minxcol=minxcol,
                            maxxcol=maxxcol,
                            minycol=minycol,
                            maxycol=maxycol,
                            fcol=fcol)  
        
    
    def import_extents(self,
                       file,
                       minxcol=5,
                       maxxcol=6,
                       minycol=7,
                       maxycol=8,
                       fcol=0):
        
        with open(file,'r') as infile:
            reader = csv.reader(infile,delimiter=';')
            next(reader,None)
            for r in reader:
                self.extents[r[fcol]]=(r[minxcol],r[maxxcol],r[minycol],r[maxycol])
                
    def export_extents(self,
                       inlist,
                       outfile):
        
        with open(os.path.join(outfile),'w') as oindex:
            for row in inlist:
                print(row)
                oindex.write(';'.join(row)+'\n')
                
                
class Sampler(IOhandler):
    def __init__(self,
                 infile,
                 points=None):
        
        super().__init__(infile)
        
        self.points = points
        self.points_by_image = {}
    
        pass
    
    def import_points(self,
                      file,
                      xcol=0,
                      ycol=1):
        
        pts = np.genfromtxt(file,delimiter=';')
        
        return pts
    
    def pts_im_query(self):   
        for e in self.extents:
            ext = self.extents[e]
            pts = self.points
            opts = pts[np.logical_and(np.logical_and(pts[:,0]>ext[0],pts[:,0]<ext[1]),
                                      np.logical_and(pts[:,1]>ext[2],pts[:,1]<ext[3]))]
            
            if opts.size > 0:
                self.points_by_image[e]=opts
        
    def get_subsets(self,
                    points,
                    odir):
        
        self.import_points(points)
        self.pts_im_query()
        
        olist =[]
        
        if len(self.points_by_image)>0:
            for f in self.points_by_image:
                n= 0
                im = DtmImage(f)
                pts = self.points_by_image[f]
                fname = os.path.split(f)[-1]
                for p in pts:
                    im.subset_from_point(p)
                    s = im.stack()
                    name = fname.split('.')[0]+'_'+str(n)+'.tif'
                    im.write_image(im.window['proj'],
                                   im.window['spat'],
                                   s,
                                   os.path.join(name))
                    olist.append()
                    
                    n=n+1

class BaseRasters(IOhandler):
    def __init__(self,
                 infile,
                 odir):
        
        super().__init__(infile)
        
        self.odir = odir
        
        rr_idx = []
        
        brr_idx = []
        
        
        rr_dir = os.path.join(odir,'rr')
        if not os.path.exists(rr_dir):
            os.mkdir(rr_dir)
        
        
        brr_dir = os.path.join(odir,'brr')   
        if not os.path.exists(brr_dir):
            os.mkdir(brr_dir)

        for e in self.extents.keys():
            i = DtmImage(e)
            
            if not np.isclose(np.sum(i.image['dtm']),0):
            
                print('rr doing')
                rr = i.rr() 
                rr_name = os.path.join(rr_dir,'rr_'+os.path.split(e)[-1])
                print('rr done')
                
                
                brr = i.binary_rr(rr=rr)
                brr_name = os.path.join(brr_dir,'brr_'+os.path.split(e)[-1])
                
    
                i.write_image(i.window['proj'],
                              i.window['spat'],
                              rr,
                              rr_name)
                
                rr_idx.append((rr_name,
                               self.extents[e][0],
                               self.extents[e][1],
                               self.extents[e][2],
                               self.extents[e][3]))
                              
                i.write_image(i.window['proj'],
                              i.window['spat'],
                              brr,
                              brr_name,
                              byte_dt = True)
                
                brr_idx.append((brr_name,
                                self.extents[e][0],
                                self.extents[e][1],
                                self.extents[e][2],
                                self.extents[e][3]))
                
                self.export_extents(rr_idx, os.path.join(rr_dir,'rr_index.txt'))
                
                self.export_extents(brr_idx, os.path.join(brr_dir,'brr_index.txt'))
            

            
        
    
    
    
        
        

    
    
    
    
    