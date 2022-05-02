#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 12:58:48 2020

@author: dav

Script for the retiling and resampling rasters prior to derivation of inputs 
for the ring finding algorothims

"""

import os
import shutil
from osgeo import gdal
import sys

class InputOutput():
    def __init__(self,
                 indir,
                 temp_dir,
                 odir,
                 epsg):
        
        
        self.indir = indir
        self.temp_dir = temp_dir 
        self.odir = odir
        self.olist = []
        
        self.idxfile=None
        
        self.epsg=epsg
        
        
    def infile_list(self):
        for file in os.listdir(self.indir):
            print(file)
            ext = file.split('.')[-1]
            if ext == 'tif' or ext == 'tiff' or ext == 'asc' or ext=='xyz':
                print(file)
                self.olist.append(os.path.join(self.indir,file))

    def get_raster_res(self):
        s =self.read_image(self.olist[0])
        
        return s[0][1]
        
    def get_raster_proj(self):
        s =self.read_image(self.olist[0])
        
        return s[1]
    
    def gen_file_list(self,
                      tdir=None,
                      olist=None): 
        
        if tdir is None:
            tdir = self.indir
        if olist is None:
            olist = self.olist
            
        with open(os.path.join(tdir,'infiles.txt'),'w') as ofile:
            for l in olist:
                ofile.write(l+'\n')
        
        return os.path.join(tdir,'infiles.txt')     

    def build_vrt(self,
                  tdir=None,
                  olist=None,
                  srs_string=None,
                  opath=None):
        if tdir is None:
            tdir = self.temp_dir
            
        if olist is None:
            olist = self.olist
            
        if opath is None:
            opath=tdir
            
            
        if srs_string is None:
            srs = '-a_srs %s ' %(self.epsg)
    
        
        flist = self.gen_file_list(olist=olist)
        
        oname = os.path.join(opath,os.path.split(tdir)[-1]+'.vrt')
        cmd = 'gdalbuildvrt %s-input_file_list %s %s' %(srs,flist,oname)
        os.system(cmd)
        
        return oname

    def read_image(self,img):
        dtm = gdal.Open(img)
        #dtm_array = dtm.ReadAsArray()
        spatial = dtm.GetGeoTransform()
        proj = dtm.GetProjection()
        return spatial, proj
                
    def retile_rasters(self, 
                       tile_size_m = 15000,
                       overlap_m=275):
        
        psize = self.get_raster_res()
        tile_size_p = int(tile_size_m/psize)
        overlap_p = int(overlap_m/psize)
        command = 'gdal_retile -v -targetDir %s -ps %s %s -overlap %s' %(self.temp_dir,
                                                                            tile_size_p,
                                                                            tile_size_p,
                                                                            overlap_p)
        co = '-of GTiff -co "BIGTIFF=IF_SAFER" -co "NUM_THREADS=ALL_CPUS" -co "COMPRESS=LZW"'
        
        #self.idxfile = os.path.join(self.temp_dir,'index.csv')
        self.idxfile = os.path.join(self.temp_dir,'index.csv')
        csvf = '-csv %s' %('index.csv')
        #csvf = '-csv %s' %(self.idxfile)
        print (csvf)
        
        
        cmd = '%s %s %s %s' %(command,co,csvf,self.build_vrt())
        
        print (cmd)
        
        os.system(cmd)

        self.overlap_index(overlap_m)

    def overlap_index(self,overlap_m):
        out_data = []
        header = ['file',
                  'minx',
                  'maxx',
                  'miny',
                  'maxy',
                  'd_minx',
                  'd_maxx',
                  'd_miny',
                  'd_maxy']
        
        out_data.append(header)

        with open(self.idxfile,'r') as index_file:
            for row in index_file.readlines():
                tilename,minx,maxx,miny,maxy = row.split(';')
                d_minx = float(minx)+overlap_m
                d_maxx = float(maxx)-overlap_m
                d_miny = float(miny)+overlap_m
                d_maxy = float(maxy)-overlap_m
                
                out_data.append((os.path.join(self.odir,tilename),
                                str(minx).strip('\n'),
                                str(maxx).strip('\n'),
                                str(miny).strip('\n'),
                                str(maxy).strip('\n'),
                                str(d_minx).strip('\n'),
                                str(d_maxx).strip('\n'),
                                str(d_miny).strip('\n'),
                                str(d_maxy).strip('\n')))
        
        with open(os.path.join(self.odir,'index_olap.csv'),'w') as oindex:
            for row in out_data:
                print(row)
                oindex.write(';'.join(row)+'\n')
              
    def resample(self, target_res=1.6):
        for file in os.listdir(self.temp_dir):
            if file.split('.')[-1]=='tif':
                ifile  = os.path.join(self.temp_dir,file)
                ofile  = os.path.join(self.odir,file)
                c = 'gdal_translate -of GTiff -r bilinear -tr %s %s' %(target_res,target_res) 
                co = '-co "TILED=YES" -co "BIGTIFF=IF_SAFER" -co "NUM_THREADS=ALL_CPUS" -co "COMPRESS=LZW' 
                cmd = '%s %s %s %s' %(c,co,ifile,ofile)
                os.system(cmd)
                
    def stopit_tidyup(self):
        shutil.rmtree(self.temp_dir) 
        os.mkdir(self.temp_dir)

    def process_rasters(self):
        '''one ting to rule them all'''
        self.infile_list()
        self.build_vrt()
        self.retile_rasters()
        self.resample()
        self.stopit_tidyup()
        pass

    def build_td_vrt(self,
                     srs=None):
        
        if srs is None:
            try:
                epsg_file = os.path.split(self.indir)[0]+'epsg.txt'
                with open(epsg_file, 'r') as e:
                    srs = e.read()
            except:
                raise Exception("Can't read epsg file")
                
        self.infile_list()
        self.build_vrt(tdir=self.indir)
            
                    
        
    
class MultiDirs():
    def __init__(self,
                 dirlist=None):
        self.inlist=[]    
        if not dirlist is None:
            self.inlist = self.loaddirlist(dirlist)
        else:
            raise Exception()
            
    def loaddirlist(self,
                    dirlist):
        with open(dirlist, 'r') as i:
            for row in i.readlines():
                d = row.strip('\n')
                if os.path.exists(d):
                    self.inlist.append()
                    
    def run_multidir(self,
                     in_prefix = 'dtm_raw',
                     temp_prefix = 'temp',
                     out_prefix = 'out'):
        
        for d in self.inlist:
            indir = os.path.join(d,in_prefix)
            tempdir = os.path.join(d,temp_prefix)
            outdir = os.path.join(d,out_prefix)
            
            InputOutput(indir,tempdir,outdir)
    
if __name__ == '__main__':
    t = InputOutput(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
    t.process_rasters()