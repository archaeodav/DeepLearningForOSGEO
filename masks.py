#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 11:18:48 2022

@author: dav
"""

import sys

#sys.path.insert(0, "/work/RingFinding/Preprocessing/")
sys.path.insert(0, "//home/dav/Documents/GitHub/RingFinding/Preprocessing/")

import RingRasters

import os
import csv
import json
import pickle

import torch

from osgeo import gdal
from osgeo import gdalconst
from osgeo import osr
from osgeo import ogr

import numpy as np

from scipy.stats import mode


from sklearn import preprocessing
from torch.utils.data import Dataset

from PIL import Image


class SubsetDtm():
    def __init__(self,
                 vrt,
                 outdir,
                 tempdir):
        # list of root directories 
        self.vrt = vrt
        
        self.outdir = outdir
        
        self.tempdir = tempdir
        
        self.spatial = None
        
        if not os.path.exists(self.tempdir):
            os.mkdir(self.tempdir)
        
        self.get_spatial()
        
        
    def get_spatial(self):
        dtm = gdal.Open(self.vrt)
        
        self.spatial = dtm.GetGeoTransform()
        
        self.proj = dtm.GetProjection()
        
        self.xsize = dtm.RasterXSize
        self.ysize = dtm.RasterYSize
        
        dtm = None
        
    def get_window(self,
                   point,
                   x_m=700,
                   y_m=700):
        
        minx = point[0]-(x_m/2)
        maxx = point[0]+(x_m/2)
        miny = point[1]-(y_m/2)
        maxy = point[1]+(y_m/2)
        
        if minx < self.spatial[0]:
            minx = self.spatial[0]
            
        if maxx > self.spatial[0]+np.abs(self.spatial[1]*self.xsize):
            maxx = self.spatial[0]+np.abs(self.spatial[1]*self.xsize)
            
        if maxy > self.spatial[3]:
            maxy= self.spatial[3]
            print (miny)
        if miny < self.spatial[3]-np.abs(self.spatial[5]*self.ysize):
            miny=self.spatial[3]-np.abs(self.spatial[5]*self.ysize)

        
        return [minx,maxy,maxx,miny]
    
    def get_buffered_bbox(self,
                          bbox,
                          buff_dist):
        
        minx = bbox[0]-buff_dist
        maxx = bbox[1]+buff_dist
        miny = bbox[2]-buff_dist
        maxy = bbox[3]+buff_dist
        
        if minx < self.spatial[0]:
            minx = self.spatial[0]
            
        if maxx > self.spatial[0]+np.abs(self.spatial[1]*self.xsize):
            maxx = self.spatial[0]+np.abs(self.spatial[1]*self.xsize)
            
        if maxy > self.spatial[3]:
            maxy= self.spatial[3]
            print (miny)
        if miny < self.spatial[3]-np.abs(self.spatial[5]*self.ysize):
            miny=self.spatial[3]-np.abs(self.spatial[5]*self.ysize)
            
            
        return [minx,maxy,maxx,miny]

    
    
    #todo refactor
    def get_x_idx(self,
                  image,
                  x):
        #method gets column from x coord
        c = (x-image['spat'][0])/image['spat'][1]
        return int(c)
    #todo refactor
    def get_y_idx(self,
                  image,
                  y):
        # method gets row from y coord
        r = (y-image['spat'][3])/image['spat'][5]
        return int(r)
    
    def subset(self,
               point = None, 
               point_name = None,
               method = 'centroid',
               bbox = None,
               xdim=550, 
               ydim=550,
               cleanup = True,
               stacked = False,
               rtn_file = True,
               rtn_img = False,
               bbox_buffer = 25):
        
        if method == 'centroid':
            projwin = self.get_window(point, x_m=xdim,y_m=ydim)
            
        elif method == 'bbox':
            projwin = self.get_buffered_bbox(bbox,bbox_buffer)
        
        name = os.path.join(self.tempdir,point_name+'.tif')
        
        out_name = os.path.join(self.outdir,point_name+'.tif')
        
        dtm = gdal.Open(self.vrt)
        dtm_o = gdal.Translate(name,
                             dtm,
                             projWin=projwin,
                             xRes=1.6,
                             yRes=1.6)
        
        dtm_o
       
        dtm = None
        dtm_o = None
        
        dtm = RingRasters.DtmImage(name)
        
        prj = dtm.image['proj']
        spt = dtm.image['spat']
        
        if not stacked is True:
            nrr = dtm.normalised_rr()    
            out = dtm.write_image(prj,
                                  spt,
                                  nrr,
                                  out_name)
        
        else:
            norm = dtm.normalise()
            rr = dtm.rr()
            dtm = None
            dtm = gdal.Open(name)
            slope_name = os.path.join(self.tempdir,point_name+'_slope.tif')
            slope = gdal.DEMProcessing(slope_name,dtm,"slope")
            slope = None
            slope = RingRasters.DtmImage(slope_name).image['dtm']
            stacked = np.dstack((norm,rr,slope))
            out = dtm.write_image(prj,
                                  spt,
                                  stacked,
                                  out_name)

            if cleanup is True:
                os.remove(slope_name)
            
        if cleanup is True:
            os.remove(name)
        
        if rtn_file is True:  
            return out,out_name
       
        elif rtn_img is True:
            return out_name,dtm
        
        else:
            return nrr
    
class PolyMask(SubsetDtm):
    def __init__(self,
                 vrt,
                 outdir,
                 tempdir,
                 country,
                 mask_dir = None,
                 classes = None):
        
        # Instatiate a SubsetDtm instance
        super().__init__(vrt,outdir,tempdir)
        
        # Dict to store training masks
        self.training_masks = {}
        
        #Get the country for the filename if nothing else
        self.country = country
        
        # Directory to contain mask images
        self.mask_dir = mask_dir
        
        
        ''' We need integer defintions in the masks, so we'll use these defaults 
        if none are supplied to the init method'''
        
        if self.classes is None:
        
            self.classes = {'Fortification':1,
                            'FortificationOther':2,
                            'Pingo':3,
                            'Kettlehole':4,
                            'Modern':5,
                            'Kame':6}
            
    def polys(self,
              shapefile,
              generate_masks = True,
              driver = "ESRI Shapefile",
              generate_im_polys = False,
              generate_instance_masks = True):
        ''' Mehtod hnadles polygon training data. Generate masks makes pixel
        mask images for pytorch etc'''
        
        # this stuff loads the vector polys
        driver = ogr.GetDriverByName(driver)
        dataSource = driver.Open(shapefile, 1)
        layer = dataSource.GetLayer()
        
        layerDefinition = layer.GetLayerDefn()
        
        fields = []
        
        # check which fields are there
        for i in range(layerDefinition.GetFieldCount()):
            fields.append(layerDefinition.GetFieldDefn(i).GetName())
            
            
        ''' if we don't have a class integer identifier make one using our default 
        OR use the supplied dict in the class init method'''
        if not "Class_int" in fields:
            new_field = ogr.FieldDefn("Class_int", ogr.OFTIntiger)
            new_field.SetWidth(8)
            layer.CreateField(new_field)
            
            for feature in layer:
                classification = feature.GetField('Class')
                feature.SetField("Class_Int", self.classes[classification])
                layer.SetFeature(feature)
        

        ''' if we don't have a instance integer identifier make one'''                
        if not "Inst_ID" in fields:
            new_field = ogr.FieldDefn("Inst_id", ogr.OFTIntiger)
            new_field.SetWidth(8)
            layer.CreateField(new_field)
            
            i = 1
            for feature in layer:
                classification = feature.GetField('Inst_id')
                feature.SetField("Inst_id", i)
                layer.SetFeature(feature)
            i+=1
            ''' This is potentiall ugly- but I think the presumption that there'll
            be fewer than 255 instances per training subset is sound and it 
            will simplify IO as we can handle everything as a byte type'''
            if i > 255:
                i = 1
                
        #TODO- add some kind of error handling to this bit to check that there
        # are the right number of classes in the data and supplied class deifnition
        
        #TODO- add something to override existing deifnitions
        
        
        dataSource = None  
        
        # Now get geometries from polygons
        dataSource = driver.Open(shapefile, 0)
        layer = dataSource.GetLayer()
                
        i = 0
        
        #iterate through features
        for feature in layer:
            points = []
            pixels = []
            # get the geometry
            geom = feature.GetGeometryRef()
            ring = geom.GetGeometryRef(0)
            pts = ring.GetPointCount()
            
            # get the vertices
            for p in range(pts):
                x,y,z = ring.GetPoint(p)
                points.append((x,y))
                
            # get the centroid of the poly
            centroid = geom.Centroid()
            
            # get the centroid of the poly
            x,y,z = centroid.GetPoint()
            centroid = (x,y)
            
            # get the encelope (bbox) of the poly
            envelope = geom.GetEnvelope()
            
            
            # Get attribute
            classification = feature.GetField('Class')
            instance= feature.GetField('Inst_id')
            
            # Derive point name
            point_name = '%s_%s_%s' %(self.country,classification,str(i))
            
            
            # Ad the point to the dict
            self.training_masks[point_name]={}
            self.training_masks[point_name]['class']=classification
            
            ''' If we're doing masks as polygons (for the matterport 
            implementation) then do them'''
            if generate_im_polys is True:
                
                # subset training image from vrt
                s = self.subset(point=centroid, 
                                point_name=point_name,
                                rtn_file=False,
                                rtn_img=True)
                
                # project pintss as pixel coordinates
                for p in points:
                    r = self.get_y_idx(s[1], p[1])
                    c = self.get_x_idx(s[1],p[0])
                    pixels.append(r,c)
            
                self.training_masks[point_name]['mask_poly']=pixels    
            
            ''' If we're using mask images for pytorch then do them '''
            if generate_masks is True:
                
                ''' subset the training image with a buffer around the bbox'''
                s = self.subset(method = 'bbox', 
                                bbox = envelope,
                                rtn_file=False,
                                rtn_img=True,
                                bbox_buffer=25)
            
            # append path to subset to dict
            self.training_masks[point_name]['image']=s[0]
            
            self.training_masks[point_name]['class_mask']=self.int_mask(s,
                                                                        point_name,
                                                                        layer,
                                                                        'Class_int',
                                                                        '_class_mask.tif')
            
            self.training_masks[point_name]['instance_mask']=self.int_mask(s,
                                                                           point_name,
                                                                           layer,
                                                                           'Inst_id',
                                                                           '_instance_mask.tif')

            i =+1
        dataSource = None
        
    def int_mask(self,
                 s,
                 point_name,
                 layer,
                 field,
                 postfix):
        
        # generate name for output file
        mask_name = os.path.join(self.mask_dir,point_name+postfix)
        
        # make a target dataset
        target_ds = gdal.GetDriverByName('GTiff').Create(mask_name, 
                                                         s[1].shape[1], 
                                                         s[1].shape[0], 
                                                         1, 
                                                         gdal.GDT_Byte)
        # set geo transform
        target_ds.SetGeoTransform(s[1].image.spatial)
        
        band = target_ds.GetRasterBand(1)
        
        tfield = "Attribute="+field
        
        # Rasterize the polygon
        gdal.RasterizeLayer(target_ds, [band], layer, options=[tfield])
        
        
        target_ds.FlushCache()
        del target_ds
        
        return mask_name
        
        
            
class MaskOutDict():
    def __init__(self,
                 odir,
                 oname='training_data.json',
                 indict = None):
        
        # if there's no supplied file to append to instatiate a dict
        if indict is None:
            self.outdict = {}
        
        # otherwise load the input file as a dict
        else:
            if type(indict) is dict:
                self.outdict = indict
            else:
                self.outdict = {}
                self.load_dict(indict)
            
    
    # loads dict from json file
    def load_dict(self,
                  infile):
        
        infiles = []
        
        data = {}
        
        if type(infile) is list:
            infiles = infile
            
        else:
            infiles.append(infile)
            
        for f in infiles:
            with open(f, 'r') as ifile:
                d = json.load(ifile)
                data.update(d)
                
        self.outdict.update(data)
        
    # writes dict to json file 
    def write_dict(self,
                   outfile):
        
        with open(outfile, 'w') as ofile:
            json.dump(self.outdict,
                      ofile,
                      sort_keys=True,
                      indent=4,
                      ensure_ascii=False)
            
            ofile.close()
    
    
    def append_data(self,
                    data):
        
        self.outdict.update(data)
        
class ptorch_masks(Dataset):
    '''https://pytorch.org/tutorials/intermediate/torchvision_tutorial.html'''
    def __init__(self,
                 mask_dict,
                 transforms):
        
        self.data = mask_dict
        
        self.images = list(self.data.keys())
        
        self.transforms = transforms
        
        
    def __len__(self,
                ):
        return len(self.images)
    
    def __getitem__(self,
                    idx):
        
        im = self.data[self.images[idx]]
        
        image = Image.open(self.data[im]['image']).convert('RGB')
        
        class_mask = Image.open(self.data[im]['class_mask'])
        
        instance_mask = Image.open(self.data[im]['instance_mask'])
        
        m = np.array(instance_mask)
        
        c = np.array(class_mask)
        
        obj_ids = np.unique(m)
        
        obj_ids = obj_ids[1:]
        
        masks = m == obj_ids[:, None, None]
        
        num_objs = len(obj_ids)
        
        boxes = []
        labels = []
        
        for i in range(num_objs):
            pos = np.where(masks[i])
            xmin = np.min(pos[1])
            xmax = np.max(pos[1])
            ymin = np.min(pos[0])
            ymax = np.max(pos[0])
            boxes.append([xmin, ymin, xmax, ymax])
            
            labels.append(mode(c[masks[i]])[0][0])
            
        # convert everything into a torch.Tensor
        boxes = torch.as_tensor(boxes, dtype=torch.float32)
        # there is only one class
        labels = torch.as_tensor(labels, dtype=torch.int64)
        masks = torch.as_tensor(masks, dtype=torch.uint8)

        image_id = torch.tensor([idx])
        area = (boxes[:, 3] - boxes[:, 1]) * (boxes[:, 2] - boxes[:, 0])
        # suppose all instances are not crowd
        iscrowd = torch.zeros((num_objs,), dtype=torch.int64)

        target = {}
        target["boxes"] = boxes
        target["labels"] = labels
        target["masks"] = masks
        target["image_id"] = image_id
        target["area"] = area
        target["iscrowd"] = iscrowd

        if self.transforms is not None:
            image, target = self.transforms(image, target)

        return image, target
            
class MaskSampler(SubsetDtm):
    def __intit__(self,
                 vrt,
                 outdir,
                 tempdir):
        super().__init__(vrt,outdir,tempdir)
        
    def get_sample(self,
                   x_size=5000,
                   y_size=5000,
                   buffer=250):
        
        pass
    
    
        
        
       
            


    
        
