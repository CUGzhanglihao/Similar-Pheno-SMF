'''
每景影像每景影像进行处理
'''

import ee
import geemap
import pandas as pd
import json
import  threading
import time
import numpy as np
import os
from tqdm import tqdm,trange

def maskL457sr(image):
  # Bit 0 - Fill
  # Bit 1 - Dilated Cloud
  # Bit 2 - Unused
  # Bit 3 - Cloud
  # Bit 4 - Cloud Shadow
  # Bit 5 - Snow
  # Bit 6 - Clear
  qaMask = image.select('QA_PIXEL').bitwiseAnd(63).eq(0)
  saturationMask = image.select('QA_RADSAT').eq(0)

  # Apply the scaling factors to the appropriate bands.
  opticalBands = image.select('SR_B3','SR_B4').multiply(0.0000275).add(-0.2)

  # Replace the original bands with the scaled ones and apply the masks.
  return opticalBands \
      .updateMask(qaMask) \
      .updateMask(saturationMask) \
      .copyProperties(image,image.propertyNames())


def GetL5NDVI(image):
    NDVI = image.normalizedDifference(['SR_B4','SR_B3']).rename('L5_NDVI')
    return image.addBands(NDVI).copyProperties(image,image.propertyNames())

# def GetL5NDVI_Trajectory(image):
#     NDVI = image.normalizedDifference(['SR_B4','SR_B3']).multiply(10000).toInt16().rename('Calibration_NDVI')
#     return image.addBands(NDVI).copyProperties(image,image.propertyNames())

def GetL7NDVI(image):
    NDVI = image.normalizedDifference(['SR_B4','SR_B3']).multiply(10000).toInt16().rename('Calibration_NDVI')
    return image.addBands(NDVI).copyProperties(image,image.propertyNames())

# This example demonstrates the use of the Landsat 8 Collection 2, Level 2
# QA_PIXEL band (CFMask) to mask unwanted pixels.

def maskL8sr(image):
  # Bit 0 - Fill
  # Bit 1 - Dilated Cloud
  # Bit 2 - Unused
  # Bit 3 - Cloud
  # Bit 4 - Cloud Shadow
  # Bit 5 - Snow
  # Bit 6 - Clear
  qaMask = image.select('QA_PIXEL').bitwiseAnd(63).eq(0)
  saturationMask = image.select('QA_RADSAT').eq(0)

  # Apply the scaling factors to the appropriate bands.
  opticalBands = image.select('SR_B4','SR_B5').multiply(0.0000275).add(-0.2)
  # NDVI = image.normalizedDifference(['SR_B5','SR_B4'])
  return opticalBands \
      .updateMask(qaMask) \
      .updateMask(saturationMask) \
      .copyProperties(image,image.propertyNames())

def GetL8NDVI(image):
    NDVI = image.normalizedDifference(['SR_B5','SR_B4']).rename('L8_NDVI')
    return image.addBands(NDVI).copyProperties(image,image.propertyNames())

# def GetL8NDVI_Trajectory(image):
#     NDVI = image.normalizedDifference(['SR_B5','SR_B4']).multiply(10000).toInt16().rename('Calibration_NDVI')
#     return image.addBands(NDVI).copyProperties(image,image.propertyNames())

from geemap import ml

#添加DOY和Year波段
def Add_DOY_Year(image):
    NDVI = image.select('Calibration_NDVI')
    DOY = ee.Number.parse(image.date().format('DDD')).toInt16()
    Year = ee.Number.parse(image.date().format('YYYY')).toInt16()
    Image_DOY = NDVI.multiply(0).add(ee.Number.parse(image.date().format('DDD'))).toInt16().rename('DOY')
    # Image_Year = NDVI.multiply(0).add(ee.Image(ee.Number.parse(image.date().format('YYYY'))).toInt16()).rename('Year')
    Props = {'DOY':DOY,'Year':Year} 
    return image.addBands(Image_DOY) \
    .set(Props) \
    .copyProperties(image,image.propertyNames())#.addBands(Image_Year)

def Get_Sensor_Calibration_NDVI_Collection(Start_Year,End_Year,SubTile):

    MultiYear_NDVI = ee.Image("projects/ee-zhanglihaobnu/assets/MultiYear_90_Quantile_Mask") \
        .select('NDVI_mean').clip(ee.Feature(SubTile).geometry())
    
    Urban_Crop_Label = ee.Image("projects/ee-zhanglihaobnu/assets/GlobalLand30_Mosaic_ClipTP_Combine_Urban_Crop") \
    .select('b1').clipToBoundsAndScale(ee.Feature(SubTile).geometry(),scale=30)

    Valid_Observations_Mask = ee.Image(MultiYear_NDVI.And(Urban_Crop_Label.eq(0)))
    # 创建一个空的图像集合
    Combine_Collection = ee.List([])
    # TP_StudyArea = geemap.shp_to_ee('E:/Zhanglihao/Shrub_Encroachment/Test/Samples_Center_Point.shp')
    L5_Image_Collection1 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2") \
                        .filterBounds(ee.Feature(SubTile).geometry()) \
                        .filterDate('1986-01-01','1994-01-01') \
                        .filter(ee.Filter.calendarRange(Start_Year,End_Year,'year')) \
                        .map(lambda img: img.updateMask(Valid_Observations_Mask).clip(ee.Feature(SubTile).geometry())) \
                        .map(maskL457sr) \
                        .map(GetL5NDVI) \
                        .select('L5_NDVI')
                        # .filterMetadata('WRS_PATH','equals',WRS_PATH) \
                        # .filterMetadata('WRS_ROW','equals',WRS_ROW) \
    L5_Image_Collection2 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2") \
                            .filterBounds(ee.Feature(SubTile).geometry()) \
                            .filterDate('1997-01-01','2001-01-01') \
                            .filter(ee.Filter.calendarRange(Start_Year,End_Year,'year')) \
                            .map(lambda img: img.updateMask(Valid_Observations_Mask).clip(ee.Feature(SubTile).geometry())) \
                            .map(maskL457sr) \
                            .map(GetL5NDVI) \
                            .select('L5_NDVI')

    L5_Image_Collection3 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2") \
                            .filterBounds(ee.Feature(SubTile).geometry()) \
                            .filterDate('2002-01-01','2012-01-01') \
                            .filter(ee.Filter.calendarRange(Start_Year,End_Year,'year')) \
                            .map(lambda img: img.updateMask(Valid_Observations_Mask).clip(ee.Feature(SubTile).geometry())) \
                            .map(maskL457sr) \
                            .map(GetL5NDVI) \
                            .select('L5_NDVI')

    L5_Image_Collection = (L5_Image_Collection1.merge(L5_Image_Collection2).merge(L5_Image_Collection3))

    #选取Landsat7影像1999--2019
    L7_Image_Collection = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2") \
                            .filterBounds(ee.Feature(SubTile).geometry()) \
                            .filterDate('1999-01-01','2020-01-01') \
                            .filter(ee.Filter.calendarRange(Start_Year,End_Year,'year')) \
                            .map(lambda img: img.updateMask(Valid_Observations_Mask).clip(ee.Feature(SubTile).geometry())) \
                            .map(maskL457sr) \
                            .map(GetL7NDVI) \
                            .select('Calibration_NDVI')

    #选取Landsat8影像
    L8_Image_Collection = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2") \
                            .filterBounds(ee.Feature(SubTile).geometry()) \
                            .filterDate('2013-01-01','2021-01-01') \
                            .filter(ee.Filter.calendarRange(Start_Year,End_Year,'year')) \
                            .map(lambda img: img.updateMask(Valid_Observations_Mask).clip(ee.Feature(SubTile).geometry())) \
                            .map(maskL8sr) \
                            .map(GetL8NDVI) \
                            .select('L8_NDVI')

    #导入DEM数据,扩充60m,防止边界区域无法算出slope和aspect
    DEM_Data = ee.Image('NASA/NASADEM_HGT/001')  \
                .select('elevation') \
                .clip(ee.Feature(SubTile)) \
                .rename('DEM')

    DEM_Slope = ee.Image(ee.Terrain.slope(DEM_Data)).clip(ee.Feature(SubTile).geometry()).rename('DEM_Slope')
    DEM_Aspect = ee.Image(ee.Terrain.aspect(DEM_Data)).clip(ee.Feature(SubTile).geometry()).rename('DEM_Aspect')
    DEM_Data = DEM_Data.clip(ee.Feature(SubTile).geometry())

    #导入随机森林回归模型
    L87_RF_FC = ee.FeatureCollection("projects/ee-zhanglihaobnu/assets/Shrub_Encroachment/L87_RF_RegMod_Final")
    L87_RF_RegMod = ml.fc_to_classifier(L87_RF_FC)

    L57_RF_FC = ee.FeatureCollection("projects/ee-zhanglihaobnu/assets/Shrub_Encroachment/L57_RF_RegMod_Final")
    L57_RF_RegMod = ml.fc_to_classifier(L57_RF_FC)

    def L5_Calibration_NDVI(image):
    
        Date_1984 = ee.Date.parse("yyyy-MM-dd HH:mm:ss","1984-01-01 00:00:00")
        Image_Date0 = image.date()
        Image_DOY = ee.Image(ee.Number.parse(image.date().format('DDD')).toInt16()).clip(image.geometry()).rename('DOY')
        Image_Year = ee.Image(ee.Number.parse(image.date().format('YYYY')).toInt16()).clip(image.geometry()).rename('Year')
        Image_Date = ee.Image(Image_Date0.difference(Date_1984,'second')).clip(image.geometry()).rename('L5_Date')
        Image_LonLat = image.pixelLonLat()
        Image_Lon = Image_LonLat.select('longitude').rename('Longitude')
        Image_Lat = Image_LonLat.select('latitude').rename('Latitude')
        
        image = image.addBands(DEM_Data).addBands(DEM_Slope).addBands(DEM_Aspect).addBands(Image_Lon) \
        .addBands(Image_Lat) .addBands(Image_DOY).addBands(Image_Year).addBands(Image_Date)
        feature_names = ['L5_NDVI','DEM','Longitude','Latitude','L5_Date','Year','DOY','DEM_Slope','DEM_Aspect']
        Origion_NDVI = ee.Image(image.select('L5_NDVI')).multiply(10000).toInt16()
        Calibration_L5_NDVI = image.select(feature_names).classify(L57_RF_RegMod)
        Calibration_L5_NDVI = (Calibration_L5_NDVI.multiply(10000)).toInt16().clip(image.geometry())
        return ee.Image(Calibration_L5_NDVI.rename('Calibration_NDVI')).addBands(Origion_NDVI) \
            .copyProperties(image,image.propertyNames())

    def L8_Calibration_NDVI(image):
    
        Date_1984 = ee.Date.parse("yyyy-MM-dd HH:mm:ss","1984-01-01 00:00:00")
        Image_Date0 = image.date()
        Image_DOY = ee.Image(ee.Number.parse(image.date().format('DDD')).toInt16()).clip(image.geometry()).rename('DOY')
        Image_Year = ee.Image(ee.Number.parse(image.date().format('YYYY')).toInt16()).clip(image.geometry()).rename('Year')
        Image_Date = ee.Image(Image_Date0.difference(Date_1984,'second')).clip(image.geometry()).rename('L8_Date')
        Image_LonLat = image.pixelLonLat()
        Image_Lon = Image_LonLat.select('longitude').rename('Longitude')
        Image_Lat = Image_LonLat.select('latitude').rename('Latitude')
        
        image = image.addBands(DEM_Data).addBands(DEM_Slope).addBands(DEM_Aspect).addBands(Image_Lon) \
        .addBands(Image_Lat) .addBands(Image_DOY).addBands(Image_Year).addBands(Image_Date)
        feature_names = ['L8_NDVI','DEM','Longitude','Latitude','L8_Date','Year','DOY','DEM_Slope','DEM_Aspect']

        Origion_NDVI = ee.Image(image.select('L8_NDVI')).multiply(10000).toInt16()

        Calibration_L8_NDVI = image.select(feature_names).classify(L87_RF_RegMod)
        Calibration_L8_NDVI = (Calibration_L8_NDVI.multiply(10000)).toInt16().clip(image.geometry())
        return ee.Image(Calibration_L8_NDVI).rename('Calibration_NDVI').addBands(Origion_NDVI) \
                .copyProperties(image,image.propertyNames())

    Calibration_L5_NDVI_Collection = L5_Image_Collection.map(L5_Calibration_NDVI)
    Calibration_L8_NDVI_Collection = L8_Image_Collection.map(L8_Calibration_NDVI)

    Landsat_NDVI_Collection = (L7_Image_Collection.merge(Calibration_L5_NDVI_Collection).merge(Calibration_L8_NDVI_Collection))
    
    return ee.ImageCollection(Landsat_NDVI_Collection)


def Export_Image_toLocal(Start_Year,SubTile_Index,Path_Dir):
    start = time.time()

    ee.Initialize()
    Dir = Path_Dir + "SubTile_"+str(SubTile_Index)+"\\"+str(Start_Year)+"\\"
    TP_SubTiles = ee.FeatureCollection('projects/ee-zhanglihaobnu/assets/Shrub_Encroachment/TP_20SubTiles_50KM')
    SubTile_Nums = TP_SubTiles.size().getInfo()
    SubTiles_List = TP_SubTiles.toList(SubTile_Nums)
    
    SubTile = ee.Feature(SubTiles_List.get(SubTile_Index))    
    if not os.path.exists(Dir):
        os.makedirs(Dir)
        print("创建文件夹："+Dir)

    End_Year = Start_Year
    Landsat_NDVI_Collection = Get_Sensor_Calibration_NDVI_Collection(Start_Year,End_Year,SubTile)
    Images_Nums = int(Landsat_NDVI_Collection.size().getInfo())
    print(Images_Nums,'Images')
    Images_ID_Filename = Dir+'Images_ID.json'
    if not os.path.exists(Images_ID_Filename):
        print('创建Images_ID列表文件:',Images_ID_Filename)
        if Images_Nums > 0:
            Images_Index = Landsat_NDVI_Collection.aggregate_array("system:index").getInfo()
            Images_ID = ['L'+(str(i).split('L')[-1]) for i in Images_Index]
        else:
            Images_ID = []
        #写出Images_ID，可用于后续快速检查文件是否下载完全。
        Df_Images_ID = pd.DataFrame({'Images_ID':Images_ID})
        Df_Images_ID.to_json(path_or_buf=Dir+'Images_ID.json',orient='records')
    else:
        if Images_Nums > 0:
            Df_Images_ID = pd.read_json(path_or_buf=Images_ID_Filename,orient='records')
            Images_ID = (Df_Images_ID['Images_ID'].values).tolist()
        else:
            Images_ID = []

    Left_Images_ID = FindFileName(Images_ID,Dir)

    # NDVI_List = Landsat_NDVI_Collection.toList(Images_Nums)
    i=0
    if len(Left_Images_ID) != 0:
            while(len(Left_Images_ID) >0):

                for Image_ID in tqdm(Left_Images_ID):
                    # image = ee.Image(NDVI_List.get(i))
                    image = ee.Image((Landsat_NDVI_Collection.filter(ee.Filter.stringContains("system:index",Image_ID))).first())
                    Date = (image.date().getInfo())['value']
                    # ID = 'L'+str(image.get("system:index").getInfo()).split('L')[-1]
                    DOY = (time.gmtime(Date/1000)).tm_yday #image.date().format('DDD').getInfo()
                    Filename = Dir + Image_ID+"_"+str(Date)+"_"+str(Start_Year)+"_"+str(DOY)+ ".tif"
                    if not os.path.isfile(Filename):
                        try:
                            geemap.download_ee_image(image,Filename,region=SubTile.geometry(),scale=30, \
                                                    dtype='int16',crs='EPSG:4326')
                        except Exception as e:
                            print("An error occurred",str(e))
                            pass
                        
                Left_Images_ID = FindFileName(Images_ID,Dir)
                i+=1
                end = time.time()
                print('Main process time:',(end-start)/60)
                print('Loop',i,'-- Year:',str(Start_Year),'-- SubTile:',str(SubTile_Index),'--',len(Left_Images_ID),'files rest!')
    
    else:
        print('Loop',i,'-- Year:',str(Start_Year),'-- SubTile:',str(SubTile_Index),'--',len(Left_Images_ID),'files rest!')
    
def FindFileName(RefList,FilePath):
    for root,dirs,files in os.walk(FilePath):
        for file in files:
            Image_ID = '_'.join((file.split('_'))[0:3])
            # 使用join函数将文件名称和文件所在根目录连接起来
            if file.endswith('.tif') and RefList.count(Image_ID)>0:
                RefList.remove(Image_ID)
    return RefList

from multiprocessing import Process
if __name__ == '__main__':

    Path_Dir = "K:\\BaiduSyncdisk\\Zhanglihao\\Shrub_Encroachment\\20240606_LandsatNDVI_Data\\"
    Years = [1986,1987,1988,1989,1990,1991,1992,1993,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,\
             2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020]

    # Start_Year = 1997
    # End_Year = Start_Year
    SubTile_Index = [0,1,2]

    process_list = []
    i=0
    for Year in Years:
        p = Process(target=Export_Image_toLocal,args=(Year,SubTile_Index,Path_Dir))
        i+=1
        p.start()
        process_list.append(p)

    for pro in process_list:
        pro.join()




