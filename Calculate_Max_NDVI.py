import pandas as pd 
import os
import json
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from sklearn.preprocessing import PolynomialFeatures,SplineTransformer
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import Ridge
from scipy.interpolate import CubicSpline
from scipy.interpolate import BSpline, make_interp_spline
from scipy import optimize as op
from scipy.stats import linregress
# import seaborn as sns
from matplotlib.pyplot import MultipleLocator
from scipy.optimize import least_squares
from tqdm import tqdm,trange
import math
import gc
from multiprocessing import Process,Pool
from osgeo import gdal


def Tifread(Tiff_Filename,Nodata_Value=np.nan):
    DataSet = gdal.Open(Tiff_Filename)
    Band = DataSet.GetRasterBand(1)
    NoDataValue = Band.GetNoDataValue()

    if DataSet == None:
        raise Exception(f"cannot find/open {Tiff_Filename}")
    Image_Width = DataSet.RasterXSize
    Image_Height = DataSet.RasterYSize
    Image_Geotrans = DataSet.GetGeoTransform()
    Image_Proj = DataSet.GetProjection()
    Image_Bandnum = DataSet.RasterCount
    Image_Data = DataSet.ReadAsArray()
    # Image_Data[Image_Data==Nodata_Value] = np.nan
    # Image_Data[np.isnan(Image_Data)] = np.nan
 
    del DataSet
    return Image_Data#, Image_Geotrans, Image_Proj

# f = m1+(m2-m7*t)*(1.0/(1.0+(exp((m3-t)/m4)))-1.0/(1.0+(exp((m5-t)/m6))))
def Double_logistic_7(x,*a):

    m1 = a[0]
    m2 = a[1]
    m3 = a[2]
    m4 = a[3]
    m5 = a[4]
    m6 = a[5]
    m7 = a[6]
    
    # ret = m1 + (m2)*(1.0/(1.0+np.exp((m3-x)/m4))-1.0/(1.0+np.exp(m5-x)/m6))
    ret = m1 + (m2-m7*x)*(1.0/(1.0+np.exp((m3-x)/m4))-1.0/(1.0+np.exp((m5-x)/m6)))

    return ret

def Double_logistic_7_residuals(params,X,Y):
    ret = Double_logistic_7(X,*params) - Y
    return ret


def Double_logistic_7_SMF3(x,x_offset,y_scale,y_offset,*a):
    
    m1 = a[0]
    m2 = a[1]
    m3 = a[2]
    m4 = a[3]
    m5 = a[4]
    m6 = a[5]
    m7 = a[6]
    
    # ret = m1 + (m2)*(1.0/(1.0+np.exp((m3-x)/m4))-1.0/(1.0+np.exp(m5-x)/m6))
    #保证变换前后背景值没有发生变化
    ret = m1+(m2-m7*(x+x_offset))*(1.0/(1.0+np.exp((m3-(x+x_offset))/m4))-1.0/(1.0+np.exp((m5-(x+x_offset))/m6)))#+m1
    ret = y_scale*ret+y_offset
    return ret


def int_to_bit(Input,Bit,N_Dims):
    result = np.zeros(N_Dims)
    for i in range(N_Dims):
        Tem = Input % Bit
        result[i] = Tem
        Input = Input//Bit
    return result


def Double_logistic_7_SMF3_residuals3(SMF_Params,X,Y,*a):#,X_Multi,Multi_Up,Multi_Down
    # x_scale = SMF_Params[0]
    x_offset = SMF_Params[0]
    y_scale = SMF_Params[1]
    y_offset = SMF_Params[2]

    # TemY_Down = Double_logistic_7_SMF3(X_Multi,x_offset,y_scale,y_offset,*a)-Multi_Down
    # TemY_Up = Double_logistic_7_SMF3(X_Multi,x_offset,y_scale,y_offset,*a)-Multi_Up

    ret = (Y-Double_logistic_7_SMF3(X,x_offset,y_scale,y_offset,*a))

    return ret


    
    return int(Mon_Index[0]+1)

def MultiProcess_Get_NDVI_Max(Pool_Num,SubTile_Index):
    #需要将每个Tile划分成不同的subtile,否则数据量太大，无法存储与计算
    Dir = 'K:\\BaiduSyncdisk\\Zhanglihao\\Shrub_Encroachment\\20240606_LandsatNDVI_Data\\'
    # SubTile_Index = 1
    SubTile_Path = 'SubTile_' + str(SubTile_Index) + '\\'

    #读取Image_Index文件，划分成100份
    Samples_File = Dir + SubTile_Path + 'SubTile_' + str(SubTile_Index) + '_Image_Index_Split.csv'
    Samples_Data = pd.read_csv(Samples_File)
    min_Samples_IDs = Samples_Data['tem_min'].values
    max_Samples_IDs = Samples_Data['tem_max'].values

    pool = Pool(Pool_Num)
    i=0
    for i in range(len(min_Samples_IDs)):
        min_Samples_ID = min_Samples_IDs[i]
        max_Samples_ID = max_Samples_IDs[i]
        pool.apply_async(func=Get_NDVI_Max_Tif,args=(SubTile_Index,min_Samples_ID,max_Samples_ID))

    pool.close()
    gc.collect()
    pool.join()


def Get_NDVI_Max_Tif(SubTile_Index,min_Samples_ID,max_Samples_ID):

    Dir = 'K:\\BaiduSyncdisk\\Zhanglihao\\Shrub_Encroachment\\20240606_LandsatNDVI_Data\\'
    # SubTile_Index = 1
    SubTile_Path = 'SubTile_' + str(SubTile_Index) + '\\'

    #先判断结果文件是否已经存在，不对已经处理过的数据进行处理
    ResultName = Dir + SubTile_Path +'ResultData\\SubTile_'+str(SubTile_Index)+'_'+str(min_Samples_ID)+'_'+str(max_Samples_ID)+'_Result.json'
    if os.path.isfile(ResultName):
        return
    
    Image_Index_Filename = Dir + SubTile_Path +'PreparedData\\NDVI_Json_Data\\SubTile_'+str(SubTile_Index)+'_'+str(min_Samples_ID)+'_'+str(max_Samples_ID)+'_Image_Index.tif'
    DOY_Filename = Dir + SubTile_Path +'PreparedData\\NDVI_Json_Data\\SubTile_'+str(SubTile_Index)+'_'+str(min_Samples_ID)+'_'+str(max_Samples_ID)+'_DOY.tif'
    NDVI_Filename = Dir + SubTile_Path +'PreparedData\\NDVI_Json_Data\\SubTile_'+str(SubTile_Index)+'_'+str(min_Samples_ID)+'_'+str(max_Samples_ID)+'_NDVI.tif'
    Mon_Filename = Dir + SubTile_Path +'PreparedData\\NDVI_Json_Data\\SubTile_'+str(SubTile_Index)+'_'+str(min_Samples_ID)+'_'+str(max_Samples_ID)+'_Mon.tif'
    Year_Filename = Dir + SubTile_Path +'PreparedData\\NDVI_Json_Data\\SubTile_'+str(SubTile_Index)+'_'+str(min_Samples_ID)+'_'+str(max_Samples_ID)+'_Year.tif'

    Image_Index_Data = Tifread(Image_Index_Filename)
    DOY_Data = Tifread(DOY_Filename)
    #将NDVI转换成0-1之间
    NDVI_Data_Origion = Tifread(NDVI_Filename)
    NDVI_Data = NDVI_Data_Origion/10000.0
    Mon_Data = Tifread(Mon_Filename)
    Year_Data = Tifread(Year_Filename)

    NDVI_Samples_DataFrame = pd.DataFrame({'ID':np.array(Image_Index_Data).flatten(),'NDVI_Data':np.array(NDVI_Data).flatten(), \
                                          'DOY':np.array(DOY_Data).flatten(),'Year':np.array(Year_Data).flatten(),'Month':np.array(Mon_Data).flatten()})
    NDVI_Samples_DataFrame.dropna(inplace=True)

    del Image_Index_Data
    del DOY_Data
    del NDVI_Data_Origion
    del NDVI_Data
    del Mon_Data
    del Year_Data
    
    # NDVI_Samples_DataFrame_Origion = NDVI_Samples_DataFrame.copy()

    # 读取物候期相似的年份数据
    GIMMS_Similar_Phenology_File = Dir + SubTile_Path +'PreparedData\\GIMMS_Similar_Phenology\\SubTile_'+str(SubTile_Index)+ '_Similar_Phenology_Data_1986_2015_'+str(min_Samples_ID)+'_'+str(max_Samples_ID)+'.tif'
    GIMMS_Similar_Phenology_Data = Tifread(GIMMS_Similar_Phenology_File)#Tifread(GIMMS_Similar_Phenology_File)
    MODIS_Similar_Phenology_File = Dir + SubTile_Path +'PreparedData\\MODIS_Similar_Phenology\\SubTile_'+str(SubTile_Index)+ '_Similar_Phenology_Data_2000_2020_'+str(min_Samples_ID)+'_'+str(max_Samples_ID)+'.tif'
    MODIS_Similar_Phenology_Data = Tifread(MODIS_Similar_Phenology_File)

    # Samples_ID_Arr = np.arange(min_Samples_ID,max_Samples_ID+1)
    Samples_ID_Arr_Real = NDVI_Samples_DataFrame['ID'].unique()

    #For Test
    # Samples_ID_Arr_Real = np.array([9725])
    # NDVI_Samples_DataFrame = NDVI_Samples_DataFrame.loc[NDVI_Samples_DataFrame['ID'].isin(Samples_ID_Arr_Real)].sort_values(by=['ID','Year','DOY'])
    # NDVI_Samples_DataFrame.to_csv(path_or_buf=Dir + SubTile_Path +'PreparedData\\SubTile_'+str(SubTile_Index)+'_'+str(min_Samples_ID)+'_'+str(max_Samples_ID)+'_Test.csv',index=False)

    Years_Arr_All = np.array([1986,1987,1988,1989,1990,1991,1992,1993,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,\
    2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020])

    Result_List = list()

    ValidCounts = 24

    Thelta = 3 #异常值的去除，均值+Thelta*标准差
    # for Iter_Num in Iter_ID_Arr:

    for Samples_ID in tqdm(Samples_ID_Arr_Real):
        Samples_ID_Index = int(Samples_ID-min_Samples_ID)
        # print(int(Samples_ID))
        GIMMS_Similar_Phenology_DataFrame = GIMMS_Similar_Phenology_Data[:,Samples_ID_Index]
        MODIS_Similar_Phenology_DataFrame = MODIS_Similar_Phenology_Data[:,Samples_ID_Index]

        Sample_Similar_Phenology_DataFrame = np.append(GIMMS_Similar_Phenology_DataFrame.flatten(),MODIS_Similar_Phenology_DataFrame.flatten(),axis=0)

        Has_Similar_Phenology = np.count_nonzero(Sample_Similar_Phenology_DataFrame)
        Is_Valid_Image_Index = np.sum(Samples_ID_Arr_Real==Samples_ID)

        if (Has_Similar_Phenology>0)&(Is_Valid_Image_Index>0):#

            NDVI_Sample_DataFrame = (NDVI_Samples_DataFrame.loc[(NDVI_Samples_DataFrame['ID']==Samples_ID)])
            Valid_Count_All = len(NDVI_Sample_DataFrame)

            if Valid_Count_All>0:
                Years = NDVI_Sample_DataFrame['Year'].unique()
                # Years = np.array([2012])
                # NDVI_Sample_DataFrame_Origion_Simulated_ID = (NDVI_Samples_DataFrame_Origion.loc[(NDVI_Samples_DataFrame_Origion['ID']==int(Samples_ID))]).sort_values(by="Year",ascending=True)
                for Year in Years:

                    NDVI_DOY_DataFrame_Year = NDVI_Sample_DataFrame.loc[NDVI_Sample_DataFrame['Year']==Year]

                    Valid_Count_Year = len(NDVI_DOY_DataFrame_Year)
                    Similar_Phenology_Year_Arr = []

                    if Valid_Count_Year>=4:
                        Max_NDVI_Origion = max(NDVI_DOY_DataFrame_Year['NDVI_Data'].values)
                        # 确定具有相似物候的年份
                        Year_Index = np.nonzero(Years_Arr_All==Year)
                        Tem_Similar_Value = (Sample_Similar_Phenology_DataFrame[Year_Index])[0]            
                        if Year<2000:
                            Tem_Start_Year = 1986
                            Tem_END_Year = 2015
                        else:
                            Tem_Start_Year = 2000
                            Tem_END_Year = 2020
                        Years_Arr = np.array(range(Tem_Start_Year,Tem_END_Year+1,1))
                        N_Dims = len(Years_Arr)
                        Similar_Phenology_Year_ArrIndex = np.nonzero(int_to_bit(Tem_Similar_Value,2,N_Dims))
                        Similar_Phenology_Year_Arr = Years_Arr[Similar_Phenology_Year_ArrIndex]                            
                    else:
                        Max_NDVI_Origion = np.nan
                    
                    if (len(Similar_Phenology_Year_Arr) > 0)&(Valid_Count_Year>=4):

                        # Similar_Phenology_Year_Arr包含自身的年份，所以查找其他年份时，不应该包含自身年份。
                        Similar_Phenology_Dataframe = (NDVI_Sample_DataFrame.loc[(NDVI_Sample_DataFrame['Year'].isin(Similar_Phenology_Year_Arr))]).copy()
                        # 将自身年份的模拟数据加上；这里这么繁琐是因为模拟的数据中每个样本点的每年的数据是原始数据的子集，不能直接将当年原始数据直接包含进去。
                        Valid_Count_SimilarYears = len(Similar_Phenology_Dataframe)

                        # 5-9月每月至少都有1个观测点
                        Sorted_Month_Dataframe = Similar_Phenology_Dataframe.loc[(Similar_Phenology_Dataframe['Month']>=5)&(Similar_Phenology_Dataframe['Month']<=9)]
                        DOY_Dis_Count = len(Sorted_Month_Dataframe['Month'].unique())

                        # 目标年份在夏季6,7,8月有1一个观测点
                        Sorted_Month_Dataframe0 = NDVI_DOY_DataFrame_Year.loc[(NDVI_DOY_DataFrame_Year['Month']>=6)&(NDVI_DOY_DataFrame_Year['Month']<=8)]
                        DOY_Dis_Count0 = len(Sorted_Month_Dataframe0['Month'].unique())

                        if (Valid_Count_SimilarYears >=ValidCounts) and (Valid_Count_Year >=4) and (DOY_Dis_Count ==5)and (DOY_Dis_Count0>0):

                            #计算冬季12，1,2月的75%-95%分位数的均值作为背景值。
                            # print(Year)
                            Year_Winter_Overlay_Sample_DataFrame = Similar_Phenology_Dataframe.loc[(Similar_Phenology_Dataframe['Month'].isin([1,2,12]))]
                            Background_Value = 0.1
                            if len(Year_Winter_Overlay_Sample_DataFrame)>0:
                                Quantile_75 = Year_Winter_Overlay_Sample_DataFrame['NDVI_Data'].quantile(0.5)
                                Quantile_95 = Year_Winter_Overlay_Sample_DataFrame['NDVI_Data'].quantile(0.95)
                                Tem_DataFrame = Year_Winter_Overlay_Sample_DataFrame.loc[(Year_Winter_Overlay_Sample_DataFrame['NDVI_Data']>=Quantile_75) \
                                                                                            & (Year_Winter_Overlay_Sample_DataFrame['NDVI_Data']<=Quantile_95)]
                                if len(Tem_DataFrame)>0:
                                    Background_Value = Tem_DataFrame['NDVI_Data'].mean()
                                else:
                                    Background_Value = 0.1
                            else:
                                Background_Value = 0.1
                            
                            #6-8月90%分位数以上的均值
                            Summer_Value = 0.1
                            SimilarYear_Summer_Sample_DataFrame = Similar_Phenology_Dataframe.loc[(Similar_Phenology_Dataframe['Month'].isin([6,7,8]))]
                            if len(SimilarYear_Summer_Sample_DataFrame)>0:
                                Quantile_90 = SimilarYear_Summer_Sample_DataFrame['NDVI_Data'].quantile(0.9)
                                Tem_DataFrame_Summer = SimilarYear_Summer_Sample_DataFrame.loc[(SimilarYear_Summer_Sample_DataFrame['NDVI_Data']>=Quantile_90)]
                                if len(Tem_DataFrame_Summer)>0:
                                    Summer_Value = Tem_DataFrame_Summer['NDVI_Data'].mean()
                                else:
                                    Summer_Value = 0.1
                            else:
                                Summer_Value = 0.1

                            if Summer_Value>=(1.2*Background_Value):
                                DOY_Data_SimilarYears = Similar_Phenology_Dataframe['DOY'].values
                                NDVI_Data_SimilarYears = Similar_Phenology_Dataframe['NDVI_Data'].values

                                #设定拟合参数的初始值
                                tem_min = max([Background_Value,0.001])
                                tem_max = max([Summer_Value,0.15])

                                tem_max0 = tem_max-tem_min

                                Tem_DOY_Data_SimilarYears = DOY_Data_SimilarYears#np.concatenate((DOY_Data_SimilarYears, Right_DOY))
                                Tem_NDVI_Data_SimilarYears = NDVI_Data_SimilarYears#np.concatenate((NDVI_Data_SimilarYears, Right_NDVI))

                                Initial_Value = [tem_min,tem_max0,150,20,270,20,0.00001]
                                Down_Bounds = [tem_min*0.5,tem_max0*0.5,90,15,210,15,-0.001]
                                Up_Bounds = [tem_min*1.5,tem_max0*2,210,30,330,30,0.001]

                                Least_result = least_squares(Double_logistic_7_residuals,Initial_Value,args=(Tem_DOY_Data_SimilarYears,Tem_NDVI_Data_SimilarYears), \
                                                            x_scale=[1,1,0.001,0.01,0.001,0.01,100],bounds=(Down_Bounds,Up_Bounds),max_nfev=50,method='trf',ftol=0.0001,xtol=0.01)
                                popt = Least_result.x
                                Orgion_DOY_Y = Double_logistic_7(DOY_Data_SimilarYears, *popt)

                                #先初步拟合，剔除异常点
                                Diff_NDVI = abs(NDVI_Data_SimilarYears - Orgion_DOY_Y)
                                Stdd_Diff_NDVI = np.std(Diff_NDVI)
                                Threshold_Diff = Thelta*Stdd_Diff_NDVI

                                #异常值的索引
                                # Nonzero_Index, = np.nonzero(Diff_NDVI>Threshold_Diff)
                                Zero_Index, = np.nonzero(Diff_NDVI<=Threshold_Diff)

                                New_DOY_Data_SimilarYears = DOY_Data_SimilarYears[Zero_Index]
                                New_NDVI_Data_SimilarYears = NDVI_Data_SimilarYears[Zero_Index]

                                if len(New_DOY_Data_SimilarYears)>=ValidCounts:

                                    New_DOY_X = NDVI_DOY_DataFrame_Year['DOY'].values
                                    New_DOY_Y = NDVI_DOY_DataFrame_Year['NDVI_Data'].values

                                    Initial_Value = popt
                                    # popt, pcov = op.curve_fit(Double_logistic_7, Tem_New_DOY_Data_SimilarYears, Tem_New_NDVI_Data_SimilarYears,Initial_Value,bounds=(Down_Bounds,Up_Bounds),max_nfev=1000)
                                    Least_result = least_squares(Double_logistic_7_residuals,Initial_Value,args=(Tem_DOY_Data_SimilarYears,Tem_NDVI_Data_SimilarYears), \
                                                                x_scale=[1,1,0.001,0.01,0.001,0.01,100],bounds=(Down_Bounds,Up_Bounds),max_nfev=10,method='trf',ftol=0.0001,xtol=0.01)
                                    popt = Least_result.x

                                    # #for test Plot data
                                    # plt.plot(Tem_DOY_Data_SimilarYears,Tem_NDVI_Data_SimilarYears,'og')
                                    # test_x = np.arange(1,366,1)
                                    # test_y = Double_logistic_7(test_x, *popt)
                                    # plt.plot(test_x,test_y,'-r')
                                    # plt.show()

                                    #4参数
                                    if (len(New_DOY_X)>=4):
                                        x0 = (np.array([0,1,0])).flatten()
                                        Down_Bounds = [-30,0.5,-0.5]
                                        Up_Bounds = [30,1.5,0.5]

                                        SMF_Fitting = least_squares(Double_logistic_7_SMF3_residuals3,x0,args=(New_DOY_X,New_DOY_Y,*popt),bounds=(Down_Bounds,Up_Bounds), \
                                                                    x_scale=[1,0.1,0.01],max_nfev=50,method='trf',ftol=0.0001,xtol=0.01)
                                        # SMF_Fitting = minimize(Double_logistic_7_SMF3_residuals3,x0,args=(New_DOY_X,New_DOY_Y,*popt),bounds=([-30,30],[0.5,1.5],[-0.1,0.1]))
                                        SMF_Coef = SMF_Fitting.x
                                        #获取拟合曲线的最大NDVI及其对应的DOY_m
                                        DOY_X = np.arange(1,366,1)
                                        x_offset = SMF_Coef[0]
                                        y_scale = SMF_Coef[1]
                                        y_offset = SMF_Coef[2]
                                    
                                        Fitted_NDVI = Double_logistic_7_SMF3(DOY_X,x_offset,y_scale,y_offset,*popt)
                                        Max_NDVI = max(Fitted_NDVI)
                                        if (Max_NDVI < 1.0):
                                            # x_scale = SMF_Coef[0]
                                            New_Fitted_NDVI = Double_logistic_7_SMF3(New_DOY_X,x_offset,y_scale,y_offset,*popt)
                                           
                                            #for test
                                            # plt.close()
                                            # plt.plot(New_DOY_X,New_DOY_Y,'og')
                                            # plt.plot(DOY_X,Fitted_NDVI,'-r')
                                            # plt.show()
                                            # plt.close()
                                            #获取校正偏差，数组类型
                                            Max_NDVI_Logistic_GIMMS = NDVI_DOY_DataFrame_Year['NDVI_Data'].values-(New_Fitted_NDVI-Max_NDVI)

                                            #小于Max origion NDVI的设为Max origion NDVI；大于等于1的去掉，然后再取中位数。
                                            Max_NDVI_Origion = max(NDVI_DOY_DataFrame_Year['NDVI_Data'].values)
                                            Max_NDVI_Logistic_GIMMS[Max_NDVI_Logistic_GIMMS < Max_NDVI_Origion] = Max_NDVI_Origion
                                            Max_NDVI_Logistic_GIMMS[Max_NDVI_Logistic_GIMMS >= 1] = np.nan
                                            Tem_Result = Max_NDVI_Logistic_GIMMS[~np.isnan(Max_NDVI_Logistic_GIMMS)]
                                            if len(Tem_Result)>0:
                                                SP_SMF = np.median(Tem_Result)
                                                m1 = popt[0]
                                                m2 = popt[1]
                                                m3 = popt[2]
                                                m4 = popt[3]
                                                m5 = popt[4]
                                                m6 = popt[5]
                                                m7 = popt[6]
                                                Tem_List = {'ID':Samples_ID,'Year':Year,'Max_Origion_NDVI':Max_NDVI_Origion,'SP_SMF':SP_SMF, \
                                                            'm1':m1,'m2':m2,'m3':m3,'m4':m4,'m5':m5,'m6':m6,'m7':m7,'x_offset':x_offset, \
                                                            'y_offset':y_offset,'y_scale':y_scale}
                                                Result_List.append(Tem_List)

            
    # 将结果输出
    Simulated_Samples_Overlay_DataFrame_Result = pd.DataFrame(Result_List)
    del Result_List
    Simulated_Samples_Overlay_DataFrame_Result.dropna(inplace=True)
    print(Simulated_Samples_Overlay_DataFrame_Result.shape)
    Out_Dir = Dir + SubTile_Path +'ResultData\\'
    if not os.path.exists(Out_Dir):
        os.makedirs(Out_Dir)
        print("创建文件夹："+Out_Dir)
    ResultName = Dir + SubTile_Path +'ResultData\\SubTile_'+str(SubTile_Index)+'_'+str(min_Samples_ID)+'_'+str(max_Samples_ID)+'_Result.json'
    Simulated_Samples_Overlay_DataFrame_Result.to_json(path_or_buf=ResultName,orient='records')


    del Simulated_Samples_Overlay_DataFrame_Result
    gc.collect()

