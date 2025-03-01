FUNCTION Arr_tranform_to_string,Input_Arr,Split=Split
  IF ~KEYWORD_SET(Split) THEN BEGIN
    Split = '_'
  ENDIF
  Str_arr = ''
  Count = SIZE(Input_Arr,/N_ELEMENTS)
  IF Count GT 0 THEN BEGIN
    FOR i=0L, Count-1L DO BEGIN
      IF i EQ 0 THEN BEGIN
        Str_arr += (Input_Arr[i]).Tostring()
      ENDIF ELSE BEGIN
        Str_arr += '_'+(Input_Arr[i]).Tostring()
      ENDELSE

    ENDFOR
  ENDIF

  RETURN,Str_arr
END
;2-d按竖列进行转换
FUNCTION Arr2d_tranform_le_32bit,Input_arr,Threshold
  Dimensions = SIZE(Input_arr,/DIMENSIONS)
  IF (N_ELEMENTS(Dimensions) NE 2) AND (Dimensions[0] NE Dimensions[1]) THEN BEGIN
    MESSAGE,'The Input array must be 2-D with the same Dimensions!',/IOERROR
  ENDIF

  IF ~KEYWORD_SET(Threshold) THEN BEGIN
    MESSAGE,'Please set the Threshold params',/IOERROR
  ENDIF

  Labels = ULONG(Input_arr LE Threshold)
  Dims = Dimensions[0]
  Tem_arr = (2UL)^ULINDGEN(Dims)
  Result = ULONG(Tem_arr##Labels)
  RETURN,Result
END


FUNCTION Arr1d_tranform_le_32bit,Input_arr,Threshold
  Dimensions = SIZE(Input_arr,/DIMENSIONS)
  IF N_ELEMENTS(Dimensions) NE 1 THEN BEGIN
    MESSAGE,'The Input array must be 1-D!',/IOERROR
  ENDIF
  IF ~KEYWORD_SET(Threshold) THEN BEGIN
    MESSAGE,'Please set the Threshold params',/IOERROR
  ENDIF

  Labels = ULONG(Input_arr LE Threshold)
  Dims = N_ELEMENTS(Input_arr)
  Tem_arr = (2UL)^ULINDGEN(Dims)
  Result = ULONG(TOTAL(Tem_arr*Labels))
  RETURN,Result
END

FUNCTION Int_to_bit,Input,bit,N_Dims

  ;  IF ~KEYWORD_SET(Input)and (Input ne 0) THEN BEGIN
  ;    MESSAGE,'Please set the Input params',/IOERROR
  ;  ENDIF

  IF ~KEYWORD_SET(bit) THEN BEGIN
    MESSAGE,'Please set the bit params',/IOERROR
  ENDIF

  IF ~KEYWORD_SET(N_Dims) THEN BEGIN
    MESSAGE,'Please set the N_Dims params',/IOERROR
  ENDIF

  result = INTARR(N_Dims)

  FOR i=0,N_Dims-1L DO BEGIN
    Tem = Input MOD bit
    result[i] = Tem
    Input = Input/bit
  ENDFOR

  RETURN,result
END
;+
; :description:
;    Describe the procedure.
;
;
; 利用MODIS NDVI计算的物候数据SOS，上升阶段的POS95，下降阶段的POS，和EOS确定同一个像元的物候期相似的年份
;
; :author: 18702
;-
PRO Zlh01_get_similar_phenology_years_2pos95_sos_eos_modis_2000_2020_old
  COMPILE_OPT IDL2

  ;  Dir = 'E:\Zhanglihao\Shrub_Encroachment\TP_MODIS_2000_2020_Phenology\'
  Dir = 'J:\Zhanglihao_Data\Landsat_NDVImax_Method\TP_MODIS_2000_2020_Phenology\'

  SOS_Dir = 'K:\fangbo_DISK_N\1NREE_250M_NDVIEVI\19_250M_MODIS_NDVI_SOS&EOS\THSOS\'
  EOS_Dir = 'K:\fangbo_DISK_N\1NREE_250M_NDVIEVI\19_250M_MODIS_NDVI_SOS&EOS\THEOS\'
  POS1_Dir = 'E:\fb\POS1_YEAR_TIF\'
  POS2_Dir = 'E:\fb\POS2_YEAR_TIF\'

  Start_Time = SYSTIME(1)

  GOTO,Step_a


  Step_a:


  ;  TemFileName = POS1_Dir + 'MODISNDVI_th95_POS1_2000.tif'
  ;  Tem_Data = READ_TIFF(TemFileName,GEOTIFF=GeoKeys)
  ;  Dimensions = SIZE(Tem_Data,/DIMENSIONS)
  ;  Bands = Dimensions[0]
  ;  Width = Dimensions[1]
  ;  Height = Dimensions[2]

  ENVI,/restore_base_save_files
  ENVI_BATCH_INIT

  ;读取影像
  StudyArea_FileName = POS1_Dir + 'MODISNDVI_th95_POS1_2000.tif'
  ENVI_OPEN_FILE,StudyArea_FileName,R_FID=StudyArea_ID
  ENVI_FILE_QUERY,StudyArea_ID,DIMS=StudyArea_DIMS,NB=StudyArea_NB,NS=StudyArea_NS,NL=StudyArea_NL
  Proj_info_StudyArea = ENVI_GET_PROJECTION(FID=StudyArea_ID)

  Start_Year = 2000
  End_Year = 2020
  N_Years = End_Year-Start_Year+1

  TP_TH20_SOS_2000_2020 = FLTARR(StudyArea_NS,StudyArea_NL,N_Years)+!VALUES.F_nan
  TP_TH50_EOS_2000_2020 = FLTARR(StudyArea_NS,StudyArea_NL,N_Years)+!VALUES.F_nan
  TP_TH95_POSL_2000_2020 = FLTARR(StudyArea_NS,StudyArea_NL,N_Years)+!VALUES.F_nan
  TP_TH95_POSR_2000_2020 = FLTARR(StudyArea_NS,StudyArea_NL,N_Years)+!VALUES.F_nan

  ;读取SOS数据
  FOR i=Start_Year,End_Year DO BEGIN
    Tem_Filename = SOS_Dir + 'MODISNDVI_th_sos_'+i.Tostring('(I4)')+'.tif'
    ENVI_OPEN_FILE,Tem_FileName,R_FID=Tem_ID
    IF Tem_ID LE 0 THEN BEGIN
      MESSAGE,'The input file not exist',Tem_Filename
    ENDIF
    Tem_Data = ENVI_GET_DATA(FID=Tem_ID,DIMS=StudyArea_DIMS,POS=0)
    TP_TH20_SOS_2000_2020[*,*,(i-Start_Year)] = REFORM(Tem_Data)
  ENDFOR

  PRINT,SYSTIME(1)-Start_Time,'Reading SOS Data Finished!'

  ;读取EOS数据
  FOR i=Start_Year,End_Year DO BEGIN
    Tem_Filename = EOS_Dir + 'MODISNDVI_th_eos_'+i.Tostring('(I4)')+'.tif'
    ENVI_OPEN_FILE,Tem_FileName,R_FID=Tem_ID
    IF Tem_ID LE 0 THEN BEGIN
      MESSAGE,'The input file not exist',Tem_Filename
    ENDIF
    Tem_Data = ENVI_GET_DATA(FID=Tem_ID,DIMS=StudyArea_DIMS,POS=0)
    TP_TH50_EOS_2000_2020[*,*,(i-Start_Year)] = REFORM(Tem_Data)
  ENDFOR

  PRINT,SYSTIME(1)-Start_Time,'Reading EOS Data Finished!'

  ;读取POS1数据
  FOR i=Start_Year,End_Year DO BEGIN
    Tem_Filename = POS1_Dir + 'MODISNDVI_th95_POS1_'+i.Tostring('(I4)')+'.tif'
    ENVI_OPEN_FILE,Tem_FileName,R_FID=Tem_ID
    IF Tem_ID LE 0 THEN BEGIN
      MESSAGE,'The input file not exist',Tem_Filename
    ENDIF
    Tem_Data = ENVI_GET_DATA(FID=Tem_ID,DIMS=StudyArea_DIMS,POS=0)
    TP_TH95_POSL_2000_2020[*,*,(i-Start_Year)] = REFORM(Tem_Data)
  ENDFOR

  PRINT,SYSTIME(1)-Start_Time,'Reading POS1 Data Finished!'

  ;读取POS2数据
  FOR i=Start_Year,End_Year DO BEGIN
    Tem_Filename = POS2_Dir + 'MODISNDVI_th95_POS2_'+i.Tostring('(I4)')+'.tif'
    ENVI_OPEN_FILE,Tem_FileName,R_FID=Tem_ID
    IF Tem_ID LE 0 THEN BEGIN
      MESSAGE,'The input file not exist',Tem_Filename
    ENDIF
    Tem_Data = ENVI_GET_DATA(FID=Tem_ID,DIMS=StudyArea_DIMS,POS=0)
    TP_TH95_POSR_2000_2020[*,*,(i-Start_Year)] = REFORM(Tem_Data)
  ENDFOR

  PRINT,SYSTIME(1)-Start_Time,'Reading POS2 Data Finished!'

  No_Valid_Index1 = WHERE(TP_TH20_SOS_2000_2020 LT 1,Count1,/L64)
  IF Count1 GT 0 THEN BEGIN
    TP_TH20_SOS_2000_2020[No_Valid_Index1] = !VALUES.F_nan
  ENDIF

  No_Valid_Index2 = WHERE(TP_TH95_POSL_2000_2020 LT 1,Count2,/L64)
  IF Count2 GT 0 THEN BEGIN
    TP_TH95_POSL_2000_2020[No_Valid_Index2] = !VALUES.F_nan
  ENDIF

  No_Valid_Index3 = WHERE(TP_TH50_EOS_2000_2020 LT 1,Count3,/L64)
  IF Count3 GT 0 THEN BEGIN
    TP_TH50_EOS_2000_2020[No_Valid_Index3] = !VALUES.F_nan
  ENDIF

  No_Valid_Index4 = WHERE(TP_TH95_POSR_2000_2020 LT 1,Count4,/L64)
  IF Count4 GT 0 THEN BEGIN
    TP_TH95_POSR_2000_2020[No_Valid_Index4] = !VALUES.F_nan
  ENDIF

  PRINT,SYSTIME(1)-Start_Time
  PRINT,Count1,Count2,Count3,Count4

  ;进行3倍标准差的剔除
  Mean_TP_TH20_SOS_2000_2020 = Mean(TP_TH20_SOS_2000_2020,DIMENSION=3,/NAN)
  Mean_TP_TH95_POSL_2000_2020 = Mean(TP_TH95_POSL_2000_2020,DIMENSION=3,/NAN)
  Mean_TP_TH95_POSR_2000_2020 = Mean(TP_TH95_POSR_2000_2020,DIMENSION=3,/NAN)
  Mean_TP_TH50_EOS_2000_2020 = Mean(TP_TH50_EOS_2000_2020,DIMENSION=3,/NAN)

  Stdev_TP_TH20_SOS_2000_2020 = Stddev(TP_TH20_SOS_2000_2020,DIMENSION=3,/NAN)
  Stdev_TP_TH95_POSL_2000_2020 = Stddev(TP_TH95_POSL_2000_2020,DIMENSION=3,/NAN)
  Stdev_TP_TH95_POSR_2000_2020 = Stddev(TP_TH95_POSR_2000_2020,DIMENSION=3,/NAN)
  Stdev_TP_TH50_EOS_2000_2020 = Stddev(TP_TH50_EOS_2000_2020,DIMENSION=3,/NAN)

  Up_3Q_TP_TH20_SOS_2000_2020 = Mean_TP_TH20_SOS_2000_2020 + Stdev_TP_TH20_SOS_2000_2020*3
  Up_3Q_TP_TH95_POSL_2000_2020 = Mean_TP_TH95_POSL_2000_2020 + Stdev_TP_TH95_POSL_2000_2020*3
  Up_3Q_TP_TH95_POSR_2000_2020 = Mean_TP_TH95_POSR_2000_2020 + Stdev_TP_TH95_POSR_2000_2020*3
  Up_3Q_TP_TH50_EOS_2000_2020 = Mean_TP_TH50_EOS_2000_2020 + Stdev_TP_TH50_EOS_2000_2020*3

  Down_3Q_TP_TH20_SOS_2000_2020 = Mean_TP_TH20_SOS_2000_2020 - Stdev_TP_TH20_SOS_2000_2020*3
  Down_3Q_TP_TH95_POSL_2000_2020 = Mean_TP_TH95_POSL_2000_2020 - Stdev_TP_TH95_POSL_2000_2020*3
  Down_3Q_TP_TH95_POSR_2000_2020 = Mean_TP_TH95_POSR_2000_2020 - Stdev_TP_TH95_POSR_2000_2020*3
  Down_3Q_TP_TH50_EOS_2000_2020 = Mean_TP_TH50_EOS_2000_2020 - Stdev_TP_TH50_EOS_2000_2020*3

  No_Valid_Index1 = WHERE((TP_TH20_SOS_2000_2020 LT Down_3Q_TP_TH20_SOS_2000_2020) OR (TP_TH20_SOS_2000_2020 GT Up_3Q_TP_TH20_SOS_2000_2020),Count1,/L64)
  IF Count1 GT 0 THEN BEGIN
    TP_TH20_SOS_2000_2020[No_Valid_Index1] = !VALUES.F_nan
  ENDIF

  No_Valid_Index2 = WHERE((TP_TH95_POSL_2000_2020 LT Down_3Q_TP_TH95_POSL_2000_2020) OR (TP_TH95_POSL_2000_2020 GT Up_3Q_TP_TH95_POSL_2000_2020),Count2,/L64)
  IF Count2 GT 0 THEN BEGIN
    TP_TH95_POSL_2000_2020[No_Valid_Index2] = !VALUES.F_nan
  ENDIF

  No_Valid_Index3 = WHERE((TP_TH50_EOS_2000_2020 LT Down_3Q_TP_TH50_EOS_2000_2020) OR (TP_TH50_EOS_2000_2020 GT Up_3Q_TP_TH50_EOS_2000_2020),Count3,/L64)
  IF Count3 GT 0 THEN BEGIN
    TP_TH50_EOS_2000_2020[No_Valid_Index3] = !VALUES.F_nan
  ENDIF

  No_Valid_Index4 = WHERE((TP_TH95_POSR_2000_2020 LT Down_3Q_TP_TH95_POSR_2000_2020) OR (TP_TH95_POSR_2000_2020 GT Up_3Q_TP_TH95_POSR_2000_2020),Count4,/L64)
  IF Count4 GT 0 THEN BEGIN
    TP_TH95_POSR_2000_2020[No_Valid_Index4] = !VALUES.F_nan
  ENDIF

  PRINT,Count1,Count2,Count3,Count4

  ;只保留同时存在SOS，POS以及EOS的像元
  No_Valid_Index = WHERE((FINITE(TP_TH20_SOS_2000_2020,/NAN)) OR (FINITE(TP_TH95_POSL_2000_2020,/NAN)) $
    OR (FINITE(TP_TH95_POSR_2000_2020,/NAN)) OR (FINITE(TP_TH50_EOS_2000_2020,/NAN)),Count,/L64)

  IF Count GT 0 THEN BEGIN
    TP_TH20_SOS_2000_2020[No_Valid_Index] = !VALUES.F_nan
    TP_TH95_POSL_2000_2020[No_Valid_Index] = !VALUES.F_nan
    TP_TH95_POSR_2000_2020[No_Valid_Index] = !VALUES.F_nan
    TP_TH50_EOS_2000_2020[No_Valid_Index] = !VALUES.F_nan
  ENDIF
  PRINT,Count

  Save,TP_TH20_SOS_2000_2020,FILENAME=Dir+'TP_TH20_SOS_2000_2020.sav'
  Save,TP_TH95_POSL_2000_2020,FILENAME=Dir+'TP_TH95_POSL_2000_2020.sav'
  Save,TP_TH95_POSR_2000_2020,FILENAME=Dir+'TP_TH95_POSR_2000_2020.sav'
  Save,TP_TH50_EOS_2000_2020,FILENAME=Dir+'TP_TH50_EOS_2000_2020.sav'

  ENVI_BATCH_EXIT
  PRINT,SYSTIME(1)-Start_Time,'S'
  RETURN

  ;逐像元算一个距离矩阵

  Step_b:
  RESTORE,FILENAME=Dir+'TP_TH20_SOS_2000_2020.sav';TP_TH20_SOS_2000_2020
  RESTORE,FILENAME=Dir+'TP_TH95_POSL_2000_2020.sav';TP_TH95_POSL_2000_2020
  RESTORE,FILENAME=Dir+'TP_TH95_POSR_2000_2020.sav';TP_TH95_POSR_2000_2020
  RESTORE,FILENAME=Dir+'TP_TH50_EOS_2000_2020.sav';TP_TH50_EOS_2000_2020

  Dimensions = SIZE(TP_TH20_SOS_2000_2020,/DIMENSIONS)
  NB = Dimensions[-1]
  NS = Dimensions[0]
  NL = Dimensions[1]

  ;保存有有效值的区域
  ;  Dir0 = 'K:\BaiduSyncdisk\Zhanglihao\Shrub_Encroachment\TP_AVHRR_1982_2015_Phenology\'
  ;  TemFileName = POS1_Dir + 'MODISNDVI_th95_POS1_2000.tif'
  ;  Tem_Data = READ_TIFF(TemFileName,GEOTIFF=Geoproj)

  Tem_Data = REFORM(TOTAL(TP_TH20_SOS_2000_2020 GT 0,3,/NAN)) GT 0
  WRITE_TIFF,Dir+'MODIS_Valid_Phenology_Mask_2000_2020.tif',Tem_Data,/LONG,GEOTIFF=Geoproj

  Per25_TP_TH20_SOS_2000_2020 = FLTARR(NS,NL)+!VALUES.F_nan
  Per75_TP_TH50_EOS_2000_2020 = FLTARR(NS,NL)+!VALUES.F_nan

  SOS_Dis_Arr_MODIS = ULONARR(NS,NL,NB)
  POS_L_Dis_Arr_MODIS = ULONARR(NS,NL,NB)
  POS_R_Dis_Arr_MODIS = ULONARR(NS,NL,NB)
  EOS_Dis_Arr_MODIS = ULONARR(NS,NL,NB)
  Nor_Arr = INTARR(NB)+1

  SOS_DOY_Threshold = 15
  POS_DOY_Threshold = 15
  EOS_DOY_Threshold = 10

  FOR i=0L,NS-1L DO BEGIN
    FOR j=0L, NL-1L DO BEGIN

      Tem_SOS = REFORM(TP_TH20_SOS_2000_2020[i,j,*])
      ;  a#nor-nor#b
      Tem_Dis1 = Tem_SOS#Nor_Arr - Nor_Arr#Tem_SOS
      SOS_Dis_Arr_MODIS[i,j,*] = REFORM(Arr2d_tranform_le_32bit(ABS(Tem_Dis1),SOS_DOY_Threshold))

      Per25_TP_TH20_SOS_2000_2020[i,j] = Quantile(Tem_SOS,0.25,NAN=!VALUES.F_nan)

      Tem_POS_L = REFORM(TP_TH95_POSL_2000_2020[i,j,*])
      ;  a#nor-nor#b
      Tem_Dis2 = Tem_POS_L#Nor_Arr - Nor_Arr#Tem_POS_L
      POS_L_Dis_Arr_MODIS[i,j,*] = REFORM(Arr2d_tranform_le_32bit(ABS(Tem_Dis2),POS_DOY_Threshold))

      Tem_EOS = REFORM(TP_TH50_EOS_2000_2020[i,j,*])
      ;  a#nor-nor#b
      Tem_Dis3 = Tem_EOS#Nor_Arr - Nor_Arr#Tem_EOS
      EOS_Dis_Arr_MODIS[i,j,*] = REFORM(Arr2d_tranform_le_32bit(ABS(Tem_Dis3),EOS_DOY_Threshold))

      Per75_TP_TH50_EOS_2000_2020[i,j] = Quantile(Tem_EOS,0.75,NAN=!VALUES.F_nan)

      Tem_POS_R = REFORM(TP_TH95_POSR_2000_2020[i,j,*])
      ;  a#nor-nor#b
      Tem_Dis4 = Tem_POS_R#Nor_Arr - Nor_Arr#Tem_POS_R
      POS_R_Dis_Arr_MODIS[i,j,*] = REFORM(Arr2d_tranform_le_32bit(ABS(Tem_Dis4),POS_DOY_Threshold))

    ENDFOR
    ;    Progressbar,i,NS,Start_Time
  ENDFOR

  Save,SOS_Dis_Arr_MODIS,FILENAME=Dir+'SOS_Dis_Arr_MODIS.sav'
  Save,POS_L_Dis_Arr_MODIS,FILENAME=Dir+'POS_L_Dis_Arr_MODIS.sav'
  Save,POS_R_Dis_Arr_MODIS,FILENAME=Dir+'POS_R_Dis_Arr_MODIS.sav'
  Save,EOS_Dis_Arr_MODIS,FILENAME=Dir+'EOS_Dis_Arr_MODIS.sav'


  WRITE_TIFF,Dir+'TP_MODIS_Per25_TP_TH20_SOS_2000_2020.tif',Per25_TP_TH20_SOS_2000_2020,/FLOAT,GEOTIFF=Geoproj
  WRITE_TIFF,Dir+'TP_MODIS_Per75_TP_TH50_EOS_2000_2020.tif',Per75_TP_TH50_EOS_2000_2020,/FLOAT,GEOTIFF=Geoproj


  PRINT, SYSTIME(1)-Start_Time,'S'
  RETURN

  Step_c:
  ;将相似年份数据输出
  RESTORE,FILENAME=Dir+'SOS_Dis_Arr_MODIS.sav'
  RESTORE,FILENAME=Dir+'POS_L_Dis_Arr_MODIS.sav'
  RESTORE,FILENAME=Dir+'POS_R_Dis_Arr_MODIS.sav'
  RESTORE,FILENAME=Dir+'EOS_Dis_Arr_MODIS.sav'

  Mean_SOS_FileName = Dir + 'TP_MODIS_Per25_TP_TH20_SOS_2000_2020.tif'
  Mean_SOS_Data = READ_TIFF(Mean_SOS_FileName,GEOTIFF=Geoproj)

  Dimensions = SIZE(EOS_Dis_Arr_MODIS,/DIMENSION)
  NB = Dimensions[-1]
  NS = Dimensions[0]
  NL = Dimensions[1]
  ;逐年的输出影像
  FOR i=0L,NB-1L DO BEGIN
    TemData0 = REFORM(SOS_Dis_Arr_MODIS[*,*,i]) AND REFORM(POS_L_Dis_Arr_MODIS[*,*,i]) $
      AND REFORM(POS_R_Dis_Arr_MODIS[*,*,i])AND REFORM(EOS_Dis_Arr_MODIS[*,*,i])

    ResultName = Dir+'MODIS_2000_2020_SimilarPhenology\MODIS_2000_2020_SimilarPhenology_'+(2000+i).Tostring()+'.tif'
    WRITE_TIFF,ResultName,TemData0,/LONG,GEOTIFF=Geoproj
    PRINT,'Finised!',ResultName
  ENDFOR

  RETURN
END