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
;+
; :description:
;    Describe the procedure.
;
;
; 利用GIMMS NDVI计算的物候数据SOS，上升阶段的POS95，下降阶段的POS，和EOS确定同一个像元的物候期相似的年份
;
;
; :author: 18702
;-
PRO Zlh01_get_similar_phenology_years_2pos95_sos_eos_gimms_1986_2015_old
  COMPILE_OPT IDL2

  Dir = 'K:\BaiduSyncdisk\Zhanglihao\Shrub_Encroachment\TP_AVHRR_1982_2015_Phenology\'

  Start_Time = SYSTIME(1)

  GOTO,Step_a

  Step_a:

  ENVI,/restore_base_save_files
  ENVI_BATCH_INIT

  ;读取影像
  SOS_FileName = Dir + 'qtp_avhrr_desnow_svg1_gdspara7_th_sos20'
  ENVI_OPEN_FILE,SOS_FileName,R_FID=SOS_ID
  ENVI_FILE_QUERY,SOS_ID,DIMS=SOS_DIMS,NB=SOS_NB,NS=SOS_NS,NL=SOS_NL
  Proj_info_SOS = ENVI_GET_PROJECTION(FID=SOS_ID)

  POS_L_FileName = Dir + 'qtp_avhrr_desnow_svg1_gdspara7_th_pos_l95'
  ENVI_OPEN_FILE,POS_L_FileName,R_FID=POS_L_ID
  ENVI_FILE_QUERY,POS_L_ID,DIMS=POS_L_DIMS,NB=POS_L_NB,NS=POS_L_NS,NL=POS_L_NL
  Proj_info_POS_L = ENVI_GET_PROJECTION(FID=POS_L_ID)

  POS_R_FileName = Dir + 'qtp_avhrr_desnow_svg1_gdspara7_th_pos_r95'
  ENVI_OPEN_FILE,POS_R_FileName,R_FID=POS_R_ID
  ENVI_FILE_QUERY,POS_R_ID,DIMS=POS_R_DIMS,NB=POS_R_NB,NS=POS_R_NS,NL=POS_R_NL
  Proj_info_POS_R = ENVI_GET_PROJECTION(FID=POS_R_ID)

  EOS_FileName = Dir + 'qtp_avhrr_desnow_svg1_gdspara7_th_eos50'
  ENVI_OPEN_FILE,EOS_FileName,R_FID=EOS_ID
  ENVI_FILE_QUERY,EOS_ID,DIMS=EOS_DIMS,NB=EOS_NB,NS=EOS_NS,NL=EOS_NL
  Proj_info_EOS = ENVI_GET_PROJECTION(FID=EOS_ID)

  Start_Year = 1982
  END_Year = 2015

  SOS_Data_1982_2015 = FLTARR(SOS_NS,SOS_NL,SOS_NB)+!VALUES.F_nan
  EOS_Data_1982_2015 = FLTARR(EOS_NS,EOS_NL,EOS_NB)+!VALUES.F_nan
  POS_L_Data_1982_2015 = FLTARR(POS_L_NS,POS_L_NL,POS_L_NB)+!VALUES.F_nan
  POS_R_Data_1982_2015 = FLTARR(POS_R_NS,POS_R_NL,POS_R_NB)+!VALUES.F_nan

  FOR i=0L,EOS_NB-1L DO BEGIN

    Tem_SOS = ENVI_GET_DATA(FID=SOS_ID,DIMS=SOS_DIMS,POS=i)
    Tem_POS_L = ENVI_GET_DATA(FID=POS_L_ID,DIMS=POS_L_DIMS,POS=i)
    Tem_POS_R = ENVI_GET_DATA(FID=POS_R_ID,DIMS=POS_R_DIMS,POS=i)
    Tem_EOS = ENVI_GET_DATA(FID=EOS_ID,DIMS=EOS_DIMS,POS=i)

    SOS_Data_1982_2015[*,*,i] = Tem_SOS
    POS_L_Data_1982_2015[*,*,i] = Tem_POS_L
    POS_R_Data_1982_2015[*,*,i] = Tem_POS_R
    EOS_Data_1982_2015[*,*,i] = Tem_EOS

  ENDFOR

  No_Valid_Index1 = WHERE(SOS_Data_1982_2015 LT 1,Count1)
  IF Count1 GT 0 THEN BEGIN
    SOS_Data_1982_2015[No_Valid_Index1] = !VALUES.F_nan
  ENDIF

  No_Valid_Index2 = WHERE(POS_L_Data_1982_2015 LT 1,Count2)
  IF Count2 GT 0 THEN BEGIN
    POS_L_Data_1982_2015[No_Valid_Index2] = !VALUES.F_nan
  ENDIF

  No_Valid_Index3 = WHERE(EOS_Data_1982_2015 LT 1,Count3)
  IF Count3 GT 0 THEN BEGIN
    EOS_Data_1982_2015[No_Valid_Index3] = !VALUES.F_nan
  ENDIF

  No_Valid_Index4 = WHERE(POS_R_Data_1982_2015 LT 1,Count4)
  IF Count4 GT 0 THEN BEGIN
    POS_R_Data_1982_2015[No_Valid_Index4] = !VALUES.F_nan
  ENDIF

  ;进行3倍标准差的剔除
  Mean_SOS_Data_1982_2015 = Mean(SOS_Data_1982_2015,DIMENSION=3,/NAN)
  Mean_POS_L_Data_1982_2015 = Mean(POS_L_Data_1982_2015,DIMENSION=3,/NAN)
  Mean_POS_R_Data_1982_2015 = Mean(POS_R_Data_1982_2015,DIMENSION=3,/NAN)
  Mean_EOS_Data_1982_2015 = Mean(EOS_Data_1982_2015,DIMENSION=3,/NAN)

  Stdev_SOS_Data_1982_2015 = Stddev(SOS_Data_1982_2015,DIMENSION=3,/NAN)
  Stdev_POS_L_Data_1982_2015 = Stddev(POS_L_Data_1982_2015,DIMENSION=3,/NAN)
  Stdev_POS_R_Data_1982_2015 = Stddev(POS_R_Data_1982_2015,DIMENSION=3,/NAN)
  Stdev_EOS_Data_1982_2015 = Stddev(EOS_Data_1982_2015,DIMENSION=3,/NAN)

  Up_3Q_SOS_Data_1982_2015 = Mean_SOS_Data_1982_2015 + Stdev_SOS_Data_1982_2015*3
  Up_3Q_POS_L_Data_1982_2015 = Mean_POS_L_Data_1982_2015 + Stdev_POS_L_Data_1982_2015*3
  Up_3Q_POS_R_Data_1982_2015 = Mean_POS_R_Data_1982_2015 + Stdev_POS_R_Data_1982_2015*3
  Up_3Q_EOS_Data_1982_2015 = Mean_EOS_Data_1982_2015 + Stdev_EOS_Data_1982_2015*3

  Down_3Q_SOS_Data_1982_2015 = Mean_SOS_Data_1982_2015 - Stdev_SOS_Data_1982_2015*3
  Down_3Q_POS_L_Data_1982_2015 = Mean_POS_L_Data_1982_2015 - Stdev_POS_L_Data_1982_2015*3
  Down_3Q_POS_R_Data_1982_2015 = Mean_POS_R_Data_1982_2015 - Stdev_POS_R_Data_1982_2015*3
  Down_3Q_EOS_Data_1982_2015 = Mean_EOS_Data_1982_2015 - Stdev_EOS_Data_1982_2015*3

  No_Valid_Index1 = WHERE((SOS_Data_1982_2015 LT Down_3Q_SOS_Data_1982_2015) OR (SOS_Data_1982_2015 GT Up_3Q_SOS_Data_1982_2015),Count1)
  IF Count1 GT 0 THEN BEGIN
    SOS_Data_1982_2015[No_Valid_Index1] = !VALUES.F_nan
  ENDIF

  No_Valid_Index2 = WHERE((POS_L_Data_1982_2015 LT Down_3Q_POS_L_Data_1982_2015) OR (POS_L_Data_1982_2015 GT Up_3Q_POS_L_Data_1982_2015),Count2)
  IF Count2 GT 0 THEN BEGIN
    POS_L_Data_1982_2015[No_Valid_Index2] = !VALUES.F_nan
  ENDIF

  No_Valid_Index3 = WHERE((EOS_Data_1982_2015 LT Down_3Q_EOS_Data_1982_2015) OR (EOS_Data_1982_2015 GT Up_3Q_EOS_Data_1982_2015),Count3)
  IF Count3 GT 0 THEN BEGIN
    EOS_Data_1982_2015[No_Valid_Index3] = !VALUES.F_nan
  ENDIF

  No_Valid_Index4 = WHERE((POS_R_Data_1982_2015 LT Down_3Q_POS_R_Data_1982_2015) OR (POS_R_Data_1982_2015 GT Up_3Q_POS_R_Data_1982_2015),Count4)
  IF Count4 GT 0 THEN BEGIN
    POS_R_Data_1982_2015[No_Valid_Index4] = !VALUES.F_nan
  ENDIF

  PRINT,Count1,Count2,Count3,Count4

  Save,SOS_Data_1982_2015,FILENAME=Dir+'SOS_Data_1982_2015.sav'
  Save,POS_L_Data_1982_2015,FILENAME=Dir+'POS_L_Data_1982_2015.sav'
  Save,POS_R_Data_1982_2015,FILENAME=Dir+'POS_R_Data_1982_2015.sav'
  Save,EOS_Data_1982_2015,FILENAME=Dir+'EOS_Data_1982_2015.sav'

  ENVI_BATCH_EXIT
  PRINT, SYSTIME(1)-Start_Time,'S'
  RETURN


  Step_b:
  ;将相似物候数据输出为tiff文件
  RESTORE,FILENAME=Dir+'SOS_Data_1982_2015.sav';SOS_Data_1982_2015
  RESTORE,FILENAME=Dir+'POS_L_Data_1982_2015.sav';POS_L_Data_1982_2015
  RESTORE,FILENAME=Dir+'POS_R_Data_1982_2015.sav';POS_R_Data_1982_2015
  RESTORE,FILENAME=Dir+'EOS_Data_1982_2015.sav';EOS_Data_1982_2015

  TP_TH20_SOS_2000_2020 = SOS_Data_1982_2015[*,*,4:*]
  TP_TH95_POSL_2000_2020 = POS_L_Data_1982_2015[*,*,4:*]
  TP_TH95_POSR_2000_2020 = POS_R_Data_1982_2015[*,*,4:*]
  TP_TH50_EOS_2000_2020 = EOS_Data_1982_2015[*,*,4:*]

  Dimensions = SIZE(TP_TH20_SOS_2000_2020,/DIMENSIONS)
  NB = Dimensions[-1]
  NS = Dimensions[0]
  NL = Dimensions[1]

  ;保存有有效值的区域
  Mask_FileName = Dir + 'Valid_Phenology_Mask.tif'
  Mask_Data = READ_TIFF(Mask_FileName,GEOTIFF=Geoproj)

  ;计算多年平均的SOS和多年平均的EOS
  ;  Mean_TP_TH20_SOS_2000_2020 = Mean(TP_TH20_SOS_2000_2020,DIMENSION=3,/NAN)
  ;  Mean_TP_TH50_EOS_2000_2020 = Mean(TP_TH50_EOS_2000_2020,DIMENSION=3,/NAN)
  ;  WRITE_TIFF,Dir+'TP_MODIS_Mean_TP_TH20_SOS_2000_2020.tif',Mean_TP_TH20_SOS_2000_2020,/FLOAT,GEOTIFF=Geoproj
  ;  WRITE_TIFF,Dir+'TP_MODIS_Mean_TP_TH50_EOS_2000_2020.tif',Mean_TP_TH50_EOS_2000_2020,/FLOAT,GEOTIFF=Geoproj

  ;计算25%分位数的SOS和75%分位数的EOS

  ;  RETURN
  Per25_TP_TH20_SOS_2000_2020 = FLTARR(NS,NL)+!VALUES.F_nan
  Per75_TP_TH50_EOS_2000_2020 = FLTARR(NS,NL)+!VALUES.F_nan

  SOS_Dis_Arr_MODIS = ULONARR(NB,NS,NL)
  POS_L_Dis_Arr_MODIS = ULONARR(NB,NS,NL)
  POS_R_Dis_Arr_MODIS = ULONARR(NB,NS,NL)
  EOS_Dis_Arr_MODIS = ULONARR(NB,NS,NL)
  Nor_Arr = INTARR(NB)+1

  SOS_DOY_Threshold = 15
  POS_DOY_Threshold = 15
  EOS_DOY_Threshold = 10

  FOR i=0L,NS-1L DO BEGIN
    FOR j=0L, NL-1L DO BEGIN

      Tem_SOS = REFORM(TP_TH20_SOS_2000_2020[i,j,*])
      ;  a#nor-nor#b
      Tem_Dis1 = Tem_SOS#Nor_Arr - Nor_Arr#Tem_SOS
      SOS_Dis_Arr_MODIS[*,i,j] = REFORM(Arr2d_tranform_le_32bit(ABS(Tem_Dis1),SOS_DOY_Threshold))

      Per25_TP_TH20_SOS_2000_2020[i,j] = Quantile(Tem_SOS,0.25,NAN=!VALUES.F_nan)

      Tem_POS_L = REFORM(TP_TH95_POSL_2000_2020[i,j,*])
      ;  a#nor-nor#b
      Tem_Dis2 = Tem_POS_L#Nor_Arr - Nor_Arr#Tem_POS_L
      POS_L_Dis_Arr_MODIS[*,i,j] = REFORM(Arr2d_tranform_le_32bit(ABS(Tem_Dis2),POS_DOY_Threshold))

      Tem_EOS = REFORM(TP_TH50_EOS_2000_2020[i,j,*])
      ;  a#nor-nor#b
      Tem_Dis3 = Tem_EOS#Nor_Arr - Nor_Arr#Tem_EOS
      EOS_Dis_Arr_MODIS[*,i,j] = REFORM(Arr2d_tranform_le_32bit(ABS(Tem_Dis3),EOS_DOY_Threshold))

      Per75_TP_TH50_EOS_2000_2020[i,j] = Quantile(Tem_EOS,0.75,NAN=!VALUES.F_nan)

      Tem_POS_R = REFORM(TP_TH95_POSR_2000_2020[i,j,*])
      ;  a#nor-nor#b
      Tem_Dis4 = Tem_POS_R#Nor_Arr - Nor_Arr#Tem_POS_R
      POS_R_Dis_Arr_MODIS[*,i,j] = REFORM(Arr2d_tranform_le_32bit(ABS(Tem_Dis4),POS_DOY_Threshold))

    ENDFOR
;    Progressbar,i,NS,Start_Time
  ENDFOR

  FOR i=0L,NB-1L DO BEGIN
    ;    TemData = ULONARR(4,NS,NL)

    ;    TemData[0,*,*] = reform(SOS_Dis_Arr_MODIS[i,*,*])
    ;    TemData[1,*,*] = REFORM(POS_L_Dis_Arr_MODIS[i,*,*])
    ;    TemData[2,*,*] = REFORM(POS_R_Dis_Arr_MODIS[i,*,*])
    ;    TemData[3,*,*] = REFORM(EOS_Dis_Arr_MODIS[i,*,*])

    ;用于检查代码
    ;    x=346
    ;    y=87
    ;
    ;    Tem_Arr_SOS = (Int_to_bit((reform(SOS_Dis_Arr_MODIS[i,x,y])),2,30))
    ;    Tem_Arr_POS_L = (Int_to_bit((REFORM(POS_L_Dis_Arr_MODIS[i,x,y])),2,30))
    ;    Tem_Arr_POS_R = (Int_to_bit((REFORM(POS_R_Dis_Arr_MODIS[i,x,y])),2,30))
    ;    Tem_Arr_EOS = (Int_to_bit((REFORM(EOS_Dis_Arr_MODIS[i,x,y])),2,30))
    ;
    ;    Valid_Index = WHERE((Tem_Arr_SOS EQ 1) AND (Tem_Arr_POS_L EQ 1)$
    ;      AND (Tem_Arr_POS_R EQ 1) AND (Tem_Arr_EOS EQ 1),Count_Valid,/L64)

    TemData0 = REFORM(SOS_Dis_Arr_MODIS[i,*,*]) AND REFORM(POS_L_Dis_Arr_MODIS[i,*,*]) AND REFORM(POS_R_Dis_Arr_MODIS[i,*,*])AND REFORM(EOS_Dis_Arr_MODIS[i,*,*])
    ;    Tem = (Int_to_bit((REFORM(TemData0[x,y])),2,30))
    ;    Tem_valid_index = where(Tem eq 1)

    ResultName = Dir+'GIMMS_1986_2015_SimilarPhenology\GIMMS_1986_2015_SimilarPhenology_'+(1986+i).Tostring()+'.tif'
    WRITE_TIFF,ResultName,TemData0,/LONG,GEOTIFF=Geoproj
    PRINT,'Finised!',ResultName
  ENDFOR

  RETURN

END