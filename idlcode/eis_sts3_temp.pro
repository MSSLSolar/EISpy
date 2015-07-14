;+
; NAME:
;       EIS_STS3_TEMP()
;
;
; PURPOSE:
;       Return the house keeping temperatures recorded closest in time to
;       a specified raster image exposure index.
;
;
; CATEGORY:
;       Wavelength correction of the orbital variation of the EIS line centre.
;       Correction based on house keeping temperatures and S. Kamio's
;       artificial neural network approach
;
;
; CALLING SEQUENCE:
;       Called by EIS_MODEL_SERIES()
;
;
; INPUTS:
;       TIMES_EIS3: the times of each EIS house keeping temperature measurement 
;       EIS_STS3: EIS house keeping temperatures, downloaded from Hinode
;       Science Data Centre Europe (SDC) in Oslo or present in the SSW distribution
;       POS_EIS3: return temperatures at this index (and other indices that
;                 are calculated from POS_EIS3) in the TIMES_EIS3 array.
;
; OPTIONAL INPUTS:
;       POS_EIS3: If not provided, the keyword TIME must be set.
;
;
; KEYWORD PARAMETERS:
;       TIME: find the house keeping temperatures that are recorded closest in
;             time to TIME, instead of providing the index of this time in the
;             POS_EIS3 input variable.
;
;
; OUTPUTS:
;       Array with temperatures.
;
;
; OPTIONAL OUTPUTS:
;       None
;
;
; SPECIAL CALLS:
;       FIND_CLOSE()
;
;
; COMMON BLOCKS:
;       None
;
;
; SIDE EFFECTS:
;       None
;
;
; RESTRICTIONS:
;       None
;       
;
; EXAMPLE:
;       See EIS_MODEL_SERIES()
;
;
; MODIFICATION HISTORY:
;       xx-xx-2009: Suguro Kamio
;       10-03-2010: Suguro Kamio - New version
;       08-06-2010: Terje Fredvik - Doc header
;-

FUNCTION eis_sts3_temp, times_eis3, eis_sts3, pos_eis3, time=time

temp = fltarr(34)

IF keyword_set(time) THEN BEGIN
   pos_eis3 = find_close(times_eis3, time)
ENDIF

IF (pos_eis3 LE 0) OR (pos_eis3 GE (n_elements(times_eis3) - 1)) THEN BEGIN
   print,'eis_sts3_temp: EIS STS3 is not avaliable at given time.'
   return,0
ENDIF

temp[0] = eis_sts3[pos_eis3].temp[1]; eis_mhc_chassis_t1
temp[1] = eis_sts3[pos_eis3].temp[5]; eis_mhc_sla_t5
temp[2] = eis_sts3[pos_eis3].temp[7]; eis_mhc_mir_base_t7
temp[3] = eis_sts3[pos_eis3].temp[10];eis_mhc_gra_mtr_t10
temp[4] = eis_sts3[pos_eis3].temp[11];eis_mhc_gra_t11
temp[5] = eis_sts3[pos_eis3].temp[13]; eis_mhc_hz_t16
temp[6] = eis_sts3[pos_eis3].temp[14]; eis_mhc_hz_t1
temp[7] = eis_sts3[pos_eis3].temp[15]; eis_mhc_hz_t2
temp[8] = eis_sts3[pos_eis3].temp[16]; eis_mhc_hz_t3
temp[9] = eis_sts3[pos_eis3].temp[17]; eis_mhc_hz_t4
temp[10] = eis_sts3[pos_eis3].temp[21];eis_mhc_hz_t8
temp[11] = eis_sts3[pos_eis3].temp[22];eis_mhc_hz_t9
temp[12] = eis_sts3[pos_eis3].temp[23];eis_mhc_hz_t10
temp[13] = eis_sts3[pos_eis3].temp[24];eis_mhc_hz_t11
temp[14] = eis_sts3[pos_eis3].temp[26];eis_mhc_hz_t13
temp[15] = eis_sts3[pos_eis3].temp[28];eis_mhc_hz_t15
j = (pos_eis3 - 5) > 0
temp[16] = eis_sts3[j].temp[1]; eis_mhc_chassis_t1
temp[17] = eis_sts3[j].temp[10];eis_mhc_gra_mtr_t10
temp[18] = eis_sts3[j].temp[11];eis_mhc_gra_t11
temp[19] = eis_sts3[j].temp[13];eis_mhc_hz_t16
temp[20] = eis_sts3[j].temp[14];eis_mhc_hz_t1
temp[21] = eis_sts3[j].temp[22];eis_mhc_hz_t9
temp[22] = eis_sts3[j].temp[23];eis_mhc_hz_t10
temp[23] = eis_sts3[j].temp[24];eis_mhc_hz_t11
temp[24] = eis_sts3[j].temp[26];eis_mhc_hz_t13
last = n_elements(times_eis3) - 1
j = (pos_eis3 + 5) < last
temp[25] = eis_sts3[j].temp[1]; eis_mhc_chassis_t1
temp[26] = eis_sts3[j].temp[10];eis_mhc_gra_mtr_t10
temp[27] = eis_sts3[j].temp[11];eis_mhc_gra_t11
temp[28] = eis_sts3[j].temp[13];eis_mhc_hz_t16
temp[29] = eis_sts3[j].temp[14];eis_mhc_hz_t1
temp[30] = eis_sts3[j].temp[22];eis_mhc_hz_t9
temp[31] = eis_sts3[j].temp[23];eis_mhc_hz_t10
temp[32] = eis_sts3[j].temp[24];eis_mhc_hz_t11
temp[33] = eis_sts3[j].temp[26];eis_mhc_hz_t13

return,temp

END 
