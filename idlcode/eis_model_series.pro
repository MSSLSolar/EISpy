;+
; NAME:
;        EIS_MODEL_SERIES()
;
;
; PURPOSE: 
;        Model the wavelength shift (due to the orbital variation of
;        the instrument temperatures) of all exposures in a EIS raster image,
;        measured in CCD pixel coordinates.
;
; CATEGORY: 
;        Wavelength correction of the orbital variation of the line centre.
;        Correction based on house keeping temperatures and S. Kamio's neural
;        network approach.
;
;
; CALLING SEQUENCE:
;        EIS_MODEL_SERIES may be called by EIS_DATA::SETHKPIXCORR when 1) the 
;        EIS_WAVE_CORR_HK wrapper is called, 2) when EIS_DATA::GETHKWAVECORR 
;        is called or when 3) EIS_PREP is called with the HKWAVECORR keyword set.
;
;
; INPUTS:
;       TIMES_EIS3: the times of each EIS house keeping temperature measurement 
;       EIS_STS3: EIS house keeping temperatures, downloaded from Hinode
;       Science Data Centre Europe (SDC) in Oslo.
;       TIMES: the time of each exposure of the EIS scientific data raster.
;
;
; OPTIONAL INPUTS:
;       None
;
;
; KEYWORD PARAMETERS:
;       None
;
;
; OUTPUTS:
;       The wavelength shift of all exposures in the raster, measured in 
;       CCD pixel coordinates
;
;
; OPTIONAL OUTPUTS:
;       None
;
; SPECIAL CALLS:
;       FIND_CLOSE(), EIS_STS3_TEMP(), EIS_TEMP_MODEL() 
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
;      See EIS_DATA::SETHKPIXCORR
;
;
; MODIFICATION HISTORY:
;      xx-xx-2009: Suguro Kamio
;      07-06-2010: Terje Fredvik - Doc header
;      16-06-2010: Terje Fredvik - Added keywords GOODEXP and BADEXP
;-

function eis_model_series,times_eis3, eis_sts3, times, slit=slit, $
                          goodexp=goodexp, badexp=badexp

IF NOT arg_present(goodexp) THEN goodexp = where(times gt 0)
if goodexp[0] eq -1 then begin
 print,'eis_model_series: no valid time'
 return,fltarr(n_elements(times))
endif
pos1 = find_close(times_eis3, min(times[goodexp]), /earlier)
pos2 = find_close(times_eis3, max(times[goodexp]), /later)
pix = fltarr(pos2 - pos1 + 1)

for i=pos1,pos2 do begin
 temp = eis_sts3_temp(times_eis3, eis_sts3, i)
 pix[i - pos1] = eis_temp_model(temp, time=times_eis3[i], slit=slit)
endfor

model = interpol(pix, times_eis3[pos1:pos2], times)
IF NOT arg_present(badexp) THEN badexp = where(times eq 0)
if badexp[0] ne -1 then model[badexp] = 0

return,model

end
