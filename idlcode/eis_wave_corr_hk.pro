;+
; NAME:
;        EIS_WAVE_CORR_HK()
;
;
; PURPOSE: 
;       Estimate the orbital variations of the EIS line centres using
;       Kamio-san's artificial neural network code and house keeping 
;       temperatures, combine with a correction of the slit tilt (which is 
;       calculated without the use of house keeping temperatures).
;
;
; CATEGORY:
;       Wavelength correction of the orbital variation of the line centre.
;       Correction based on house keeping temperatures and S. Kamio's neural
;       network approach
;
;
; CALLING SEQUENCE:
;       eis_wave_corr_hk,fileOrObj,wvl_corr[,dw_tilt,dw_t,lam=lam,wvl_cube=wcvl_cube]
;
;
; INPUTS: 
;       fileOrObj: EIS level 0 or level 1 file or an eis_data object
;       (level 0 or level1).
;
;
; OPTIONAL INPUTS:
;       None
;
;
; KEYWORD PARAMETERS:
;       WVL_CUBE : if set to a named variable. pass back a 3D wavelength cube
;                 which is the corrected 1D wavelength scale at each
;                 pixel. wvl_cube is in the format used by the cfit line
;                 fitting package. If set, the LAM keyword must also be set if
;                 the line of interest is not Fe XII 195.12 A.
;       LAM : the wavelength for which the wavelength correction is to be
;             calculated for. We assume that the wavelength correction is the
;             same for all lines on both detectors, when measured in pixels,
;             but the correction will vary slightly with wavelength when
;             measured in Angstroms. However, the variation is very small,
;             e.g. is a typical maximum difference between the correction for
;             Fe XII 195.12 AA and Fe XV 284.16 AA approximately 0.1 km/s.
;
;
; OUTPUTS:
;       WVL_CORR : The (total) wavelength correction (in Angstrom) for each pixel in
;                  the raster image (temporal variations + slit tilt) 
;
;
; OPTIONAL OUTPUTS:
;      W_TILT : If set to a named variable, pass back a 1D array containing the
;               slit tilt component of the total wavelength correction (Angstrom)
;      DW_T :   If set to a named variable, pass back a 1D array containing the
;               orbital variation (temporal variation) component of the total 
;               wavelength correction (Angstrom)
;
; SPECIAL CALLS:
;       An EIS_DATA object is created and methods of this object are called,
;       the methods of this object will call EIS_MODEL_SERIES(),
;       EIS_STS3_TEMP(), EIS_TEMP_MODEL(), FIND_CLOSE() and FPP1_DOPP_SERIES().
;
;
; COMMON BLOCKS:
;       None
;
;
; SIDE EFFECTS:
;       If EIS house keeping data is not already stored on the local hard
;       disk, files will be downloaded from the Oslo archive and saved locally
;
;
; RESTRICTIONS:
;       If the EIS house keeping data is not already stored on the local hard disk
;       an Internet connection is required in order to download the data from
;       the Oslo archive. 
;
;
; EXAMPLE:
;       IDL> eis_wave_corr_hk,fileOrObj, wvl_corr, dw_tilt, dw_t,lamcube=lamcube
;       IDL> ana = mk_analysis(lamcube, data,wts,adef,miss,r,res)
;
;
; MODIFICATION HISTORY:
;       15-04-2010: Terje Fredvik
;       08-06-2010: Doc header
;-

PRO eis_wave_corr_hk, fileOrObj, wvl_corr, dw_tilt, dw_t, $
                      wvl_cube=wvl_cube, lam=lam

  IF n_params() EQ 0 THEN $
     message,'Usage: eis_wave_corr_hk, fileOrObj, wvl_corr [,dw_tilt] [,dw_t] [,wvl_cube=wvl_cube] [,lam=lam]'
  
  default, lam, 195.12
  getcube = arg_present(wvl_cube)
   
  obj = ((size(fileOrObj,/tname)) EQ 'OBJREF') ? 1 : 0
  IF obj THEN BEGIN 
     IF obj_class(fileOrObj) EQ 'EIS_DATA' THEN o = fileOrObj ELSE $
        message,'Input object must be of class EIS_DATA.'
  ENDIF ELSE o = obj_new('eis_data',fileOrObj)
  wavecorr = o->gethkwavecorr(lam,wvl_cube=getcube)
  
  wvl_corr = wavecorr.im
  
  IF arg_present(dw_tilt) THEN dw_tilt = wavecorr.tilt
  IF arg_present(dw_t) THEN dw_t = wavecorr.time
    
  IF getcube THEN wvl_cube = wavecorr.cube
  
  IF ~obj THEN obj_destroy,o

END
