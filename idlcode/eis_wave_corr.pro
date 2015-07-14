
;+
; 
;             NAME : eis_wave_corr
; 
;          PURPOSE : EIS wavelength corrections for orbital variations
;
; CALLING SEQUENCE : eis_wave_corr,file
; 
;           INPUTS : file - the full EIS level 1 file name.
; 
;  OPTIONAL INPUTS : none
; 
;         KEYWORDS : Some optional inputs can be passed to eis_wave_corr_old, see that
;                    routine for information.
; 
;          OUTPUTS : wvl_corr - the wavelength correction for each pixel
;
;                  : dw_tilt - the slit tilt component
;
;                  : dw_t - the temporal component
;
;                  : also saves these arrays to a genx file for future
;                    use.
;
; OPTIONAL OUTPUTS : none
;
;    COMMON BLOCKS : none
; 
;          EXAMPLE : IDL> eis_wave_corr,file,wvl_corr,dw_tilt,dw_t
;                    IDL> d.wvl = d.wvl - wvl_corr[i,j]
;
;    SPECIAL CALLS : eis_wave_corr_hk - the real work gets done here!
;
;          WRITTEN :  HPW : 22-JUN-2010 : based on earlier code
;                     HPW : 25-SEP-2010 : added call to eis_wave_corr_hk & /old
; 
;-

pro eis_wave_corr,file,wvl_corr,dw_tilt,dw_t,$
                  old=old,$
                  ymin=ymin,ymax=ymax,$
                  show=show,$
                  nospline=nospline,$
                  do_fits=do_fits

  ;; --- call the old way if necessary
  if keyword_set(old)  then begin
    eis_wave_corr_old,file,wvl_corr,dw_tilt,dw_t,ymin=ymin,ymax=ymax,$
        show=show,$
        nospline=nospline,$
        do_fits=do_fits
    return
  endif

  ;; --- wavelength correction is computed here
  eis_wave_corr_hk,File,wvl_corr,dw_tilt,dw_t

  ;; --- save the output
  break_file,file,disk,dir,name,ext
  opf = str_replace(name,'_l0_','_l1_')
  opf = opf + '.wave'

  method = 'HouseKeeping Data'
  ymin   = -1
  ymax   = -1
  
  text = [" ---------- EIS_WAVE_CORR OUTPUT ---------- ",$
          " computed with "+method,$
          " ymin (pixels) = "+trim(ymin),$
          " ymax (pixels) = "+trim(ymax),$
          " file = "+opf+".genx",$
          " read with ",$
          "  IDL> restgen,wvl_corr,dw_tilt,dw_t,file='"+opf+".genx'",$
          "  wvl_corr = 2D wavelength correction array",$
          "  dw_t = 1D time correction",$
          "  dw_tilt = 1D slit tilt correction",$
          " use with ",$
          "  d.wvl = d.wvl - wvl_corr[i,j]",$
          " ----------------------------------------- "]

  savegen,wvl_corr,dw_tilt,dw_t,file=opf+'.genx',text=text

return
end
