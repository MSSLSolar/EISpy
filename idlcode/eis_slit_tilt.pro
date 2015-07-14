;+
; 
;             NAME : eis_slit_tilt
; 
;          PURPOSE : Calculates the slit tilt as a function position
;                    along the CCD.
;
; CALLING SEQUENCE : dw = eis_slit_tilt(yws,ny,long=long,short=short)
; 
;           INPUTS : yws - y window start, this parameter is in the
;                    header.
; 
;                    ny - the number of elements along y.
; 
;  OPTIONAL INPUTS : d - the structure returned from eis_getwindata.
; 
;         KEYWORDS : long  - slit tilt for long band
;                    short - slit tilt for short band
;                    slit  - slit tilte for 1" (slit=1) or 2" (slit=2)
;                    date  - date of observations, account for focus adjustment 
; 
;          OUTPUTS : An estimate of the slit tilt in Angstroms.
;
; OPTIONAL OUTPUTS : locations - CCD Y positions (for plotting)
;                    info - information on the calculation
;
;    COMMON BLOCKS : none
; 
;          EXAMPLE : dw  = eis_slit_tilt(256,512,/short,locations=y)
;                    wvl = wvl - dw
;                    plot,y,dw
;                    ALTERNATIVE USAGE:
;                    d  = eis_getwindata(file,wave)
;                    dw = eis_slit_tilt(d,locations=y)
;
;    SPECIAL CALLS : Needs to have the file eis_slit_title.genx in the
;                    same directory.
;
;          WRITTEN : HPW : 17-MAY-2007
;                    HPW : 23-MAY-2007 : added long/short
;                    HPW : 25-JUN-2007 : modified for new calculations
;                    HPW : 07-JUL-2007 : long/short bug fixes
;                    HPW : 09-OCT-2007 : added 2" slit
;                    HPW : 27-JAN-2010 : modified for new coefficents
;                    HPW : 28-JUL-2014 : stopped slit input from being modified
; 
;-

function eis_slit_tilt,input,ny,long=long,short=short,locations=locations,$
                       info=info,slit=input_slit,date=date

  ;; --- check inputs
  if datatype(input) eq 'STC' then begin
    ny    = input.ny
    yws   = input.hdr[0].yws
    det   = strcompress(input.hdr[0].TWBND,/remove_all)
    short = 0 
    long  = 0
    slit  = str_replace(input.hdr[0].slit_id,'"')
    if det eq 'A' then long = 1 else short = 1
    date  = input.time_ccsds[0]
  endif else begin
    yws = input
    if not(keyword_set(input_slit)) then begin
      message,'SLIT NOT SPECIFIED, RETURING RESULTS FOR 1" SLIT',/info
      slit = '1'
    endif else slit = strmid(trim(input_slit),0,1)
    if not(keyword_set(date)) then begin
      now  = str2arr(systime(0),' ')
      date = now[2]+'-'+now[1]+'-'+now[4]
      message,'DATE NOT SPECIFIED, RETURING RESULTS FOR TODAY',/info
      message,' '+date,/info
    endif
      
  endelse

  case 1 of 
    keyword_set(short): i = 0
    keyword_set(long):  i = 1
    else: message,"EITHER LONG OR SHORT KEYWORD MUST BE USED"
  endcase

  ;; --- find the slit tilt file
  file = find_with_def('eis_slit_tilt.pro',!path)
  file = str_replace(file,'.pro','.txt')

  ;; --- read the fit
  src = rd_tfile(file,/auto,/convert)

  ;; --- find the data for this date
  FOCUS_DATE = '24-AUG-2008 00:00'
  if anytim(date) lt anytim(FOCUS_DATE) then begin
    dateS = 'Before Focus'
    fit   = src[*,0:3]
  endif else begin
    dateS = 'After Focus'
    fit = src[*,4:7]
  endelse

  ;; --- find the data for this wavelength
  if keyword_set(short) then begin
    waveS = 'SW'
    fit   = fit[*,0:1] 
  endif else begin
    waveS = 'LW'
    fit   = fit[*,2:3]
  endelse

  ;; --- find the data for this slit
  if slit eq '1' then begin
    slitS = '1"'
    fit   = reform(fit[*,1]) 
  endif else begin
    slitS = '2"'
    fit   = reform(fit[*,0])
  endelse

  ;; -- calcultate the tilt in Angstroms as a function of pixel position
  y  = yws + findgen(ny)
  dw = poly(y,fit)
  
  dispersion = 0.0223
  dw = dw*dispersion

  ;; -- some optional output
  locations = y
  bands = ['SW','LW']
  slit  = ['1"','2"']
  info  = {date: dateS, fit: fit, $
           yws: yws, ny:  ny, band: waveS,slit: slitS}

return,dw
end
