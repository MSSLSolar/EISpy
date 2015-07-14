;+
; NAME:
;        EIS_TEMP_MODEL()
;
;
; PURPOSE:
;       Apply the input house keeping temperatures to a set of coefficients
;       found by S. Kamio's artificial neural network approach.
;
;
;
; CATEGORY:
;       Wavelength correction of the orbital variation of the line centre.
;       Correction based on house keeping temperatures and S. Kamio's neural
;       network approach
;
;
; CALLING SEQUENCE:
;       Called by EIS_MODEL_SERIES()
;
;
; INPUTS:
;       TEMP: a set of temperatures recorded at a specific time; i.e. the output
;       from EIS_STS3_TEMP(). 
;
;
;
; OPTIONAL INPUTS:
;       None
;
;
; KEYWORD PARAMETERS:
;       TIME: the time of the scientific observation. Different sets of
;             coefficients may be applied at different dates.
;       QUIET: if set, do not print a warning to the terminal when the TIME
;              keyword is not set.
;       SLIT:  if set to 2 it is assumed that the scientific observation was
;             obtained by the 2" slit.
;
;
; OUTPUTS:
;       The wavelength shift of a single exposure in a raster image, measured
;       in CCD pixel coordinates
;
;
; OPTIONAL OUTPUTS:
;       None
;
;
; SPECIAL CALLS:
;       None      
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
;       10-03-2010: Sugoro Kamio - New version
;       08-06-2010: Terje Fredvik - Doc header
;
;-

function eis_temp_model,temp,time=time,quiet=quiet,slit=slit
;;  
;; New version 10.03.2010
;;  
if keyword_set(time) eq 0 then begin
if keyword_set(quiet) eq 0 then print,'eis_temp_model: Time is not specified. Assuming current time.'
get_utc,utc
time = anytim(utc)
endif

slit2_offset = -8.20D ; offset between two slits (pixel)

t1 = 9.1229760D+08; = anytim('2007-11-29T00:00:00') heater config
t2 = 9.3553920D+08; = anytim('2008-08-24T00:00:00') slit focus adjust
t3 = 9.4057920D+08; = anytim('2008-10-21T08:00:00') grating focus adjust

if (time lt t1) then begin
; Hinode launch -- heater config on 2007-11-29


pixel_ref = 1.34524D+03

c = [$
4.10562D-01,$
2.51204D+00,$
-7.03979D-01,$
1.21183D+00,$
-1.46165D+00,$
-2.03801D+00,$
-5.09189D+00,$
-3.31613D+00,$
2.28654D-01,$
3.72455D+00,$
8.19741D-01,$
1.17212D+00,$
3.19226D+00,$
2.21462D+00,$
-2.76307D+00,$
-7.75230D+00,$
2.27707D+00,$
8.62746D-02,$
-3.87772D+00,$
8.50736D-01,$
2.50457D-01,$
-4.62109D+00,$
-1.49986D+00,$
-9.98911D-01,$
-5.24012D+00,$
-4.88090D+00,$
8.41629D-01,$
1.53231D+00,$
-5.56888D+00,$
5.46359D+00,$
5.00476D+00,$
6.83911D+00,$
2.10491D+00,$
6.89056D+00 ]

endif

if (time ge t1) and (time lt t2) then begin
; heater config on 2007-11-29 -- slit focus on 2008-08-24
pixel_ref = 1.34915D+03

c = [$
-7.60169D+00,$
-1.46383D+00,$
3.64224D+00,$
6.22838D+00,$
1.02071D+00,$
-5.87856D+00,$
-7.07813D+00,$
-3.29145D+00,$
-2.68002D+00,$
6.44214D+00,$
-5.64250D+00,$
9.41400D+00,$
1.02490D+01,$
1.00514D+00,$
1.54987D+01,$
-2.43897D+01,$
6.93774D+00,$
7.99804D+00,$
-4.24839D+00,$
1.94191D+00,$
-4.11472D+00,$
2.67682D+00,$
2.63193D+00,$
-1.58034D+00,$
-1.36976D+01,$
-1.78314D+00,$
-3.97698D+00,$
-5.86437D+00,$
2.30465D+00,$
1.23473D+01,$
-1.35947D+00,$
1.85987D+00,$
4.27904D+00,$
-4.35809D+00 ]

endif

if time ge t2 then begin
; slit focus on 2008-08-24 -- present

pixel_ref = 1.34281D+03

if (time ge t2) and (time lt t3) then pixel_ref += 4.88D
; offset slit focus on 2008-08-24 -- grating focus on 2008-10-21
c = [$
-9.69118D-01,$
2.12159D+00,$
-2.99428D+00,$
2.61100D+00,$
1.41035D+00,$
-9.76397D-01,$
-1.61651D+01,$
-9.94312D-01,$
1.04603D+00,$
8.57033D-01,$
2.07951D+00,$
4.80522D+00,$
8.65133D+00,$
-2.37848D-02,$
1.09901D+00,$
-5.51204D+00,$
1.58325D+00,$
1.97708D+00,$
-3.42620D+00,$
1.76606D+00,$
6.50817D+00,$
-7.19983D+00,$
-3.21551D+00,$
-6.81840D-01,$
-5.75801D+00,$
-1.08458D-01,$
-3.76701D+00,$
-3.05294D+00,$
-4.01884D+00,$
1.00570D+01,$
4.61089D-01,$
6.69429D+00,$
-6.84122D-01,$
4.38880D+00 ]

endif

if keyword_set(slit) eq 0 then slit = 0
if slit eq 2 then pixel_ref = pixel_ref + slit2_offset

return, total(c * (temp - 15.0D) / 10.0D) + pixel_ref
end
