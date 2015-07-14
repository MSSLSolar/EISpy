;+
; NAME:
;        FPP1_DOPP_SERIES()
;
;
; PURPOSE:
;        Return the line-of-sight velocity of the spacecraft 
;
;
; CATEGORY:
;        Wavelength correction of the orbital variation of the EIS line centre
;
;
; CALLING SEQUENCE:
;        Called by EIS_DATA::SETHKPIXCORR
;
;
; INPUTS:
;        TIMES_FPP1: the times of each SOT/FPP measurement
;        PFF1: The velocity of the spacecraft in the line of sight direction,
;              measured in m/s. Based on SOT/FPP data.
;        TIMES: the time of each exposure of the EIS scientific data raster.
;
; OPTIONAL INPUTS:
;        None
;
;
; KEYWORD PARAMETERS:
;        None
;
;
; OUTPUTS:
;        Array with line-of-sight velocity of the spacecraft at the times
;        specified by TIMES.
;
;
; OPTIONAL OUTPUTS:
;        None
;
;
; SPECIAL CALLS:
;        FIND_CLOSE()
;
;
; COMMON BLOCKS:
;        None
;
;
; SIDE EFFECTS:
;        None
;
;
; RESTRICTIONS:
;        None
;
;
; EXAMPLE:
;        See EIS_DATA::SETHKPIXCORR
;
;
; MODIFICATION HISTORY:
;      xx-xx-2009: Suguro Kamio
;      09-06-2010: Terje Fredvik - Doc header
;-

function fpp1_dopp_series, times_fpp1, fpp1, times

last = n_elements(times) - 1
pos1 = find_close(times_fpp1, min(times), /earlier)
pos2 = find_close(times_fpp1, max(times), /later)

return,interpol(fpp1[pos1:pos2].dop_vel_rcv, times_fpp1[pos1:pos2], times)

end
