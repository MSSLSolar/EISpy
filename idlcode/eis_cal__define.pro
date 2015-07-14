;+
; NAME:
;       EIS_CAL__DEFINE
;
; PURPOSE:
;       EIS_CAL__DEFINE defines the class EIS_CAL. Objects of this
;       class contains calibration data and routines 
;       for the EIS instrument on SOLAR-B. The
;       EIS_CAL class inherits the superclass EIS_DATA
;
; CATEGORY:
;       Hansteen/Wikstol Data analysis SW
;
; CALLING SEQUENCE:
;       The EIS_CAL__DEFINE procedure is not called directly. An
;       object of class EIS_CAL is created with the following
;       statement:
;                   eis_cal = obj_new('eis_cal')
;       To fill the object with information (see EIS_DATA__DEFINE for
;       more information about the contents of the object), use the
;       following statement: eis_cal->seteisdefaults 
;       This method is also called on object instatation.
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;       Effective area data can be read in by using the keyword genaeff
;                   eis_cal = obj_new('eis_cal',/gen_aff)
;       The time of the observation and the slit used are necessary to
;       get the proper coeffiecients to transform pixels into wavelengths.
;       They should passed to the eis_cal object through the time_obs and
;       slit_ind keywords:
;        eis_cal = obj_new('eis_cal',time_obs=time_obs,slit_ind=slit_ind)
;       Otherwise, they would be set to the current time and 1" slit by
;       default. In this case a warning message would be issued, unless
;       the quite keyword is set:
;                   eis_cal = obj_new('eis_cal',/quiet)
;
; OUTPUTS:
;       Objects of type EIS_CAL
;
; CALLS:
;
; METHODS:
;       seteisdefaults: sets defaults based on J. Mariska's "EIS Instrument 
;                       Info" document version version 0.85
;       genaeff,detector: reads effective area data for detector 
;       ergs_to_photon(spectrum): converts ergs per lambda to photons per 
;                                 pixel
;       photon_to_dn,lambda(photons): converts photons to data numbers
;
; PROCEDURE:
;
; RESTRICTIONS:
;
;
; MODIFICATION HISTORY:
;   xx-Mar-2004: Oivind Wikstol - Completed version 1.0 of the cal object
;   26-Jan-2006: Viggo Hansteen - Fleshed out the definition. Added several 
;                                 new methods for absolute calibration and 
;                                 wavelength scale computation.
;   22-Jun-2006: Oivind Wikstol - Added dn_to_photon and photon_to_ergs 
;                                 functions.
;   23-Nov-2006: Viggo Hansteen - Changed method seteisdefaults so that it i
;                                 first applies hard coded defaults, then 
;                                 looks in file 'calresp.txt' in 
;                                 self.respdir=getenv('EIS_RESPONSE') 
;                                 for updates 
;  10-Feb-2007: Viggo Hansteen  - Added help methods.
;  17-Feb-2007: Harry P Warren  - Bug fix of quadrant method. 
;  15-Mar-2007: Viggo Hansteen  - cleanup method that actually does some 
;                                 cleanup!
;  28-Mar-2007: Harry P Warren  - Fixed DC correction method
;  29-Sep-2007: A. Gardini      - Removed an unnecessary cleanup method.
;                                 Other changes already made in June 2007.
;  09-Oct-2007: A. Gardini      - Freed pointers in cleanup.
;  18-Oct-2007: A. Gardini      - Set window limits in hot_pixels.
;  23-Oct-2007: A. Gardini      - Added the dusty_pixels method.
;  29-Nov-2007: Viggo Hansteen  - Exchanged shutter open/shutter close 
;                                 exposure time with mhc_dur 
;                                 (now retrieved by getexp()).
;  11-Feb-2008: A. Gardini      - Modified hot_pixels to manage the hot pixel
;                                 files defined on all the pixels in the CCD.
;   5-Mar-2008: A. Gardini      - Substituted the executions of
;                                 calresponse.txt with its inclusion.
;  26-Mar-2008: A. Gardini      - Changed calresponse.txt to calresponse.pro
;                                 because of problems with Windows OS.
;  19-Feb-2008: A. Gardini      - Added the warm_pixels method. Not uploaded.
;   1-Apr-2008: A. Gardini      - Uploaded the warm_pixels method. Updated
;                                 hot_pixels to manage top and bottom parts
;                                 in different directories.
;   4-Apr-2008: A. Gardini      - Set the selection of WP dir like HP dir.
;                                 Corrected the hp, wp generation for an
;                                 offset in the bottom maps (they cover
;                                 [1:512] in Y instead of [0:511]).
;  16-May-2008: P.R. Young      - Previously when removing the dark current
;                                 any negative values were forced to be zero.
;                                 Now negative values are allowed.
;  30-May-2008: A. Gardini      - Corrected the selection of the warm pixels
;                                 file when the exposure is longer than 100s.
;   4-Jul-2008: A. Gardini      - Added HELP and INIT_EIS_CAL_METHODS methods
;  23-Sep-2008: A. Gardini      - Added gainerr, aeffrelerr and methods
;  30-Apr-2009: A. Gardini      - Fixed 1 pixel shift error.
;  30-Apr-2009: A. Gardini      - Set the time_obs, slit_inf, quiet keywords.
;  15-Jul-2009: V. Hansteen     - Added yoff_lw, yoff_sw functions for
;                                 correcting offset between SW and LW spectra
;  18-Feb-2010  V. Hansteen     - Implemented John Mariska's 
;                                 "EIS pointing" ver 0.95 document.
;  29-Apr-2010: T. Fredvik      - Added keyword float to lamb2pix: if set
;                                 return a float pixel value instead
;                                 of an integer.
;  30-Jun-2010: V. Hansteen     - Added new methods getyoff_cmirr,
;                                 getcmirr_move_dates,getlaunch_date,gettau_sensitivity, and
;                                 sensitivity_loss
;  19-Apr-2011: V. Hansteen     - Updated sensitivty tau to 1894 days.
;
;  $Id: eis_cal__define.pro 4175 2014-09-17 15:31:47Z viggoh $
;
;-
function eis_cal::init, caldir=caldir,dc_file=dc_file,ff_file=ff_file,$
         gen_aeff=gen_aeff,time_obs=time_obs,slit_ind=slit_ind,quiet=quiet
  self->init_eis_cal_methods ; initialization of the help method
  self.help=obj_new('hw_help')
; set the default caldir
  if n_elements(caldir) eq 0 then begin
    dirsep = get_delim()
    self.caldir = getenv('EIS_CAL_DATA')+dirsep
  endif
  self.respdir = getenv('EIS_RESPONSE')+dirsep
  self.calrespfile='cal_response.txt'
  if not n_elements(time_obs) then begin
    if not keyword_set(quiet) then begin
      message,"WARNING Observation time not set: using current time.",/info
      message,"        Otherwise, set it by seteisdefaults.",/info
    endif
    get_utc, time_obs
  endif
  if not n_elements(slit_ind) then begin
    if not keyword_set(quiet) then begin
      message,"WARNING Slit indicator not set: using default value.",/info
      message,"        Otherwise, set it by seteisdefaults.",/info
    endif
    slit_ind=self.slit_ind
  endif
  self->seteisdefaults, time_obs, slit_ind
  if n_elements(gen_aeff) ne 0 then begin
    self->genaeff,'A'
    self->genaeff,'B'
  endif
; set dc_file
  if n_elements(dc_file) ne 0 then self.dc_file = dc_file else self.dc_file = ''
; set ff_file
  if n_elements(ff_file) ne 0 then self.ff_file = ff_file else self.ff_file = ''
;
  return,1
end

pro eis_cal::cleanup
  obj_destroy,self.help
  if ptr_valid(self.aux) then begin 
    obj_destroy,*self.aux
    ptr_free,self.aux
  endif
  if ptr_valid(self.hdr) then begin
    obj_destroy,*self.hdr
    ptr_free,self.hdr
  endif
  if ptr_valid(self.cal) then begin
    obj_destroy,*self.cal
    ptr_free,self.cal
  endif
  for i=0,self.nwin-1 do ptr_free,self.w[i]
  ptr_free,self.aeff.lambda_a
  ptr_free,self.aeff.eff_a
  ptr_free,self.aeff.lambda_b
  ptr_free,self.aeff.eff_b
end

pro eis_cal::display_all,init=init
  if n_elements(init) eq 1 then return
  self.help->display_all
  return
end

pro eis_cal::display_methods,init=init
  if n_elements(init) eq 1 then return
  self.help->display_methods
  return
end

;; --- Modified 16-FEB-2007 HPW to work for full ccd case

function eis_cal::quadrant,pos,init=init
  if n_elements(init) eq 1 then return,-1
; compute which ccd quadrant based on pos[0]

  left  = [50L,1074L,2198L,3222L]
  right = left + 1024 - 1
  match = where((pos[0]-left)*(right-pos[0]) ge 0,c)
  if c eq 0 then message,'NO SUCH QUADRANT: CCD POS[0] = '+trim(pos[0]),/info
  
  return,match[0]
end
;------------------------------
;; dark current subtraction
;------------------------------
;
function eis_cal::dc, data, iwin, dcfile=dcfile,verbose=verbose         ,init=init
  if n_elements(init) eq 1 then return,-1
; dark current subraction of window iwin from object data
  if n_elements(dcfile) ne 0 then self.dc_file=dcfile
  if n_elements(verbose) eq 0 then verbose=0
  if self.dc_file eq '' then begin  ; set to dc file name dc_defaultxxx.idlsave
    dirsep = get_delim()
    sdir = file_search(self.caldir,/mark_directory) ; append dirsep to end if not there already
    sdir = sdir+'dc'+dirsep+'dc_default*.idlsave'
    dclist = file_search(sdir,count = nfiles)
    if nfiles gt 0 then begin 
      self.dc_file = dclist[nfiles-1] 
    endif else begin
      message,'No dark current file found, using Q&D default routine',/info
      dc_data={version:'quickndirty'}
    endelse
  endif
  if iwin eq 0 then begin
    message,'using dark current file '+self.dc_file,/info
  endif
;
  data->getwin,iwin,wd,pos
  quad=self->quadrant(pos)
  aux_data=data->getaux_data()
;  exp_dur=((*aux_data.ti_2)-(*aux_data.ti_1))/512.
  exp_dur=data->getexp()
  restore,self.dc_file
  if datatype(dc_data) eq 'UND' then begin
    dc_data={version:'full_ccd'}
    if iwin eq 0 then begin
      message,'CAREFUL, full ccd dark current subtraction has no variation with',/info
      message,'         exposure time nor ccd temperature',/info
    endif
  endif
  
  case dc_data.version of
  
    'quickndirty': begin   

      slit_id = strtrim(data->getslit_id(),2)

      if slit_id eq '1"' or slit_id eq '2"' then begin
        sz = size(wd)
        if sz[1] eq 1024 then begin
          ;; --- full ccd, use low sensitivity areas. this is generally very accurate
          p1 = ([ 39, 944, 39, 926])[iwin] 
          p2 = ([ 84, 989, 84, 971])[iwin] 
          if iwin eq 0 or verbose then $
          message,"SLIT - FULL CCD, USING CONTINUUM REGION TO ESTIMATE THE BACKGROUND",/info
          for i=0,data->getnexp()-1 do begin
            wd_i = reform(wd[p1:p2,*,i])
            dc   = median(wd_i)  
            wd[*,*,i] = (wd[*,*,i] - dc); >0
            if verbose eq 1 then $
            message,"WINDOW = "+trim(iwin,'(i2.2)')+" NEXP = "+trim(i,'(i3.3)')+" DC = "+trim(dc,'(i3.3)'),/info
          endfor
        endif else begin
          ;; --- window, use the median of some low intensity pixels. this is generally an underestimate.
          if iwin eq 0 or verbose then $
          message,"SLIT - WINDOW, USING LOWEST 2% OF INTENSITIES TO ESTIMATE THE BACKGROUND",/info
          npx = 0.02*float(sz[1])*float(sz[2]) 
          for i=0,data->getnexp()-1 do begin
            wd_i = wd[*,*,i]
            idx  = (sort(wd_i))[0:npx-1]
            dc   = median(wd_i[idx])
            wd[*,*,i] = (wd_i - dc); >0
            if verbose eq 1 then $
            message,"WINDOW = "+trim(iwin,'(i2.2)')+" NEXP = "+trim(i,'(i3.3)')+" DC = "+trim(dc,'(i3.3)'),/info
          endfor
        endelse
      endif else begin
        if iwin eq 0 or verbose then $
        message,"SLOT - USING PEDESTAL",/info
        dc = ([549.0,500.0,556.0,547.0])(quad)
        for i=0,data->getnexp()-1 do begin
          wd[*,*,i] = (wd[*,*,i] - dc); >0
          if verbose eq 1 then $
          message,"WINDOW = "+trim(iwin,'(i2.2)')+" NEXP = "+trim(i,'(i3.3)')+" DC = "+trim(dc,'(i3.3)'),/info
        endfor
      endelse

    end

    '1st_order': begin
       if iwin eq 0 then begin
         message,'CAREFUL, 1st order dark current subtraction has no variation with',/info
         message,'         ccd temperature',/info
       endif
       for i=0,data->getnexp()-1 do begin
         wd[*,*,i]=(wd[*,*,i]-(dc_data.coeff[0,quad]+dc_data.coeff[1,quad]*exp_dur[i])); >1
       endfor
                 end
    'full_ccd': begin  
       if pos[0] gt 2147 then begin
         dc=ccda
         pos[0]=pos[0]-2148
       endif else dc=ccdb
       for i=0,data->getnexp()-1 do begin
         wd[*,*,i] = (wd[*,*,i] - dc[pos[0]:pos[0]+pos[1]-1, pos[2]:pos[2]+pos[3]-1]) ; > 0
       endfor
                end
    else: begin
        message,'no such dc.version '+dc_data.version,/info        
          end
  endcase
  return, wd
end


function eis_cal::ff,data,iwin,init=init
  if n_elements(init) eq 1 then return,-1
; flat fielding 
  if self.ff_file eq '' then begin
    dirsep = get_delim()
    sdir = file_search(self.caldir,/mark_directory) ; append dirsep to end if not there already
    sdir = sdir+'ff'+dirsep+'ff_default*.idlsave'
    fflist = file_search(sdir,count = nfiles)
    if nfiles gt 0 then begin
      self.ff_file = fflist[nfiles-1]  
    endif else begin
      ok = dialog_message('No flat field found. Returning.', /information)
      return, -1
    endelse
  endif
  if iwin eq 0 then begin
    message,'using flat field file '+self.ff_file,/info
  endif
;
  data->getwin,iwin,wd,pos
  restore, self.ff_file  ; gives us variables ccda and ccdb     
  if pos[0] gt 2147 then begin
    ff=ccda 
    pos[0]=pos[0]-2148
  endif else ff=ccdb
;  ff=smooth(ff,3)
;  if iwin eq 0 then message,'WARNING!! Flat fielding does nothing (yet).',/info
  for i=0,data->getnexp()-1 do begin
    wd[*,*,i] = wd[*,*,i]/(ff[pos[0]:pos[0]+pos[1]-1, pos[2]:pos[2]+pos[3]-1]>0.9)
  endfor
  return, wd
end


function eis_cal::hot_pixels,data,iwin,init=init
  if n_elements(init) eq 1 then return,-1
; 
; first find the last hp dir that preceeds the observation
;
  dirsep = get_delim()
; if self.hp_dir eq '' then begin
  date=anytim2tai((data->gethdr())->getdate_obs())
  path=getenv('EIS_CAL_DATA')+dirsep+'hp'+dirsep 
  pathlen=strlen(path)
  hp_dirs=file_search(path+'*') 
  if n_elements(hp_dirs) eq 0 then begin
    if iwin eq 0 then message,'no hot pixel directory found ',/info
    return,-1
  endif
  hp_dates=anytim2tai(strmid(hp_dirs,pathlen,10))
;
  dum=min(hp_dates-date,dirindex,/abs)
  self.hp_dir=hp_dirs[dirindex]+dirsep
; endif
  if iwin eq 0 then begin
    message,'Hot pixels dir: '+self.hp_dir,/info
  endif
;
; find right quadrant and choose the file prefix accordingly
  data->getwin,iwin,wd,pos
  quad=self->quadrant(pos)
  case quad of 
  0:pref='coords_b_left'
  1:pref='coords_b_right'
  2:pref='coords_a_left'
  3:pref='coords_a_right'
  else:
  endcase
; check if the hp files are the first ones which cover the CCDs' middle strip
  cm = 0
  hp_dir=hp_dirs[dirindex]+dirsep
  fm=file_search(hp_dir+pref+'*',count=cm)
  fb=file_search(hp_dir+pref+'*bottom*',count=cb)
  ft=file_search(hp_dir+pref+'*top*',count=ct)
  if cm eq 1 and cb eq 0 and ct eq 0 then begin
   left  = [50L,1074L,2198L,3222L]
   restore,fm[0] ; restores ccd_data - UINT  = Array[1024, 512]
   a = pos[0]-left[quad]
   b = pos[0]+pos[1]-left[quad]-1
   c = pos[2]-256L
   d = pos[2]+pos[3]-257L
   ccd_data_info = size(ccd_data)
   mc = 0L
   md = ccd_data_info[2]-1L
   if c lt mc then begin
     message,"WARNING Window extending below of the hot pixels matrix",/info
     message,"        Part exceeding not treated",/info
     c=mc
   endif
   if d gt md then begin
     message,"WARNING Window extending above of the hot pixels matrix",/info
     message,"        Part exceeding not treated",/info
     d=md
   endif
   hp=ccd_data[a:b,c:d]
   return, hp
  endif
; alternatively, now find the CCD's bottom and top files
  cb = 0
  index = dirindex
  repeat begin
    hp_dir=hp_dirs[index]+dirsep
    fb=file_search(hp_dir+pref+'*bottom*',count=cb)
    index = index + 1
  endrep until cb gt 0 or index eq n_elements(hp_dirs)
  if cb eq 0 and dirindex gt 0 then begin
   index = dirindex-1
   repeat begin
     hp_dir=hp_dirs[index]+dirsep
     fb=file_search(hp_dir+pref+'*bottom*',count=cb)
     index = index - 1
   endrep until cb gt 0 or index lt 0
  endif
  if cb eq 1 and iwin eq 0 and hp_dir ne self.hp_dir then $
    message,'Hot pixel CCD bottom files from '+hp_dir,/info
;
  ct = 0
  index = dirindex
  repeat begin
    hp_dir=hp_dirs[index]+dirsep
    ft=file_search(hp_dir+pref+'*top*',count=ct)
    index = index + 1
  endrep until ct gt 0 or index eq n_elements(hp_dirs)
  if ct eq 0 and dirindex gt 0 then begin
   index = dirindex-1
   repeat begin
     hp_dir=hp_dirs[index]+dirsep
     ft=file_search(hp_dir+pref+'*top*',count=ct)
     index = index - 1
   endrep until ct gt 0 or index lt 0
  endif
;
  if ct eq 1 and iwin eq 0 and hp_dir ne self.hp_dir then $ 
     message,'Hot pixel CCD top files from '+hp_dir,/info
  if cb eq 0 and ct eq 0 then begin
   if iwin eq 0 then message,'No hot pixel files found',/info
   return,-1
  endif
; now construct the hot pixel matrix
  left  = [50L,1074L,2198L,3222L]
  hp=uintarr(pos[1],pos[3])
  a = pos[0]-left[quad]
  b = pos[0]+pos[1]-left[quad]-1L
  if cb eq 1 then begin
   restore,fb[0]
   c = pos[2]
   d = pos[2]+pos[3]-1L
   ccd_data_info = size(ccd_data)
   md = ccd_data_info[2]-1L
   if c le md then begin
     if d gt md then d=md
;There is an offset in the bottom maps: they cover [1:512] instead of [0:511]
     hp[0:b-a,0:d-c]=ccd_data[a:b,c-1:d-1]
   end
  endif else if iwin eq 0 then $
      message,"WARNING Bottom half of the hot pixels matrix missing",/info
  if ct eq 1 then begin
   restore,ft[0]
   c = pos[2]-512L
   d = pos[2]+pos[3]-513L
   if d ge 0 then begin
     if c lt 0 then c=0
     hp[0:b-a,pos[3]-1-(d-c):pos[3]-1]=ccd_data[a:b,c:d]
   end
  endif else if iwin eq 0 then $
      message,"WARNING Top half of the hot pixels matrix missing",/info
  return, hp
end

function eis_cal::warm_pixels,data,iwin,init=init
  if n_elements(init) eq 1 then return,-1
; 
; first find the last wp dir that preceeds the observation
;
  dirsep = get_delim()
  date=anytim2tai((data->gethdr())->getdate_obs())
  path=getenv('EIS_CAL_DATA')+dirsep+'wp'+dirsep 
  pathlen=strlen(path)
  wp_dirs=file_search(path+'*') 
  if n_elements(wp_dirs) eq 0 then begin
    if iwin eq 0 then message,'no warm pixel directory found ',/info
    return,-1
  endif
  wp_dates=anytim2tai(strmid(wp_dirs,pathlen,10))
;
; dirindex = max(where (wp_dates le date)) > 0
  dum=min(wp_dates-date,dirindex,/abs)
;
  self.wp_dir=wp_dirs[dirindex]+dirsep
  if iwin eq 0 then begin
    message,'Warm pixels dir: '+self.wp_dir,/info
  endif
;
; find right quadrant and choose the file prefix accordingly
  data->getwin,iwin,wd,pos
  quad=self->quadrant(pos)
  case quad of 
  0:pref='coords_b_left'
  1:pref='coords_b_right'
  2:pref='coords_a_left'
  3:pref='coords_a_right'
  else:
  endcase
; find the max exposure time and choose the file suffix accordingly
  exp_dur=data->getexp()
  exp_dur=max(data->getexp())
  suff=[30,60,100]
  isuff=min(where (suff gt exp_dur-1))
  min_exp_dur=min(data->getexp())
  min_isuff=min(where (suff gt min_exp_dur-1))
  if min_isuff lt isuff and iwin eq 0 then begin
    message,"WARNING  Exposures have different duration.",/info
    message,"Warm pixels taken according to the longest.",/info
  endif
  if isuff eq -1 then begin
    if iwin eq 0 then begin
      message,"Some exposures last longer than warm pixel's ones.",/info
      message,"The longest warm pixels exposure has been selected.",/info
    endif
    isuff = n_elements(suff)-1
  endif
; check if the wp files are the first ones which cover the CCDs' middle strip
  cm = 0
  wp_dir=wp_dirs[dirindex]+dirsep
  fm=file_search(wp_dir+pref+'*middle*',count=cm)
  if cm eq 1 then begin
   left  = [50L,1074L,2198L,3222L]
   restore,fm[0] ; restores ccd_data - UINT  = Array[1024, 512]
   a = pos[0]-left[quad]
   b = pos[0]+pos[1]-left[quad]-1
   c = pos[2]-256L
   d = pos[2]+pos[3]-257L
   ccd_data_info = size(ccd_data)
   mc = 0L
   md = ccd_data_info[2]-1L
   if c lt mc then begin
     message,"WARNING Window extending below of the warm pixels matrix",/info
     message,"        Part exceeding not treated",/info
     c=mc
   endif
   if d gt md then begin
     message,"WARNING Window extending above of the warm pixels matrix",/info
     message,"        Part exceeding not treated",/info
     d=md
   endif
   wp=ccd_data[a:b,c:d]
   return, wp
  endif
; alternatively, now find the CCD's bottom and top files
  cb = 0
  index = dirindex
  repeat begin
    wp_dir=wp_dirs[index]+dirsep
    is = isuff
    repeat begin
    fb=file_search(wp_dir+pref+'*bottom*'+'_'+strtrim(suff[is],2)+'s*',count=cb)
    is = is + 1
    endrep until cb gt 0 or is eq n_elements(suff)
    index = index + 1
  endrep until cb gt 0 or index eq n_elements(wp_dirs)
  if cb eq 0 and dirindex gt 0 then begin
   index = dirindex-1
   repeat begin
     wp_dir=wp_dirs[index]+dirsep
     is = isuff
     repeat begin
     fb=file_search(wp_dir+pref+'*bottom*'+'_'+strtrim(suff[is],2)+'s*',count=cb)
     is = is + 1
     endrep until cb gt 0 or is eq n_elements(suff)
     index = index - 1
   endrep until cb gt 0 or index lt 0
  endif
  if cb eq 1 and iwin eq 0 and wp_dir ne self.wp_dir then $
    message,'Warm pixel CCD bottom files from '+wp_dir,/info
;
  ct = 0
  index = dirindex
  repeat begin
    wp_dir=wp_dirs[index]+dirsep
    is = isuff
    repeat begin
    ft=file_search(wp_dir+pref+'*top*'+'_'+strtrim(suff[is],2)+'s*',count=ct)
    is = is + 1
    endrep until ct gt 0 or is eq n_elements(suff)
    index = index + 1
  endrep until ct gt 0 or index eq n_elements(wp_dirs)
  if ct eq 0 and dirindex gt 0 then begin
   index = dirindex-1
   repeat begin
     wp_dir=wp_dirs[index]+dirsep
     is = isuff
     repeat begin
     ft=file_search(wp_dir+pref+'*top*'+'_'+strtrim(suff[is],2)+'s*',count=ct)
     is = is + 1
     endrep until ct gt 0 or is eq n_elements(suff)
     index = index - 1
   endrep until ct gt 0 or index lt 0
  endif
  if ct eq 1 and iwin eq 0 and wp_dir ne self.wp_dir then $ 
     message,'Warm pixel CCD top files from '+wp_dir,/info
  if cb eq 0 and ct eq 0 then begin
   if iwin eq 0 then message,'No warm pixel files found',/info
   return,-1
  endif
; now construct the warm pixel matrix
  left  = [50L,1074L,2198L,3222L]
  wp=uintarr(pos[1],pos[3])
  a = pos[0]-left[quad]
  b = pos[0]+pos[1]-left[quad]-1L
  if cb eq 1 then begin
   restore,fb[0]
   c = pos[2]
   d = pos[2]+pos[3]-1L
   ccd_data_info = size(ccd_data)
   md = ccd_data_info[2]-1L
   if c le md then begin
     if d gt md then d=md
;There is an offset in the bottom maps: they cover [1:512] instead of [0:511]
     wp[0:b-a,0:d-c]=ccd_data[a:b,c-1:d-1]
   end
  endif else if iwin eq 0 then $
      message,"WARNING Bottom half of the warm pixels matrix missing",/info
  if ct eq 1 then begin
   restore,ft[0]
   c = pos[2]-512L
   d = pos[2]+pos[3]-513L
   if d ge 0 then begin
     if c lt 0 then c=0
     wp[0:b-a,pos[3]-1-(d-c):pos[3]-1]=ccd_data[a:b,c:d]
   end
  endif else if iwin eq 0 then $
      message,"WARNING Top half of the warm pixels matrix missing",/info
  return, wp
end

function eis_cal::dusty_pixels,data,iwin ,init=init
  if n_elements(init) eq 1 then return,-1
; 
; find the dusty pixel file 
  dirsep = get_delim()
  if self.dp_dir eq '' then begin
    ;self.dp_dir='./'
    self.dp_dir=getenv('EIS_CAL_DATA')+dirsep+'dp'+dirsep
  endif
;
  f=file_search(self.dp_dir+'dusty_pixels.sav',count=count)
  if count eq 0 then message,'no dusty pixel file found in '+self.dp_dir
;
  restore,f[0] ; restores dp_data - UINT  = Array[4246,1024] 
  data->getwin,iwin,wd,pos
  dp=dp_data[pos[0]:pos[0]+pos[1]-1,pos[2]:pos[2]+pos[3]-1]
  return, dp
end

pro eis_cal::seteisdefaults,time_obs,slit_ind,init=init
  if n_elements(init) eq 1 then return
; set various instrument parameters, read $EIS_RESPONSE/cal_response.txt      
;
; absolute calibration
;
  self.aeff_version=3
  self.aeff_fileroot=self.respdir + 'EIS_EffArea_'
  self.sr_factor=(725./1.496e8)^2
  self.ergs_to_photons=6.626e-27*2.998e10*1.e8
  self.gain=6.3 ;Hiro 6.45 8/2/2007 changed to 6.3 by Hiro & Louisa 15/4/2007
  self.phot_to_elec=12398.5/3.65
  self.tau_sensitivity = 1894.0 ; updated from 2300.0 on 19/04/2011
;
;pointing info
;
  self.xmidmirrpos=1800    ; home position of the fine mirror, reset to 1800 from 1200 on 08/02/2007
  self.xmidcmirrpos=43703l ; home position of the coarse mirror (NB! unsigned int -> stored in a long)
  self.ymidslitpos=512     ; mid position on the slit (should be 511.5?)
  self.xoffsetut=-129.6    ; offset between AOCS-SOT FOV and EIS midpoint measured in units of 
  self.yoffsetut=-36.3     ; No longer XRT pixels, but arcsec - taken from The List of Solar-B "Hinode" 
                           ; wide FITS Keywords, Dec 25 2006 Masumi Shimojo 
  xrtpixel=1.004           ; 
  self.fmirr_ss = 0.1248   ; fine mirror step size in arcsec (Khalid says 0.12299)
                           ; number should be multiplied by 2 to get image motion 
                           ; (Hiro 0.1203 from observation 8/2/2007)
  self.cmirr_ss = 0.032862 ; coarse mirror step size in arcsec according to Hiro's measurment
                           ; after the first coarse mirror move on 22/01/2007
  self.yoff_sw = [-17.0,8.72e-2] ; yoffsets as linear function of wavelength Kamio-san and Hiro-san
  self.yoff_lw = [-3.49,7.68e-2] ; "Y-offsets of EIS spectra" 18/8/2008
  self.yoff_cmirr = [0.0,-13.7]
  self.cmirr_move_dates = [20060922L,20070122L]

  self.launch_date = '2006-09-22T21:36'

  self.slit_pointing_offset = [0.,0.,8.,0.] ; offset of 2 arcsec slit
;
;define wavelength scale courtesy of Charlie Brown cbrown@ssd5.nrl.navy.mil 16 November 2006.
;
  self.npix = 2148
  nscan = 50
  self.lambda0.a = 199.9389   ; reference wavelength, pixel 0 (including scan pixels detector A)
  self.lambda0.b = 166.131    ; reference wavelength, pixel 0 (including scan pixels detector B)
  self.disp.a = 0.022332      ; dispersion, detector A (long wavelength)
  self.disp.b = 0.022317      ; dispersion, detector B (short wavelength)
  self.dispsq.a = -1.329e-8   ; dispersion, quadratic term, in wavelength equation, detector A (long)
  self.dispsq.b = -1.268e-8   ; dispersion, quadratic term, in wavelength equation, detector B (short)

  ccd_length = findgen(self.npix) - nscan + self.npix - 2*nscan
  self.lambda.scale_a = self.lambda0.a + self.disp.a*ccd_length + self.dispsq.a*ccd_length*ccd_length
  ccd_length = findgen(self.npix) - nscan  
  self.lambda.scale_b = self.lambda0.b + self.disp.b*ccd_length + self.dispsq.b*ccd_length*ccd_length
  
;----------------------------------------------------------------------------------------------------
;  Read calrespfile containing up to date calibration data
;
  problem=0 ; problem? what problem?
;
; define variables in case they are not to be found in an (old) calrespfile
;
  xmidmirrpos=self.xmidmirrpos
  ymidslitpos=self.ymidslitpos
  xoffsetut=self.xoffsetut
  yoffsetut=self.yoffsetut
  fmirr_ss = self.fmirr_ss
  yoff_sw = self.yoff_sw
  yoff_lw = self.yoff_lw
;                                                     
  calrespfile=self.calrespfile
  info=file_info(getenv('EIS_RESPONSE')+path_sep()+calrespfile)
  if not info.exists then begin
    message,'there was some problem reading '+getenv('EIS_RESPONSE')+path_sep()+calrespfile,/info
    message,'using hard coded calibration parameters',/info
    return
  endif
;
; The lines below are substituted by the instruction that follows
; The execute instruction cannot be executed without an IDL license
;
; commands=rd_tfile(getenv('EIS_RESPONSE')+path_sep()+calrespfile,/nocomment)
; for i=0,n_elements(commands)-1 do begin
;   result=execute(commands[i])
;   if not result then problem=1
; endfor
  @eis_cal_response
;
; absolute calibration
;
  self.aeff_version=aeff_version
  self.aeff_fileroot=self.respdir + aeff_fileroot
  self.sr_factor=sr_factor
  self.ergs_to_photons=ergs_to_photons
  self.gain=gain
  self.gainerr=gainerr
  self.aeffrelerr=aeffrelerr
  self.phot_to_elec=phot_to_elec
;
;pointing info
;
  self.xmidmirrpos=xmidmirrpos
  self.ymidslitpos=ymidslitpos
  self.xoffsetut=xoffsetut
  self.yoffsetut=yoffsetut
  self.fmirr_ss = fmirr_ss 
;
;define wavelength scale
;
  self.npix = npix
  nscan=50
  self.lambda0.a = lambda0_a  ; reference wavelength, pixel 0 (including scan pixels detector A)
  self.lambda0.b = lambda0_b  ; reference wavelength, pixel 0 (including scan pixels detector B)
  self.disp.a = disp_a       ; dispersion, detector A (long wavelength)
  self.disp.b = disp_b       ; dispersion, detector B (short wavelength)
  self.dispsq.a = dispsq_a   ; dispersion, quadratic term, in wavelength equation, detector A (long)
  self.dispsq.b = dispsq_b   ; dispersion, quadratic term, in wavelength equation, detector B (short)

  ccd_length = findgen(self.npix) - nscan + self.npix - 2*nscan
  self.lambda.scale_a = self.lambda0.a + self.disp.a*ccd_length + self.dispsq.a*ccd_length*ccd_length
  ccd_length = findgen(self.npix) - nscan  
  self.lambda.scale_b = self.lambda0.b + self.disp.b*ccd_length + self.dispsq.b*ccd_length*ccd_length
;
;define wavelength scale according to data in eis_get_ccd_translation
;
  nscan=0
  coeffs = eis_get_ccd_translation(slit_ind , time = time_obs) 
  self.lambda0.b = coeffs[0].a0
  self.disp.b    = coeffs[0].a1
  self.dispsq.b  = coeffs[0].a2
  self.lambda0.a = coeffs[1].a0
  self.disp.a    = coeffs[1].a1
  self.dispsq.a  = coeffs[1].a2
  ccd_length = findgen(self.npix) - nscan + self.npix - 2*nscan
  self.lambda.scale_a = self.lambda0.a + self.disp.a*ccd_length + self.dispsq.a*ccd_length*ccd_length
  ccd_length = findgen(self.npix) - nscan
  self.lambda.scale_b = self.lambda0.b + self.disp.b*ccd_length + self.dispsq.b*ccd_length*ccd_length
;
  if problem then begin
    message,'there was some problem reading '+getenv('EIS_RESPONSE')+path_sep()+calrespfile,/info
    message,'using hard coded calibration parameters',/info
  endif

  return
end

pro eis_cal::genaeff,detector,init=init
  if n_elements(init) eq 1 then return
; read and store John Ms effective area files

  ;; --- search for files and select the latest version
  detector  = strupcase(detector)
  all_files = file_search(self.aeff_fileroot+detector+'.*',count=nfiles)
  all_files = all_files[sort(all_files)]
  file = all_files[n_elements(all_files)-1]
  message,'CALIBRATION FILE = '+file,/info
  break_file,file,disk,dir,name,ext,/last
  version = fix(str_replace(ext,'.',''))

  ;; --- 
  case detector of 
  'A': begin
     self.aeff.version_a=version
     self.aeff.file_a=file
     data = read_ascii (file,comment_symbol='#')
     self.aeff.lambda_a = ptr_new(fltarr(n_elements(reform(data.field1[0, *]))))
     *self.aeff.lambda_a = reform(data.field1[0,*])
     self.aeff.eff_a = ptr_new(fltarr(n_elements(reform(data.field1[1, *]))))
     *self.aeff.eff_a = reform(data.field1[1,*])
        end
  'B': begin
     self.aeff.version_b=version
     self.aeff.file_b=file
     data = read_ascii (file,comment_symbol='#')
     self.aeff.lambda_b = ptr_new(fltarr(n_elements(reform(data.field1[0, *]))))
     *self.aeff.lambda_b = reform(data.field1[0,*])
     self.aeff.eff_b = ptr_new(fltarr(n_elements(reform(data.field1[1, *]))))
     *self.aeff.eff_b = reform(data.field1[1,*])
       end
  else: begin
     message,'No such detector known, no effective areas loaded',/info
         end
  endcase
  return
end

function eis_cal::ergs_to_photon,lambda,spectrum,init=init
  if n_elements(init) eq 1 then return,-1
; includes geometry and dispersion
  disp =self->getdispersion(lambda)
  return,self.sr_factor/self.ergs_to_photons*lambda*disp*spectrum
end

function eis_cal::photon_to_ergs,lambda,photons,exp=exp,slitsize=slitsize, $
                     days_since_launch=days_since_launch,init=init
  if n_elements(init) eq 1 then return,-1
; includes geometry and dispersion
  if n_elements(exp) eq 0 then exp=fltarr(nrast)+1.0
  if n_elements(slitsize) eq 0 then slitsize=1.0 
  if n_elements(days_since_launch) eq 0 then days_since_launch=0.0
  sz=size(photons)
  intensity=photons
  case sz[0] of 
    3: begin
      nrast=sz[3]
      nslit=sz[2]
       end
    2: begin
      nrast=1
      nslit=sz[2]
       end
    1: begin
      nrast=1
      nslit=1
       end
    else: begin
      message,'max number of dimensions in photons is 3',/info
      return,intensity*0.-1
          end
  endcase
  if slitsize gt 0.0 then begin
    disp = self->getdispersion(lambda)
    lam=lambda*disp#(intarr(nslit)+1.)
  endif else begin
    disp = fltarr(n_elements(lambda))+1.0
    lam = 0.5*(max(lambda)+min(lambda))*disp#(intarr(nslit)+1)
  endelse
  for irast=0,nrast-1 do begin
    intensity[*,*,irast]=self.ergs_to_photons/self.sr_factor/lam*photons[*,*,irast] $
                                /exp[irast]/abs(slitsize)
  endfor
  
  sens_loss=self->sensitivity_loss(abs(days_since_launch))
;
; negative number of days since launch when undoing previous
; correction
;
  if days_since_launch gt 0.0 then intensity=intensity/sens_loss else intensity=intensity*sens_loss

  return,intensity
end

function eis_cal::photon_to_dn,lambda,photons,init=init
  if n_elements(init) eq 1 then return,-1
; includes gain, NB not eff area
  return,self.phot_to_elec/lambda/self.gain*photons
end

function eis_cal::dn_to_photon,lambda,dn,init=init
  if n_elements(init) eq 1 then return,-1
; includes gain, NB not eff area 
  sz=size(dn)
  photon=float(dn)
  case sz[0] of 
    3: begin
      nrast=sz[3]
      nslit=sz[2]
       end
    2: begin
      nrast=1
      nslit=sz[2]
       end
    1: begin
      nrast=1
      nslit=1
       end
    else: begin
      message,'max number of dimensions in dn is 3',/info
      return,photon*0.-1.0
          end
  endcase
  lam=lambda#(intarr(nslit)+1.)
  for irast=0,nrast-1 do begin
    photon[*,*,irast]=self.gain/self.phot_to_elec*lam*float(dn[*,*,irast])
  endfor
  return,photon
end

function eis_cal::lamb2pix, lambda,init=init, float=float
  if n_elements(init) eq 1 then return,-1
; pixel number of given lambda
  if max(lambda) gt 230 then detector='A' else detector='B'
  case detector of
  'A': begin
    pix=(-self.disp.a+sqrt(self.disp.a^2-4.*self.dispsq.a*(self.lambda0.a-lambda)))/2./self.dispsq.a
       end
  'B': begin
    pix=(-self.disp.b+sqrt(self.disp.b^2-4.*self.dispsq.b*(self.lambda0.b-lambda)))/2./self.dispsq.b
       end
  else:
ENDCASE
  
  return, (keyword_set(float)) ? pix : fix(pix)
end

function eis_cal::getdispersion, lambda,init=init
  if n_elements(init) eq 1 then return,-1
; dispersion at given lambda, if n_param() eq 0 returns low precision disp
  if n_elements(lambda) eq 0 then return,0.0223 ; low precision dispersion "for government work"
  pix = self->lamb2pix(lambda)
  if max(lambda) gt 230 then detector='A' else detector='B'
  case detector of 
  'A': begin
    disp = self.disp.a + self.dispsq.a*pix
       end
  'B': begin
    disp = self.disp.b + self.dispsq.b*pix
       end
  else: begin
     message,'This cannot happen',/info
         end
  endcase
  return,disp
end

function eis_cal::getdisp,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.disp
end

function eis_cal::getdispsq,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.dispsq
end

function eis_cal::getnpix,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.npix
end

function eis_cal::getlambda0,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.lambda0
end

function eis_cal::getxmidmirrpos,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.xmidmirrpos
end

function eis_cal::getxmidcmirrpos,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.xmidcmirrpos
end

function eis_cal::getymidslitpos,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.ymidslitpos
end

function eis_cal::getxoffsetut,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.xoffsetut
end

function eis_cal::getyoffsetut,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.yoffsetut
end

function eis_cal::getxoffsetsu,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.xoffsetsu
end

function eis_cal::getyoffsetsu,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.yoffsetsu
end

function eis_cal::getslit_pointing_offset,slit_ind,init=init
  if n_elements(init) eq 1 then return,-1
  if slit_ind eq -1 then begin
    message,'invalid slit_ind, pointing offset set to 0',/info
    return,0.0
  endif
  return,self.slit_pointing_offset[slit_ind]
end

function eis_cal::getfmirr_ss,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.fmirr_ss
end

function eis_cal::getcmirr_ss,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.cmirr_ss
end

function eis_cal::getgain,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.gain
end

function eis_cal::getgainerr,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.gainerr
end

function eis_cal::getaeffrelerr,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.aeffrelerr
end

pro eis_cal::setgain, val, init=init
  if n_elements(init) eq 1 then return
  self.gain = val
end

pro eis_cal::setgainerr, val, init=init
  if n_elements(init) eq 1 then return
  self.gainerr = val
end

pro eis_cal::setaeffrelerr, val, init=init
  if n_elements(init) eq 1 then return
  self.aeffrelerr = val
end

function eis_cal::getyoff_sw,lam,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.yoff_sw[1]*lam+self.yoff_sw[0]
end

function eis_cal::getyoff_lw,lam,init=init
  if n_elements(init) eq 1 then return,-1
  return,self.yoff_lw[1]*lam+self.yoff_lw[0]
end

function eis_cal::getyoff_cmirr,date
  epoch=max(where(long(date) gt self->getcmirr_move_dates()))
  return,(self.yoff_cmirr)[epoch]
end

function eis_cal::getcmirr_move_dates
  return,self.cmirr_move_dates
end

function eis_cal::getlaunch_date
  return,self.launch_date
end

function eis_cal::gettau_sensitivity
  return,self.tau_sensitivity
end

pro eis_cal::settau_sensitivity,tau_sensistivity,init=init
  if n_elements(init) eq 1 then return
  self.tau_sensitivity=tau_sensitivity
end

function eis_cal::getcorrected_sensitivity
  return,self.corrected_sensitivity
end

pro eis_cal::setcorrected_sensitivity,corrected_sensitivity,init=init
  if n_elements(init) eq 1 then return
  self.corrected_sensitivity=corrected_sensitivity
end

function eis_cal::sensitivity_loss,days_since_launch
  if n_elements(days_since_launch) eq 0 then begin 
    days_since_launch=0.
    message,'Days since launch not given, set to zero',/info
  endif
  return,exp(-days_since_launch/self->gettau_sensitivity())
end

pro eis_cal::settau_sensitivity,tau_sensitivity,init=init
  if n_elements(init) eq 1 then return
  self.tau_sensitivity=tau_sensitivity
end

pro eis_cal::setyoff_sw,yoff_sw,init=init
  if n_elements(init) eq 1 then return
  self.yoff_sw=yoff_sw
end

pro eis_cal::setyoff_lw,yoff_lw,init=init
  if n_elements(init) eq 1 then return
  self.yoff_lw=yoff_lw
end

function eis_cal::getaeff,detector,init=init
  if n_elements(init) eq 1 then return,-1
          ; return eff area for detector (in struct)
  case strupcase(detector) of
  'A':begin
     if not ptr_valid(self.aeff.eff_a) then self->genaeff,'A'
     aeff = create_struct('lambda',*self.aeff.lambda_a,'eff',*self.aeff.eff_a)
     end
  'B':begin
     if not ptr_valid(self.aeff.eff_b) then self->genaeff,'B'
     aeff = create_struct('lambda',*self.aeff.lambda_b,'eff',*self.aeff.eff_b)
     end
  else: message,'No such detector: '+strupcase(detector),/info
  endcase
  return,aeff
end

pro eis_cal::help
  help,self,/obj
  return
end

pro eis_cal::init_eis_cal_methods

  self->init_methods
  self->display_all,/init
  self->display_methods,/init
  a=self->quadrant(/init)
  a=self->dc(/init)
  a=self->ff(/init)
  a=self->hot_pixels(/init)
  a=self->warm_pixels(/init)
  a=self->dusty_pixels(/init)
  self->seteisdefaults,/init
  self->genaeff,/init
  a=self->ergs_to_photon(/init)
  a=self->photon_to_ergs(/init)
  a=self->photon_to_dn(/init)
  a=self->dn_to_photon(/init)
  a=self->lamb2pix(/init)
  a=self->getdispersion(/init)
  a=self->getdisp(/init)
  a=self->getdispsq(/init)
  a=self->getnpix(/init)
  a=self->getlambda0(/init)
  a=self->getxmidmirrpos(/init)
  a=self->getxmidcmirrpos(/init)
  a=self->getymidslitpos(/init)
  a=self->getxoffsetut(/init)
  a=self->getyoffsetut(/init)
  a=self->getxoffsetsu(/init)
  a=self->getyoffsetsu(/init)
  a=self->getfmirr_ss(/init)
  a=self->getcmirr_ss(/init)
  a=self->getgain(/init)
  a=self->getgainerr(/init)
  a=self->getaeffrelerr(/init)
  self->setgain,/init
  self->setgainerr,/init
  self->setaeffrelerr,/init
  a=self->getaeff(/init)

  return
end

pro eis_cal__define
  struct = {eis_cal, $
            caldir:' ', $
            respdir:' ', $
            calrespfile:' ',$
            dc_file:'', $
            ff_file:'', $
            hp_dir:'', $
            wp_dir:'', $
            dp_dir:'', $
            aeff_version: 0, $
            aeff_fileroot: ' ',$
            sr_factor: 0.0, $
            ergs_to_photons: 0.0, $
            phot_to_elec: 0.0, $
            gain: 0.0, $
            gainerr: 0.0, $
            aeffrelerr: 0.0, $
            npix: 0, $
            xmidmirrpos:0.0,$
            xmidcmirrpos:0.0,$
            ymidslitpos:0.0,$
            xoffsetut:0.0,$
            yoffsetut:0.0,$
            xoffsetsu:0.0,$
            yoffsetsu:0.0,$
            fmirr_ss:0.0,$ 
            cmirr_ss:0.0,$
            yoff_sw:fltarr(2),$
            yoff_lw:fltarr(2),$
            yoff_cmirr:fltarr(2),$
            cmirr_move_dates:lonarr(2),$
            launch_date:' ', $
            corrected_sensitivity: 0, $
            tau_sensitivity:0.0, $
            slit_pointing_offset:fltarr(4),$
            disp:create_struct(name='disp', $
                                     'a', 0.0d, 'b', 0.0d ), $
            dispsq:create_struct(name='dispsq', $
                                     'a', 0.0d, 'b', 0.0d ), $
            lambda0:create_struct(name='lambda0', $
                                     'a', 0.0d, 'b', 0.0d ), $
            aeff:create_struct(name='aeff', $
                                    'version_a', ' ', $
                                    'lambda_a', ptr_new(), $
                                    'eff_a', ptr_new(), $
                                    'file_a', ' ', $
                                    'version_b', ' ', $
                                    'lambda_b', ptr_new(), $
                                    'eff_b', ptr_new(), $
                                    'file_b', ' '), $
            inherits eis_data}
end
