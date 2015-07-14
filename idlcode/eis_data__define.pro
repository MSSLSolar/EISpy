;+
; NAME:
;       EIS_DATA__DEFINE
;
; PURPOSE:
;       EIS_DATA__DEFINE defines the class EIS_DATA. Objects of this
;       class contains data from the EIS instrument on SOLAR-B. The
;       EIS_DATA class inherits the superclass HW_DATA
;
; CATEGORY:
;       Hansteen/Wikstol Data analysis SW
;
; CALLING SEQUENCE:
;       The EIS_DATA__DEFINE procedure is not called directly. An
;       object of class EIS_DATA is created with the following
;       statement:
;                   eis_data = obj_new('eis_data')
;       To fill the object with information (see HW_DATA__DEFINE for
;       information about the contents of the object), use the
;       following statement:
;                     eis_data-> read, filename [, eis_hdr = eis_hdr]
;
;       One may also create the object and fill it with information
;       with one statement:
;                     eis_data = obj_new('eis_data', filename).
;       In this case the init-procedure of eis_data__define will be
;       run, this in turn calls EIS_DATA__READ
;
;
; INPUTS:
;       file:  Name of file (string) containing EIS data (and
;       headers).
;       Note: Giving filename as input when declaring the EIS_DATA
;       object is optional. Filename can also be sent to EIS_DATA__READ.
;
; KEYWORD PARAMETERS:
;       eis_hdr: Set this keyword to obtain also the header
;       information in the form of an object of type EIS_HDR.
;
; OUTPUTS:
;       Objects of type EIS_DATA (and optionally EIS_HDR)
;
; CALLS:
;       EIS_DATA__READ;
; COMMON BLOCKS:
;
;
; PROCEDURE:
;       The procedure opens an object of class EIS_DATA. The EIS_DATA
;       class inherits the superclass HW_DATA. If filename is given as
;       input when declaring the object, then EIS_DATA__READ is also
;       run by the init-procedure below. In that case the parameters
;       of the object is filled with information.
;
; RESTRICTIONS:
;
;
; MODIFICATION HISTORY:
;       28-Mar-2001: Oivind Wikstol.
;       21-Sep-2001: Oivind Wikstol  - Added documentation.
;       25-Feb-2004: Oivind Wikstol. - Added structure array wd_def[nwin].
;       06-May-2004: Oivind Wikstol  - Added save method
;       11-May-2004: Oivind Wikstol  - Added keyword unit to init-method
;       17-Jun-2004: Oivind Wikstol  - Added filename string and get
;                                      method to eis_data
;       02-Jul-2004:                 - Added dir (directory) to es_data 
;                                      (string) and get method
;       15-Feb-2005:                 - Added number of exp. pr. raster pos.
;       08-Mar-2005: Oivind Wikstol  - Added sit_and_stare variable object
;                                     (0 for raster 1 for sit-and-stare mode)
;       08-Apr-2005: Oivind Wikstol  - Added obsmode:
;                                      e.g. 'sit-and-stare' or 'scanning'
;       03-Jan-2006: Oivind Wikstol  - Added fits_reformat structure and 
;                                      fitslev keyword
;       14-Jul-2006: Oivind Wikstol  - Changd wavelength definitions - moved
;                                      dispersion to cal object.
;	07-Nov-2006: Oivind Wikstol  - Added slit_id
;       23-Nov-2006: Viggo Hansteen  - Added sec_from_obs_start and ti2tai 
;                                      methods
;       16-Dec-2006: Viggo Hansteen  - Added ccsds_packet_time struct and get;                                      method
;       10-Feb-2007: Viggo Hansteen  - Added help methods.
;       09-Mar-2007: Viggo Hansteen  - Added setccd[a,b]_temp methods.
;       15-Mar-2007: Viggo Hansteen  - cleanup method that actually
;                                      does some cleanup!
;       29-Sep-2007: A. Gardini      - Pointers' redefinition and cleanup.
;                                      Other changes made on June 2007.
;       07-Oct-2007: A. Gardini      - Addition of Doppler vel pointers and
;                                      set methods
;       23-Oct-2007: A. Gardini      - Addition of dusty pixels
;       29-Nov-2007: Viggo Hansteen  - Added new version of getexp() method 
;                                      that returns mhc_dur if it exists, 
;                                      this should be a more accurate
;                                      exposure time
;        8-Jan-2008: A. Gardini      - Commented FF references.
;        2-Feb-2008: Mike Marsh      - Added phot tag to calstat structure 
;                                      definition, and setcalstat.
;        8-Feb-2008: A. Gardini      - Added the error keyword in readerr.
;       19-Feb-2008: A. Gardini      - Added the wp tag to calstat structure.
;                                      Uploaded on 1-Apr-2008
;       20-May-2008: A. Gardini      - Added the retain tag to calstat.
;        3-Jul-2008: A. Gardini      - Added HELP and INIT_METHODS methods.
;       13-Nov-2008: A. Gardini      - Added setnexp and setaux_data methods.
;       30-Apr-2009: A. Gardini      - Set keywords in eis_cal.
;       15-Jul-2009: V. Hansteen     - Added getwin function, in
;                                      addition to added functionality
;                                      of getvar, getlam. Added
;                                      correction for yoffsets between
;                                      detectors to getycen.
;       16-Jul-2009: A. Gardini      - Renamed getwin like getwindx.
;        3-Aug-2009: A. Gardini      - Cleaned getvar.
;       20-Jan-2010: V. Hansteen     - Improved mapping, getwindx
;                                      functions
;       18-Feb-2010  V. Hansteen     - Implemented John Mariska's 
;                                      "EIS pointing" ver 0.95 document.
;       14-Apr-2010  T. Fredvik      - Implemented Kamio-san's wavelength
;                                      correction procedures
;       16-Jun-2010  T. Fredvik      - Added check for missing columns in 
;                                      ::sethkpixcorr
;       30-Jun-2010  V. Hansteen     - getycen now calls
;                                      eis_cal::getyoff_cmirr
;       30-Jun-2010  V. Hansteen     - Added new methods
;                                      date_obs2date,
;                                      days_since_launch
;       03-Jul-2010  V. Hansteen     - Added new methods getfovy,
;                                      getfovx. New options 'raster' in get{x,y}cen.
;       10-Aug-2010  T. Fredvik      - Modified ::gethkpixcorr to deal with
;                                      observations with nexp_prp > 1
;       23-Aug-2010  T. Fredvik      - Fixed bug in ::sethkpixcorr when
;                                      calculating slit tilt. Removed
;                                      ::gethkpixcorr, moved code to
;                                      ::sethkpixcorr and ::gethkwavecorr.
;       23-Jun-2011  T. Fredvik      - ::sethkpixcorr: removed uncecessary
;                                      call of eis_get_iwin when extracting
;                                      time array.
;       02-May-2013  T. Fredvik     - ::gethkwavecorr: fixed bug when
;                                      nexp_prp gt 1
;       19-Jun-2013  H. Warren       - ::readerr: modified behavior for changing
;                                      the file name.
;       17-Oct-2013  T. Fredvik      - ::gethkwavecorr: fixed the bug I 
;                                      introduced when trying to fix a bug
;                                      that wasn't a bug... When nexp_prp gt 1
;                                      the wavelength correction arrays now
;                                      have one extra dimension.
;       11-Nov-2103  T. Fredvik      - ::sethkpixcorr: prevent crashing if
;                                      wavelength correction fails (by added 
;                                      keyword INFO when calling message.pro to
;                                      display error message). ::gethkwavecorr: 
;                                      return -1 if sethkpixcorr
;                                      failed  
;       28-Jul-2014  P. Young        - modified sethkpixcorr as the
;                                      wrong input was being given to
;                                      eis_slit_tilt, although the
;                                      result was still the same.
;       22-Oct-2014  T. Fredvik      - in ::readerr, exclude the path when
;                                      creating the name of the error file. 
;
; $Id: eis_data__define.pro 4286 2014-10-22 10:10:33Z tfredvik $
;
;-
;

function eis_data::init, file, datasource=datasource, hdr = hdr,  unit = unit
  self->setcomment,'EIS_data'
  self->init_methods ; initialization of the help method
  self.help=obj_new('hw_help')
  if n_elements(datasource) ne 0 then self.datasource=datasource $
  else self.datasource='fits'
  if n_elements(unit) eq 0 then unit = 'DN'
  self.unit[0] = unit
  self.ccd_sz = [4296, 1024]
  self.home_inst = getenv('EIS_INST')  ; institute where software is run
;
  if n_params() ge 1 then begin
    case self.datasource of
    'ccsds': self-> readccsds, file, hdr = hdr
    'fits' : self-> readfits, file, hdr = hdr
    'xrt'  : self-> readforeign, file, hdr = hdr
    'sot'  : self-> readforeign, file, hdr = hdr
    else: begin
           message,'Unknown datasource '+datasource,/info
           return,-2
         end
    endcase
    self.aux=ptr_new(obj_new('eis_aux'))
; find various calibration data
    time_obs=(self->gethdr())->getdate_obs()
    case self->getdatasource() of
    'fits': slit_ind=self->getinfo('slit_ind')
    'ccsds': slit_ind=self->getslit_ind((*self.hdr[0]->getexp_info()).slit_nr)
    else:
    endcase
    self->setslit_ind,slit_ind
    self.cal=ptr_new(obj_new('eis_cal',time_obs=time_obs,slit_ind=slit_ind))
; check if sensitivity calibration has been done, retrieve e-folding
; time if so...
    if (self->getcalstat()).sens then begin
      (self->getcal())->settau_sensitivity,float(self->getinfo('tau_sens'))
      (self->getcal())->setcorrected_sensitivity,1
    endif
; force reading of error data if fits lever > 0
    if self->getfitslev() gt 0 then self->readerr
;
    *self.cal->setnwin,self.nwin
    *self.cal->setnexp,self.nexp
;define wavelength scale
    lambda=*self.cal->getlambda()
    self.lambda.scale_b=lambda.scale_b
    self.lambda.scale_a=lambda.scale_a
    slit_width=[1.,266.,2.,40.] ; in arcsec
    self->setslit_width,slit_width
    case self->getdatasource() of
    'fits': begin
      slit_ind=self->getinfo('slit_ind')
      fmir_step=self->getinfo('fmir_ss')
      if fmir_step eq 0 then begin
        if slit_ind ne -1 then begin
          self.dx_size=slit_width[slit_ind]
        endif else begin
          message,'invalid slit_ind in hdr, dx_size assumed and set to 1',/info
          self.dx_size=1.0
        endelse
      endif else begin
        self.dx_size=((*self.cal)->getfmirr_ss())*2.0*fmir_step
      endelse
      end
    'ccsds': begin
      fmir_step=(*self.hdr[0]->getexp_info()).fmir_step
      if fmir_step eq 0 then begin
        slit_nr=(*self.hdr[0]->getexp_info()).slit_nr
        if slit_nr ne -1 then self.dx_size=slit_width[self->getslit_ind(slit_nr)] else begin
          message,'invalid slit number in hdr, dx_size assumed and set to 1',/info
          self.dx_size=1.0
        endelse
      endif else begin
        self.dx_size=((*self.cal)->getfmirr_ss())*2.0*fmir_step
      endelse
      end
    else: begin
      self.dx_size=1.0
      end
    endcase
; set default line and continuum wavelength definitions
    for iwin=0,self->getnwin()-1 do begin
      self->setline_px,iwin,[0,(self->getxw())[iwin]-1]
      self->setcont_px,iwin,[0,0]
   ENDFOR
    return, 1
  endif
  message, 'No datafile specified!', /info
  message, 'data->read, file, datasource=datasource',/info
  return,-1
end

pro eis_data::cleanup
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
  ptr_free,self.exp ; Inherited from hw_data
  ptr_free,self.aux_data.ti_1
  ptr_free,self.aux_data.ti_2
  ptr_free,self.aux_data.mhc_dur
  ptr_free,self.aux_data.exp_dur
  ptr_free,self.aux_data.fmirr
  ptr_free,self.aux_data.hslstat
  ptr_free,self.aux_data.xrt_flfl
  ptr_free,self.aux_data.xrt_fl_x
  ptr_free,self.aux_data.xrt_fl_y
  ptr_free,self.aux_data.aec_hepc
  ptr_free,self.aux_data.aec_lepc
  ptr_free,self.aux_data.mhcfmsg
  ptr_free,self.aux_data.xcen
  ptr_free,self.aux_data.ycen
  ptr_free,self.aux_data.ccda_temp
  ptr_free,self.aux_data.ccdb_temp
  ptr_free,self.aux_data.mhc_hz_t10
  ptr_free,self.aux_data.mhc_hz_t15
  ptr_free,self.aux_data.v_sat
  ptr_free,self.aux_data.v_sun
  ptr_free,self.aux_data.v_eth
  ptr_free,self.aux_data.hkpixcorrtilt
  ptr_free,self.aux_data.hkpixcorrtime
  return
end

pro eis_data::display_all,init=init
  if n_elements(init) eq 1 then return
  self.help->display_all
  return
end

pro eis_data::display_methods,init=init
  if n_elements(init) eq 1 then return
  self.help->display_methods
  return
end

pro eis_data::save,file=file,doplan=doplan,init=init
; save current state of data object to fits file
  if n_elements(init) eq 1 then return
  case self.datasource of
    'ccsds': begin
        if n_elements(file) ne 0 then begin
          eis_mkfits,  self,self-> gethdr(),fitsfile=file,doplan=doplan
        endif else eis_mkfits,self,(self-> gethdr()),doplan=doplan
             end
    'fits': begin
        if n_elements(file) ne 0 then begin
          eis_modfits,self,self-> gethdr(),fitsfile=file
        endif else eis_modfits,self,self-> gethdr()
      end
  endcase
  return
end

function eis_data::getaux_data,init=init
; return auxilary data structure, including exposure times etc
  if n_elements(init) eq 1 then return, -1
  return, self.aux_data
end

pro eis_data::setaux_data, aux, init=init
  if n_elements(init) eq 1 then return
  self.aux_data = aux
end

function eis_data::getlambda,lamtype=lamtype,init=init
; lamtype = 'wav' or 'pixel'
  if n_elements(init) eq 1 then return, -1
  if n_elements(lamtype) eq 0 then lamtype='wav'
  if lamtype eq 'wav' then begin
    return,self.lambda
  endif
  if lamtype eq 'pix' then begin
    pixel=self.lambda
    pixel.scale_b=findgen(2148)
    pixel.scale_a=findgen(2148)+2148
    return, pixel
  endif
end

function eis_data::getwindx,input,init=init
  if n_elements(init) eq 1 then return, -1
  if n_params() eq 0 then begin
    message,'no information given!!',/info
    return,-1
  endif
  iwin=intarr(n_elements(input))
  for iw=0,n_elements(input)-1 do begin
    if datatype(input[iw]) eq 'STR' then begin
      iwin[iw]=(where((strupcase(self->getline_id())) eq $
                       strupcase(input[iw]),c))[0]
      if c eq 0 then begin
        message,'Line_id not found : '+input[iw],/info
        iwin[iw]=-1
      endif
    endif else begin
      if input[iw] ge 0 and input[iw] le (self->getnwin())-1 then begin
        iwin[iw]=input[iw]
      endif else begin
;   else e.g. input=195
        nwin=self->getnwin()
        winmax=fltarr(nwin)
        winmin=fltarr(nwin)
        for i=0,nwin-1 do begin 
          winmax[i]=max(self->getlam(i))
          winmin[i]=min(self->getlam(i))
       endfor
        prod=(winmax-input[iw])*(input[iw]-winmin)
        iwin[iw]=(where(prod gt 0,c))[0]
        if c eq 0 then begin
          message,'wavelength not found '+trim(input[iw],'(f10.2)'),/info
          iwin[iw]=-1
        endif 
      endelse
    endelse
  endfor
  return,iwin
end

function eis_data::getvar,iwin,init=init
  if n_elements(init) eq 1 then return,-1
  if n_elements(iwin) eq 0 then iwin=0
  iwin=(self->getwindx(iwin))[0]
  if ptr_valid(self.w[iwin]) then return,*self.w[iwin] else return,-1
end

function eis_data::getlam,iwin,init=init
  if n_elements(init) eq 1 then return, -1
  if n_params() eq 0 then begin
    message,'no window nr input',/info
    iwin=-1
  endif
  iwin=(self->getwindx(iwin))[0]
  if iwin eq -1 then return,-1
  xs=(self->getxs())[iwin]
  xw=(self->getxw())[iwin]
  case xs gt 2148 of
    0: return,self.lambda.scale_b[xs:xs+xw-1]
    1: return,self.lambda.scale_a[xs-2148:xs-2148+xw-1]
  endcase
end

function eis_data::gettau_sensitivity
  return,(self->getcal())->gettau_sensitivity()
end

pro eis_data::settau_sensitivity,tau_sensistivity
  (self->getcal())->settau_sensitivity,tau_sensistivity
end

pro eis_data::correct_sensitivity,tau_sensitivity,reverse=reverse
  if n_elements(tau_sensitivity) ne 0 then self->settau_sensitivity,tau_sensitivity
  if n_elements(reverse) eq 0 then reverse=0
  if (self->getcalstat()).abs then begin
      if reverse then begin 
        if (self->getcalstat()).sens then begin
          message,'Uncorrecting for sensitivity with tau = ' $
             +string((self->getcal())->gettau_sensitivity(),format='(f7.1)')+' days.',/info
          for iwin=0,self->getnwin()-1 do begin
            self->getwin,iwin,intensity,pos
            err=self->geterr(iwin)
;
            days_since_launch=self->days_since_launch()
            sens_loss=(self->getcal())->sensitivity_loss(abs(days_since_launch))
            intensity=intensity*sens_loss
            err=err*sens_loss
            self->setvar,intensity,iwin
            self->seterr,err,iwin
;
            (self->getcal())->setcorrected_sensitivity,0
            self->setcalstat, 'sens', 0
;
            delvarx,err,intensity
          endfor    
        endif else message,'Sensitivity correction has not been done, therefore cannot be reversed',/info
      endif else begin
        if (self->getcalstat()).sens then begin
          message,'Sensitivity correcion already done!',/info
          message,'If redoing, run with correct_sensitivity,/reverse first!',/info
        endif else begin
          message,'Correcting for sensitivity with tau = ' $
             +string((self->getcal())->gettau_sensitivity(),format='(f7.1)')+' days.',/info
          for iwin=0,self->getnwin()-1 do begin
            self->getwin,iwin,intensity,pos
            err=self->geterr(iwin)
;
            days_since_launch=self->days_since_launch()
            sens_loss=(self->getcal())->sensitivity_loss(abs(days_since_launch))
            intensity=intensity/sens_loss
            err=err/sens_loss
            self->setvar,intensity,iwin
            self->seterr,err,iwin
;
            (self->getcal())->setcorrected_sensitivity,1
            self->setcalstat, 'sens', 1
;
            delvarx,err,intensity
         endfor
       endelse
     endelse
     endif else begin
     message,'Not absolute calibrated data, cannot correct for sensitivity loss',/info
     message,'Re-run eis_prep,/correct_sensivity instead',/info
   endelse
end

function eis_data::getcalstat,init=init
  if n_elements(init) eq 1 then return, -1
  return, self.calstat
end

pro eis_data::setcalstat, cal, val ,init=init
; set calibration status ('dc','hp','wp','dp','cr','abs','phot','retain') to 0 or 1
  if n_elements(init) eq 1 then return
  cal=strupcase(cal)
  case cal of
    'DC':self.calstat.dc = val
    'HP':self.calstat.hp = val
    'WP':self.calstat.wp = val
    'DP':self.calstat.dp = val
    'CR':self.calstat.cr = val
;   'FF':self.calstat.ff = val
    'ABS':self.calstat.abs = val
    'SENS':self.calstat.sens = val
    'PHOT':self.calstat.phot = val
    'RETAIN':self.calstat.retain = val
;   'WVL':self.calstat.wvl = val
   endcase
end

function eis_data::getobsmode,init=init
  if n_elements(init) eq 1 then return, -1
  return, self.obsmode
end

function eis_data::getccd_sz,init=init
  if n_elements(init) eq 1 then return, -1
  return, self.ccd_sz
end

function eis_data::getfilename,init=init
  if n_elements(init) eq 1 then return, -1
  return, self.filename
end

function eis_data::getdir,init=init
  if n_elements(init) eq 1 then return, -1
  return, self.dir
end

function eis_data::getline_id,iwin,init=init
  if n_elements(init) eq 1 then return, -1
  if ptr_valid(self.hdr) then begin
    if n_params() eq 0 then return,(*self.hdr)[0]->getline_id() $
    else if self->getwindx(iwin) ne -1 then $
      return,((*self.hdr)[0]->getline_id())[self->getwindx(iwin)]
  endif
  return,' '
end

function eis_data::getslit_id,slit_ind,init=init
  if n_elements(init) eq 1 then return, -1
  if n_elements(slit_ind) eq 0 then slit_ind=self->getslit_ind()
  slit_names=['1"','266"','2"','40"']
  slit_id = 'Unknown'
  if slit_ind ge 0 and slit_ind le 3 then $
    slit_id=slit_names[slit_ind] else begin
    message,'Unknown slit_ind '+strtrim(string(slit_ind),2)+ $
         ' encountered, slit_id set to "unknown"',/info
  endelse
  return, slit_id
end

function eis_data::getnslit,init=init
  if n_elements(init) eq 1 then return, -1
  return, self.nslit
end

function eis_data::getdx_size,init=init
  if n_elements(init) eq 1 then return, -1
  return, self.dx_size
end

function eis_data::getwd_def,init=init
  if n_elements(init) eq 1 then return, -1
  return, self.wd_def
end

function eis_data::getnexp_prp,init=init
  if n_elements(init) eq 1 then return, -1
  return, self.nexp_prp
end

pro eis_data::setnexp, nexp, init=init
  if n_elements(init) eq 1 then return
  self.nexp = nexp
end

function eis_data::getsit_and_stare,init=init
  if n_elements(init) eq 1 then return, -1
  return, self.sit_and_stare
end

pro eis_data::setsit_and_stare, sit_and_stare,init=init
  if n_elements(init) eq 1 then return
  self.sit_and_stare = sit_and_stare
end

pro eis_data::setline_px, iwin, var,init=init
; define line position in window
  if n_elements(init) eq 1 then return
  self.wd_def[iwin].line_px = var
  return
end

pro eis_data::setcont_px, iwin, var,init=init
; define continuum position in window
  if n_elements(init) eq 1 then return
  self.wd_def[iwin].cont_px = var
  return
end

function eis_data::getfitslev,init=init
  if n_elements(init) eq 1 then return, -1
  return,self.fitslev
end

pro eis_data::setfitslev, fitslev,init=init
  if n_elements(init) eq 1 then return
  self.fitslev = fitslev
  return
end


function eis_data::getccsds_packet_time,init=init
  if n_elements(init) eq 1 then return, -1
  return,self.ccsds_packet_time
end

function eis_data::getfits_reformat,init=init
  if n_elements(init) eq 1 then return, -1
  return, self.fits_reformat
end

pro eis_data::setfits_reformat, param, value,init=init
  if n_elements(init) eq 1 then return
  case param of
    'date_rf0':self.fits_reformat.date_rf0 = value
    'date_rf1':self.fits_reformat.date_rf1 = value
    'orig_rf0':self.fits_reformat.orig_rf0 = value
    'orig_rf1':self.fits_reformat.orig_rf1 = value
    'ver_rf0':self.fits_reformat.ver_rf0 = value
    'ver_rf1':self.fits_reformat.ver_rf1 = value
   endcase
  return
end


function eis_data::getinfo,tag,init=init
; get value of 'tag' from fits header
  if n_elements(init) eq 1 then return, -1
  return, (*self.hdr)->getinfo(tag)
end

pro eis_data::setccda_temp,ccda_temp,init=init
  if n_elements(init) eq 1 then return
  *(self.aux_data).ccda_temp=ccda_temp
end

pro eis_data::setccdb_temp,ccdb_temp,init=init
  if n_elements(init) eq 1 then return
  *(self.aux_data).ccdb_temp=ccdb_temp
end

pro eis_data::setmhc_hz_t10,mhc_hz_t10,init=init
  if n_elements(init) eq 1 then return
  if not ptr_valid((self.aux_data).mhc_hz_t10) then $
    self.aux_data.mhc_hz_t10=ptr_new(fltarr(self.nexp))
  *(self.aux_data).mhc_hz_t10=mhc_hz_t10
end

pro eis_data::setmhc_hz_t15,mhc_hz_t15,init=init
  if n_elements(init) eq 1 then return
  if not ptr_valid((self.aux_data).mhc_hz_t15) then $
    self.aux_data.mhc_hz_t15=ptr_new(fltarr(self.nexp))
  *(self.aux_data).mhc_hz_t15=mhc_hz_t15
end

pro eis_data::setv_sat,v_sat,init=init
  if n_elements(init) eq 1 then return
  if not ptr_valid((self.aux_data).v_sat) then $
    self.aux_data.v_sat=ptr_new(fltarr(self.nexp))
  *(self.aux_data).v_sat=v_sat
end

pro eis_data::setv_sun,v_sun,init=init
  if n_elements(init) eq 1 then return
  if not ptr_valid((self.aux_data).v_sun) then $
    self.aux_data.v_sun=ptr_new(fltarr(self.nexp))
  *(self.aux_data).v_sun=v_sun
end

pro eis_data::setv_eth,v_eth,init=init
  if n_elements(init) eq 1 then return
  if not ptr_valid((self.aux_data).v_eth) then $
    self.aux_data.v_eth=ptr_new(fltarr(self.nexp))
  *(self.aux_data).v_eth=v_eth
end

pro eis_data::setxcen,xcen,init=init
  if n_elements(init) eq 1 then return
  *(self.aux_data).xcen=xcen
end

pro eis_data::setycen,ycen,init=init
  if n_elements(init) eq 1 then return
  *(self.aux_data).ycen=ycen
end

pro eis_data::sethkpixcorrtilt,hkpixcorrtilt,init=init
  if n_elements(init) eq 1 then return
  IF NOT ptr_valid((self.aux_data).hkpixcorrtilt) THEN $
     self.aux_data.hkpixcorrtilt = ptr_new(fltarr(n_elements(hkpixcorrtilt)))
  *(self.aux_data).hkpixcorrtilt=hkpixcorrtilt
END

pro eis_data::sethkpixcorrtime,hkpixcorrtime,init=init
  if n_elements(init) eq 1 then return
  IF NOT ptr_valid((self.aux_data).hkpixcorrtime) THEN $
     self.aux_data.hkpixcorrtime = ptr_new(fltarr(n_elements(hkpixcorrtime)))
  *(self.aux_data).hkpixcorrtime=hkpixcorrtime
END

function eis_data::getfovx
  cdelt1=-1.*self->getdx_size()
  fmir_step=self->getinfo('fmir_ss')
  nraster=1
  if fmir_step gt 0 then nraster=self->getnexp()/self->getnexp_prp()
  return,-1.*nraster*cdelt1
end

function eis_data::getfovy
  cdelt2=1.
  return,(self->getyw())[0]*cdelt2
end

function eis_data::getxcen,init=init,raster=raster
  if n_elements(init) eq 1 then return, -1
  if n_elements(raster) eq 0 then raster=0
  case self->getdatasource() of
  'fits': begin
    xcen=*(self->getaux_data()).xcen
    if xcen[0] lt -9000 then xcen=fltarr(self->getnexp())
        end
  'ccsds': begin
    xcen=fltarr(self->getnexp())
         end
  else: begin
    xcen=fltarr(self->getnexp())
        end
  endcase
  xcen=xcen+(self->getinfo('cmirr')- $
    (self->getcal())->getxmidcmirrpos())*((self->getcal())->getcmirr_ss())
  xcen=xcen+(self->getcal())->getslit_pointing_offset(self->getslit_ind())
  if raster then begin
      xcen=xcen[0]+((self->getcal())->getxmidmirrpos() $
      -(*(self->getaux_data()).fmirr)[0])*((self->getcal())->getfmirr_ss())*2. $
      -(self->getfovx())/2.
  endif
  return,xcen
end

function eis_data::getycen,iwin,init=init,raster=raster 
  if n_elements(init) eq 1 then return, -1
  if n_elements(raster) eq 0 then raster=0
  case self->getdatasource() of
  'fits': begin
    ycen=*(self->getaux_data()).ycen
    if ycen[0] lt -9000 then ycen=fltarr(self->getnexp())
    if n_elements(iwin) ne 0 then begin
      iwin=self->getwindx(iwin)
      lam0=(max(self->getlam(iwin))+min(self->getlam(iwin)))/2.0
      if lam0 lt 220. then ycen=ycen+(self->getcal())->getyoff_sw(lam0) $
      else ycen=ycen+(self->getcal())->getyoff_lw(lam0)
    endif
    ycen=ycen+(self->getcal())->getyoff_cmirr(self->getdate())
    if raster then begin
      ycen=ycen[0]-(self->getcal())->getymidslitpos() $
          +(self->getys())[0]+(self->getfovy()+1)/2.-1.0
    endif
          end
  'ccsds': begin
    ycen=fltarr(self->getnexp())
         end
  else: begin
    ycen=fltarr(self->getnexp())
        end
  endcase
  return,ycen
end

function eis_data::getxycen,iwin,init=init
  if n_elements(init) eq 1 then return, -1
  return,{xcen:self->getxcen(),ycen:self->getycen(iwin)}
end

function eis_data::getxpos,init=init
; compute x position as function of exposure nr NB under development
  if n_elements(init) eq 1 then return, -1
  xcen=self->getxcen()
; fmirr_ss multiplied by 2 to get image motion
  xpos=xcen+((self->getcal())->getxmidmirrpos()- $
      (*(self->getaux_data()).fmirr))*((self->getcal())->getfmirr_ss())*2.
;  xpos=xpos-cal->getxoffsetut() ; fixed at FITS generation to EIS
;  pointing
  return,xpos
end

function eis_data::getypos,iwin,init=init
; compute y position as function of exposure nr NB under development
  if n_elements(init) eq 1 then return, -1
  ycen=self->getycen(iwin)
  cal=obj_new('eis_cal',/quiet)
; assume all windows have same height, ie. use [0], base pixel scale on first pointing...
  ypos=ycen[0]-cal->getymidslitpos()+((self->getys())[0]+indgen((self->getyw())[0]))
;  ypos=ypos-cal->getyoffsetut() ;fixed at FITS generation to EIS pointing
  obj_destroy,cal
;  ypos=ycen
  return,ypos
end

function eis_data::sec_from_obs_start,ti,init=init
; input ti assumed ICU in mdp units of 1/512 s, ti_1,ti_2,...
  if n_elements(init) eq 1 then return, -1
  obs_start=self->getinfo('obt_time')
  ti=ulong(ti) ; recasting to deal with old, erroneous data (it can do no harm...)
  dt=(ti-obs_start)
  sub=where(dt lt 0)
  if sub[0] ne -1 then begin ; assume rollover
    message,'Warning ti < date_obs, assuming rollover.',/info
    dt(sub)=ulong(double(ti(sub))+double(ulong(-1))+1-double(obs_start))
  endif
  return,dt/512.
end

function eis_data::getslit_width,slit_ind,init=init
  if n_elements(init) eq 1 then return, -1
if n_elements(slit_ind) eq 0 then slit_ind=self->getslit_ind()
if slit_ind eq -1 then begin
  message,'invalid slit_ind, slit width assumed set to 1',/info
  return,1.0
endif
return,self.slit_width[slit_ind]
end

pro eis_data::setslit_width,slit_width,init=init
  if n_elements(init) eq 1 then return
  if n_elements(slit_width) ne 4 then begin
    message,'slit_width must have 4 elements',/info
    return
  endif
  self.slit_width=slit_width
end

pro eis_data::setslit_ind,slit_ind,init=init
  if n_elements(init) eq 1 then return
  self.slit_ind=slit_ind
end

function eis_data::getslit_ind,slit_nr,init=init
; slit index given slit nr
  if n_elements(init) eq 1 then return, -1
  if n_elements(slit_nr) eq 0 then return,self.slit_ind
;
  date=self->getdate()
;
  slit_ind=-1
  if long(date) le 20080824L then begin
    if slit_nr ge uint('bfff'x) and slit_nr le uint('c080'x) then slit_ind=1
    if slit_nr ge uint('404f'x) and slit_nr le uint('40ef'x) then slit_ind=3
    if slit_nr ge uint('800f'x) and slit_nr le uint('8090'x) then slit_ind=2
    if slit_nr le uint('0080'x) and slit_nr ge 0 or slit_nr ge uint('ffdf'x) then slit_ind=0
  endif else begin
    if slit_nr ge uint('c0c1'x) and slit_nr le uint('c142'x) then slit_ind=1
    if slit_nr ge uint('410c'x) and slit_nr le uint('41ac'x) then slit_ind=3
    if slit_nr ge uint('80cd'x) and slit_nr le uint('81b3'x) then slit_ind=2
    if slit_nr ge uint('00b4'x) and slit_nr le uint('0154'x) then slit_ind=0
    if slit_nr eq uint('ffff'x) then begin
      message,'slit nr eq 0xffff, assuming engineering study, setting slit_ind=0',/info
      slit_ind=0
    endif
  endelse 
  return,slit_ind
end

function eis_data::getti_1,init=init
  if n_elements(init) eq 1 then return, -1
  return,*((self->getaux_data()).ti_1)
end

function eis_data::getti_2,init=init
  if n_elements(init) eq 1 then return, -1
  return,*((self->getaux_data()).ti_2)
end

function eis_data::getdate,init=init
  if n_elements(init) eq 1 then return, -1
  dum=strsplit(self->getfilename(),path_sep(),/extract)
  filename=dum[n_elements(dum)-1]
  date=(strsplit(filename,'_',/extract))[2]
  return,date
end

function eis_data::date_obs2date,date_obs
  if n_elements(date_obs) eq 0 then date_obs=self->getdate_obs()
  date=strmid(strjoin(strsplit(date_obs,'-T:',/extract),''),0,8)
  return,date
end

function eis_data::getdate_obs,init=init
  if n_elements(init) eq 1 then return, -1
  return,(self->gethdr())->getdate_obs()
end

function eis_data::days_since_launch,date
  if n_elements(date) eq 0 then date=self->getdate_obs()
  jd_date=(anytim2jd(date)).int
  jd_launch=(anytim2jd((self->getcal())->getlaunch_date())).int
  return,jd_date-jd_launch
end                             

function eis_data::ti2tai,ti,init=init
; time ti, ICU in mdp units of 1/512 s, converted to atomic time units (tai) (sort of....)
  if n_elements(init) eq 1 then return, -1
  if n_elements(ti) eq 0 then ti=*(self->getaux_data()).ti_1
  return,anytim2tai(self->getinfo('date_obs'))+self->sec_from_obs_start(ti)
end

function eis_data::ti2utc,ti,init=init
  if n_elements(init) eq 1 then return, -1
  if n_elements(ti) eq 0 then ti=*(self->getaux_data()).ti_1
  return,anytim2utc(self->ti2tai(ti),/time_only,/ccsds,/truncate)
end

function eis_data::getcal,init=init
  if n_elements(init) eq 1 then return, -1
  return,*self.cal
end

function eis_data::geterr,win,init=init
  if n_elements(init) eq 1 then return, -1
  win=self->getwindx(win)
  return,(self->getcal())->getvar(win)
end

pro eis_data::seterr,err,win,init=init
  if n_elements(init) eq 1 then return
  (self->getcal())->setvar,err,win
end

pro eis_data::saveerr,file=file,init=init
  if n_elements(init) eq 1 then return
  eis_modfits,*self.cal,self-> gethdr(),fitsfile=file,/noaux
end

pro eis_data::readerr,file=file,error=error,init=init
  if n_elements(init) eq 1 then return
  if n_elements(file) eq 0 then begin
    file=self->getfilename()
      ;; Position of the start of the filename, excluding the path.
    fnameix = last_nelem(strsplit(file,'/'))
    if strpos(file,'_er_') eq -1 then $
       strput,file,'er',stregex(strmid(file,fnameix),'_l._')+1+fnameix
    file=stregex(file,'^.*.(fits|fits.gz)',/extract)
    message,'Restoring errfile "'+file+'"',/info
  endif
  (self->getcal())->readfits,file,/noaux,error=error
end

function eis_data::getexp,init=init
  if n_elements(init) eq 1 then return, -1
  if ptr_valid(self.aux_data.mhc_dur) then begin
    if (*self.aux_data.mhc_dur)[0]/1.e6 lt 4000.0 then return,*self.aux_data.mhc_dur/1.e6 $
    else return,*self.exp
  endif else begin
    message,'no MHC duration exposure time found, returning requested exposure time',/info
    return,*self.exp
  endelse
end

;; **********************************************************************
;; START methods used when correcting the data for orbital variation of the line
;; centre. We use Kamio-san's code to correct the data using house keeping
;; temperatures instead of the acutal observed data. External calls:
;; eis_model_serie and fpp1_dopp_series. 
;; ***********************************************************************


PRO eis_data::sethkpixcorr,  init=init
  ;;  
  ;; Apply Kamio-san's procedures to calculate the wavelength correction due
  ;; to the temperature variations within the instrument (measured in pixels)
  ;; as a function of raster exposure number. Store this array in
  ;; self.pixcorrtime. Also calculate the slit tilt and store in
  ;; self.pixcorrtilt. To retrieve these arrays and a 2D wavelength correction
  ;; array combining the two, in Angstrom, see ::gethkwavecorr.
  ;;
  IF n_elements(init) EQ 1 THEN return
  
  lamFeXII = 195.12
  ;; REMOVE THE NEXT TWO LINES WHEN FINDING A METHOD THAT RETURNS TIMES!
  wd = self->getwindata(0)
  time = anytim(wd.time_ccsds)
  
  
  yws = self->getinfo('yws')   ; y window start
  ny = median((self->getyw())) ; number of elements along y
  date = self->getinfo('date_obs')
  slit_ind = self->getinfo('slit_ind')
  IF slit_ind EQ 0 THEN slit=1 ELSE slit=2   ; PRY, 28-Jul-2014
  
  lamshort = self->getinfo('TWAVE1') ;; Principal wavelength of first data window
  lamlong = self->getinfo('TWAVE'+TRIM(self->getinfo('nwin'))) ;; lam of last win
 ;
 ; PRY, 28-Jul-2014, I've modified the lines below to give
 ; slit=slit instead of
 ; slit=self->getinfo('slit_ind'), which was the wrong input (although
 ; the result was the same)
 ;
  lamcorrtiltshort = eis_slit_tilt(yws, ny, date=date, $
                                   slit=slit, /short)
  lamcorrtiltlong = eis_slit_tilt(yws, ny, date=date, $
                                  slit=slit, /long)
  
  ;; A rather nasty hack... eis_slit_tilt returns slit tilt in Angstroms, not in
  ;; pixels, and there is one tilt array for detector A, one for B. For the time
  ;; being we transform the tilt arrays from wavelenght to pixels, and then
  ;; back to wavelenght again in ::gethkwavecorr... If eis_slit_tilt in a future
  ;; update returns tilt as a funciton of wavelength, not only detector, this
  ;; way of doing it may be useful. Now it might be a bit confusing...
  pixcorrtilt = [[self->dlam2dpix(lamcorrtiltshort,lamshort)],$
                 [self->dlam2dpix(lamcorrtiltlong,lamlong)]]
  

  
  self->getstatus_data, self->getinfo('date_obs'), time_eis3, eis3, erreis, /eis
  self->getstatus_data, self->getinfo('date_obs'), time_fpp1, fpp1, errsot, /sot
  err = (erreis NE '' AND errsot NE '') ? erreis + ', ' + errsot : erreis + errsot
  
  IF err EQ '' THEN BEGIN 
     goodexp = where(self->check_ti() NE 0, complement=badexp, ngoodexp, $
                   ncomplement=nbadexp)
     IF self->getsit_and_stare() EQ 0 THEN BEGIN
        IF ngoodexp NE 0 THEN goodexp = reverse(self->getnexp()-1-goodexp)
        IF nbadexp NE 0 THEN badexp = reverse(self->getnexp()-1-badexp)
     ENDIF
     
     pixel = eis_model_series(time_eis3, eis3, time, $
                              slit = self->getinfo('slit_ind'), $
                              goodexp=goodexp, badexp=badexp)
     
     dopp = fpp1_dopp_series(time_fpp1 , fpp1, time)  
     cal = self->getcal()
     
     dispersion = cal->getdispersion(lamFeXII)
     shift_dopp = dopp / 3e8 * lamFeXII / dispersion
     pixel += shift_dopp
     
     pixcorrtime = pixel - (self->getcal())->lamb2pix(lamFeXII,/float)
     
     IF nbadexp GT 0 THEN pixcorrtime[badexp] = 0
     
     self->sethkpixcorrtilt, pixcorrtilt
     self->sethkpixcorrtime, pixcorrtime
  ENDIF ELSE message,'Could not perform wavelength correction on EIS data.' + err,/info
END


FUNCTION eis_data::dlam2dpix, lamcorr, lam
  ;; Help method called by sethkpixcorr. Converts an array with wavelength
  ;; differences for a line with wavelength lam, to an CCD pixel position
  ;; difference array. This operation will most likely only be done to convert
  ;; the slit tilt array found by eis_slit_tilt in sethkpixcorr
  
  lamcorr = double(lamcorr) + lam
   
  sz = size(lamcorr)  
   
  pix = dblarr(sz[1])
  cal = self->getcal()
  FOR i=0,sz[1]-1 DO pix[i] = cal->lamb2pix(lamcorr[i],/float)
   
  pix -= cal->lamb2pix(lam,/float)
  
  lamcorr -= lam
  
 return, pix
END


PRO eis_data::getstatus_data, date, time, data, err, eis=eis, sot=sot
  ;; Check if a house keeping status data file (EIS or SOT) for this
  ;; month is stored locally. First search in EIS_WAVE_CORR_HK_DATA (if this
  ;; environment variable is set), then search in the eis_wave_corr_hk_data
  ;; directory in $EIS_DATA then finally search in the current working
  ;; directory. If the file is found, and it has the same size as the version
  ;; of the file that is stored at the Hinode Science Data Centre Europe Oslo,
  ;; restore the file. If the file is not found or it has a different size
  ;; than the file at SDC, download the file from Oslo, then restore the file.
  date = anytim(date,/ccsds)
  yyyy = strmid(date,0,4)
  mm = strmid(date,5,2)
  
  statusdir = [getenv('EIS_WAVE_CORR_HK_DATA'), $
               concat_dir(getenv('EIS_DATA'),'wave_corr')]
  statusdir = [statusdir[where(statusdir NE '')],'']
  
  instrtxt = (keyword_set(eis)) ? 'eis3_' : 'fpp1_'
  
  err = ''
  i = 0
  keepsearching = 1
  
  WHILE keepsearching AND i LT n_elements(statusdir) DO BEGIN 
     file = concat_dir(statusdir[i], instrtxt + yyyy + mm+'.sav') 
     IF file_test(file) THEN BEGIN 
        IF self->same_size(file) THEN BEGIN
           ;; i.e. 1) the local and SDC file has the same size, or 2) could
           ;; not open socket to server and we therefore assume that the local
           ;; file is ok.
           restore,file = file
           keepsearching = 0
        ENDIF ELSE keepsearching = 1
     ENDIF
     i++
  ENDWHILE 
  
  IF keepsearching THEN BEGIN
     self->download_status_data, downloaddir, file, timearr, dataarr, err 
     IF err EQ '' THEN restore,concat_dir(downloaddir, file)
  ENDIF
  
  timediff = anytim(date) - last_nelem(time)
  IF timediff GT 0 THEN err += ' No housekeeping data exists for given date (' + $
                             strmid(instrtxt,0,3) + ' housekeeping data ends ~'+ $
                             trim(round(timediff/3600.)) + $
                             ' hours prior to the date_obs of the present file)'
END


FUNCTION eis_data::same_size, file
  ;; Help method called by getstatus_data. Check if the local version of file
  ;; is the same as the version in Oslo
  o = obj_new('http')
  o->open,'sdc.uio.no', err=err
  url = 'http://sdc.uio.no/eis_wave_corr_hk_data/'
  IF err EQ '' THEN same_size = $
     o->same_size(url+file_break(file), file) ELSE BEGIN 
     message,err,/info
     message,'Assuming that the local house keeping status file is ok.',/info
     same_size = -1 
  ENDELSE
  
  obj_destroy,o
  
  return, same_size
END


PRO eis_data::download_status_data, downloaddir, file, time, status, err
  ;; Help method called by getstatus_data if status data is not found on
  ;; harddisk. Download house keeping status data from Hinode Science Data
  ;; Centre Europe, Oslo.
  
  ;; First give user some friendly advice
  downloaddir = getenv('EIS_WAVE_CORR_HK_DATA')
  complaintxt = (downloaddir EQ '') ? ', or setting the environment '+$
                'variable EIS_WAVE_CORR_HK_DATA to the directory where you have '+$
                'downloaded the house keeping data.' : '.'
  instrumenttxt = strupcase(strmid(file,0,3))
   
  message,'The required ' + instrumenttxt + ' house keeping status data was not'+$
          ' found on disk, trying to download from Oslo Hinode archive.' + $
          ' Please consider running ssw_update' + complaintxt, /info
  
  ;; Now that we're done with the pleasantries, start the downloading!
  o = obj_new('http')
  o->open,'sdc.uio.no', err=err
  IF err EQ '' THEN $
     o->copy,'http://sdc.uio.no/eis_wave_corr_hk_data/' + file, $
             out_dir=downloaddir,err=err,/verbose ELSE message,err,/info
  obj_destroy,o
END


FUNCTION eis_data::gethkwavecorr, lam, wvl_cube=wvl_cube
  ;; Return a structure consisting of: IM: a 2D wavelength correction array
  ;; (measured in Aangstroms) constructed from the arrays TIME: the wavelength 
  ;; correction due to temperature variations, and TILT: the slit tilt. If
  ;; keyword wvl_cube is present, lamcorr will also include a 3D wavlenght cube
  ;; WVL_CUBE where wvl_cube[*,x,y] is the corrected 1D wavlength array for pixel
  ;; (x,y) in the raster image. wvl_cube is in the format used by the cfit line 
  ;; fitting package.
  
  IF n_params() EQ 0 THEN lam = 195.12
  
  
  IF ~ptr_valid((self->getaux_data()).hkpixcorrtime) THEN self->sethkpixcorr
  
  aux_data = self->getaux_data()
  
  ;; If and only if sethkpixcorr was successful then aux_data.hkpixcorrtime
  ;; is a valid pointer. Return lamcorr=-1 if sethkpixcorr failed.
  IF ~ptr_valid(aux_data.hkpixcorrtime) THEN lamcorr = -1 ELSE BEGIN 
     pixcorrtime = *(aux_data.hkpixcorrtime)
     lamcorrtime = self->dpix2dlam(pixcorrtime, lam)
     
     pixcorrtilt = *(aux_data.hkpixcorrtilt)
     long = (lam GT 230) ? 1 : 0
     conversionlam = (long) ? self->getinfo('TWAVE'+TRIM(self->getinfo('nwin'))) : $
                     self->getinfo('TWAVE1')
     pixcorrtilt = pixcorrtilt[*,long]
     
     lamcorrtilt = self->dpix2dlam(pixcorrtilt,conversionlam)
     
     sztilt = size(lamcorrtilt)
     sztime = size(lamcorrtime)
     
     IF sztime[0] EQ 1 THEN BEGIN
        lamim = dblarr(sztime[1],sztilt[1])
        FOR y=0,sztilt[1]-1 DO lamim[*,y] = lamcorrtime + lamcorrtilt[y] 
     ENDIF ELSE BEGIN ;; when nexp_prp greater than 1
        lamim = dblarr((size(lamcorrtime))[1],sztilt[1],(size(lamcorrtime))[2])
        FOR y=0,sztilt[1]-1 DO lamim[*,y,*] = lamcorrtime + lamcorrtilt[y]
     ENDELSE
     
     IF keyword_set(wvl_cube) THEN $
        lamcorr = {im:lamim, time:lamcorrtime, tilt:lamcorrtilt, lam:lam, $
                   cube:self->getwvl_cube(lam, lamim)} $ 
     ELSE $
        lamcorr = {im:lamim, time:lamcorrtime, tilt:lamcorrtilt, lam:lam}
  ENDELSE 
  
  return, lamcorr
  
END



FUNCTION eis_data::dpix2dlam, pix, lam_i
  ;; Help method called by gethkwavecorr. Converts an array with CCD pixel
  ;; position differences for a line with wavelength lam, to wavlength
  ;; difference array (measured in Aangstroms).
  
   cal = self->getcal()
  
   px = cal->lamb2pix(lam_i,/float) + pix
   
   IF lam_i GT 230 THEN $ 
      lam = (cal->getlambda0()).a + (cal->getdisp()).a*px + $
            (cal->getdispsq()).a*px*px $
   ELSE $ 
      lam = (cal->getlambda0()).b + (cal->getdisp()).b*px + $
            (cal->getdispsq()).b*px*px
   
   return, lam - lam_i
END


FUNCTION eis_data::getwvl_cube, lam, lamim
  ;; Help method called by gethkwavecorr. Returns a 3D/4D wavelenght cube where
  ;; wvl_cube[*,x,y,z] is the corrected 1D wavlength array for pixel (x,y) and
  ;; exposure number z in the raster image. wvl_cube is in the format used by 
  ;; the cfit line fitting package.
  iwin = eis_get_iwin(self->getfilename(),lam)
  IF iwin NE -1 THEN BEGIN 
     wd = self->getwindata(iwin)     
     lamarr = wd.wvl
     ;; Determine the size of the wvl_cube, 3D if lamim is 2D, 4D if lamim is
     ;; 3D (the latter is the case for nexp_prp gt 1)
     ndim = (size(lamim))[0]
     sz = [(size(lamarr))[1], (size(lamim))[1:ndim]]

     wvl_cube = fltarr(sz)
     ;; The three stars notation works also if wvl_cube i 3D
     FOR j = 0,sz[0]-1 DO wvl_cube[j,*,*,*] = lamarr[j] - lamim
     return, wvl_cube
  ENDIF ELSE return, -1
END
    

;; **********************************************************************
;; END methods used when correcting the data for orbital variation of the line
;; centre.
;; ***********************************************************************


function eis_data::mk_eis_map,win,mom=mom
  id=strtrim(self->getslit_id(),2)
  if id eq '40"' or id eq '266"' then begin
    win=self->getwindx(win)
    xw=(self->getxw())[win]
    yw=(self->getyw())[0]
    im=fltarr(self->getnexp()*xw,yw)
    dx=1.
    for i=0,self->getnexp()-1 do begin
      im[xw*i:xw*(i+1)-1,*]=(self->getvar(win))[*,*,i]
    endfor
    im=rotate(im,5)
  endif else begin
    im=rotate(total(self->getvar(win),1),1)
    dx=(self->getxpos())[0]-(self->getxpos())[1]
  endelse
  xc=self->getxcen(/raster)
  yc=self->getycen(win,/raster)
  dy=((self->getypos())[1]-(self->getypos())[0])
  time=self->getinfo('DATE_OBS')
  eis_map=make_map(im,xc=xc,yc=yc,dx=dx,dy=dy,time=time)
  add_prop,eis_map,id='EIS '+self->getline_id(win),/replace
  dur=max(self->ti2tai(self->getti_2()))-min(self->ti2tai(self->getti_1()))
  add_prop,eis_map,dur=dur,/replace
  ang=pb0r(time,/arcsec)
  add_prop,eis_map,b0=ang[1]
  add_prop,eis_map,rsun=ang[2]
  return,eis_map
end

pro eis_data::help,var,_extra=ex
  help,self,/obj,_extra=ex 
  return
end

pro eis_data::init_methods

  self->init_hw_methods
  self->display_all,/init
  self->display_methods,/init
  self->save,/init
  a=self->getaux_data(/init)
  self->setaux_data,/init
  a=self->getlambda(/init)
  a=self->getwindx(/init)
  a=self->getlam(/init)
  a=self->getcalstat(/init)
  self->setcalstat,/init
  a=self->getobsmode(/init)
  a=self->getccd_sz(/init)
  a=self->getfilename(/init)
  a=self->getdir(/init)
  a=self->getline_id(/init)
  a=self->getslit_id(/init)
  a=self->getti_1(/init)
  a=self->getti_2(/init)
  a=self->getdate(/init)
  a=self->getnslit(/init)
  a=self->getdx_size(/init)
  a=self->getwd_def(/init)
  a=self->getnexp_prp(/init)
  self->setnexp,/init
  a=self->getsit_and_stare(/init)
  self->setsit_and_stare,/init
  self->setline_px,/init
  self->setcont_px,/init
  a=self->getfitslev(/init)
  self->setfitslev,/init
  a=self->getccsds_packet_time(/init)
  a=self->getfits_reformat(/init)
  self->setfits_reformat,/init
  a=self->getinfo(/init)
  self->setccda_temp,/init
  self->setccdb_temp,/init
  self->setmhc_hz_t10,/init
  self->setmhc_hz_t15,/init
  self->setv_sat,/init
  self->setv_sun,/init
  self->setv_eth,/init
  self->setxcen,/init
  self->setycen,/init
  a=self->getxcen(/init)
  a=self->getycen(/init)
  a=self->getxycen(/init)
  a=self->getxpos(/init)
  a=self->getypos(/init)
  a=self->sec_from_obs_start(/init)
  a=self->getslit_ind(/init)
  a=self->ti2tai(/init)
  a=self->ti2utc(/init)
  a=self->getcal(/init)
  a=self->geterr(/init)
  self->seterr,/init
  self->saveerr,/init
  self->readerr,/init
  a=self->getexp(/init)
  a=self->check_ti(/init)
  a=self->dn_to_ergs(/init)
  self->read,/init
  self->readccsds,/init
  self->readfits,/init
  self->readforeign,/init
  return
end

pro eis_data__define
  nwin = 25
  wdstruct = create_struct(name = 'wd_def',  'line_px', intarr(2), $
                           'cont_px', intarr(2))
  struct = {eis_data, $
;           cal:ptr_new(obj_new()), $
            cal:ptr_new(), $ ; This definition prevents the dangling pointer
            obsmode:'', $
            slit_id:'', $
            slit_width:fltarr(4), $
            slit_ind:0, $
            nslit:0, $
            nexp_prp:0, $
            dx_size:0.0, $
            dy_size:0.0, $
            sit_and_stare:0, $
            ccd_sz:intarr(2), $
            filename:'', $
            dir:'', $
            lambda:create_struct(name='lambda', $
                                    'scale_b', dblarr(2148), $
                                    'scale_a', dblarr(2148)), $
            aux_data:create_struct(name='auxdata', $
                                    'ti_1', ptr_new(), $
                                    'ti_2', ptr_new(), $
                                    'mhc_dur', ptr_new(),$
                                    'exp_dur', ptr_new(), $
                                    'fmirr', ptr_new(), $
                                    'hslstat', ptr_new(), $
                                    'xrt_flfl', ptr_new(), $
                                    'xrt_fl_x', ptr_new(), $
                                    'xrt_fl_y', ptr_new(), $
                                    'aec_hepc', ptr_new(), $
                                    'aec_lepc', ptr_new(), $
                                    'mhcfmsg', ptr_new(), $
                                    'xcen', ptr_new(), $
                                    'ycen', ptr_new(), $
                                    'ccda_temp', ptr_new(), $
                                    'ccdb_temp',ptr_new(), $
                                    'mhc_hz_t10',ptr_new(), $
                                    'mhc_hz_t15',ptr_new(), $
                                    'v_sat',ptr_new(), $
                                    'v_sun',ptr_new(), $
                                    'v_eth',ptr_new(), $
                                    'hkpixcorrtilt', ptr_new(), $
                                    'hkpixcorrtime', ptr_new()), $
            calstat:create_struct(name='calstat', $
                                      'dc', 0, $
;                                     'ff', 0, $
                                      'hp', 0, $
                                      'wp', 0, $
                                      'dp', 0, $
                                      'cr', 0, $
;                                     'wvl', 0, $
                                      'abs', 0, $
                                      'sens',0, $
                                      'phot',0, $
                                      'retain',0), $
            fits_reformat:create_struct(name='fits_reformat', $
                                             'date_rf0','NA', $
                                             'orig_rf0','NA', $
                                             'ver_rf0','NA', $
                                             'date_rf1','NA', $
                                             'orig_rf1','NA', $
                                             'ver_rf1','NA'), $
            ccsds_packet_time:create_struct(name='ccsds_packet_time', $
                                            'first',0ul, $
                                            'last',0ul), $
            fitslev:0, $
            wd_def:replicate(wdstruct, nwin), $
            inherits hw_data}
end

