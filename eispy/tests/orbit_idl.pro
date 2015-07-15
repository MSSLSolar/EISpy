pro orbit_idl

  ;This will provide the

t1 = 9.1229760D+08; = anytim('2007-11-29T00:00:00') heater config
t2 = 9.3553920D+08; = anytim('2008-08-24T00:00:00') slit focus adjust
t3 = 9.4057920D+08; = anytim('2008-10-21T08:00:00') grating focus adjust
times = [t1 - 3*30*24*3600., $   ; - 3 months
         t1 + 6*30*24*3600., $
         t2 + 1*30*24*3600., $
         t3 + 6*30*24*3600]

openw, lun_1, './orbit_slit1.txt', /get_lun
openw, lun_2, './orbit_slit2.txt', /get_lun
comment = '# Output from eis_model_series: Time, correction [unknown units]'
printf, lun_1, comment
printf, lun_2, comment

foreach time_i, times do begin
   print, anytim(time_i, /ccsds)
   ; donwload callibration file
   hk_filename = 'eis3_' + strmid(time2file(time_i, /date),0,6) + '.sav'
   sock_copy, 'http://sdc.uio.no/eis_wave_corr_hk_data/' + hk_filename, out_dir='/tmp/'
   restore, '/tmp/'+hk_filename ; time, data
   time_i_array = findgen(30) * 45 + time_i
   result = eis_model_series(time, data, time_i_array, slit=1)
   for j=0, n_elements(result)-1 do printf, lun_1, anytim(time_i_array[j], /ccsds), result[j]
   result = eis_model_series(time, data, time_i_array, slit=2)
   for j=0, n_elements(result)-1 do printf, lun_2, anytim(time_i_array[j], /ccsds), result[j]
endforeach

free_lun, lun_1, lun_2



end
