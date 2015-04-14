FUNCTION nan_rebin,array,newsize

  ;; This function rebins the array containing NAN values to a new
  ;; length with correct NAN handling

  ;; 2009/05/31 Reto Stockli (Blue Marble Research)

  cursize=n_elements(array)

  temp=fltarr(cursize/newsize,newsize)
  temp[*,*]=array
  valid=finite(temp)
  newarray=total(temp,1,/nan)/float(total(valid,1)>1)
  a=where(total(valid,1) EQ 0,count)

  IF count GT 0 THEN newarray[a]=!values.f_nan

  RETURN, newarray

END
