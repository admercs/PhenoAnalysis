;+
; :Description:
;   This code performs file locking and unlocking under Unix using IDL  
;
; :Categories:
;   System
;
; :Params:
;   directory: in, required, type="string"
;      directory where to lock or unlock the file
;   filename: in, required, type="string"  
;      file name of file to lock or unlock
;
; :Keywords:
;   lock: in, optional, type=boolean     
;      set this keyword (/lock) to lock the file
;   unlock: in, optional, type=boolean
;      set this keyword (/unlock) to unlock the file
;   verbose: in, optional, type=boolean   
;      set this keyword (/verbose) to print diagnostic messages
;   id: in, optional, type="integer or long"       
;      a file id (e.g. the x/y process on a 2D grid
;   force: in, optional, type=boolean
;      force the use of this file if locked after maxtries have
;      passed. This means that the corresponding lock
;      file is removed with potentially uncool consequences for the
;      process that has locked the file. It may however remedy the
;      situation if the other process has left the file locked and has
;      exited already. 
;
; :Returns: 
; 
; :Uses:
;   uname 
;   touch
;   ln
;   rm 
;   libidl
;
; :Bugs:
;   Warning: This code is not portable to Windows (shared object call on linux and osx)
;
; :Todo:
;
; :Requires:
;   IDL 7.1
;
; :Examples:
;    call the code before opening a file for writing with the /lock
;    keyword. Then call it again after you have finished writing to
;    the file with the /unlock keyword
;
;    file_lock,'mydir','myfile',/lock
;
;    ...
;
;    (open the file myfile, work on the file and close it)
;
;    ...
;
;    file_lock,'mydir','myfile',/unlock  
;
; :History:
;   2011/12/23 : added force keyword to avoid timeouts with leftover
;   file locks
;
;   2011/03/18 : added standard idldoc-compatible rst-style header
;
;   2011/07/25 : added better treatment of temporary (random) files
;
;   including a file ID and a process ID plus a part with the hostname
;
; :Author:
;   Reto Stockli (Blue Marble Research)
;    
; :Copyright:
;   This program is free software: you can redistribute it and/or modify
;   it under the terms of the GNU General Public License as published by
;   the Free Software Foundation, either version 3 of the License, or
;   (at your option) any later version.
;
;   This program is distributed in the hope that it will be useful,
;   but WITHOUT ANY WARRANTY; without even the implied warranty of
;   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;   GNU General Public License for more details.
;
;   You should have received a copy of the GNU General Public License
;   along with this program.  If not, see http://www.gnu.org/licenses/.
;
;-
PRO file_lock,directory,filename,lock=lock,unlock=unlock,verbose=verbose,id=id, force=force

  COMMON LOCK_DATA, lock_file, temp_file, lock_seed, process_id, machine_name

  ;; flags
  maxtries = 10000    ;; number of retries if a file is locked
  maxfiles = 100      ;; maximum number of simultaneously lockable files in a process

;;  verbose = 1

  IF keyword_set(lock) THEN BEGIN

     IF n_elements(lock_seed) EQ 0L THEN BEGIN
        print,"Initializing lock file names"

        ;; 1. initialize random numbers by millisecond system time
        lock_seed = long64(systime(1,/seconds)*1000.d0)

        ;; if several processes initialize at the same millisecond, the
        ;; systime is not enough: also use the process id
        spawn,['uname','-s'],ostype,error,/noshell
        IF error[0] NE '' THEN BEGIN
           print,'File Lock Error : spawning the operating system type was not successful.'
           print,ostype
           print,error
           stop
        ENDIF

        ostype = strtrim(ostype[0],2)
        process_id = 0L
        CASE ostype OF
           'Linux'  : process_id = call_external(!dlm_path+'/libidl.so', 'getpid', /cdecl)
           'Darwin' : process_id = call_external(!dlm_path+'/libidl.dylib', 'getpid', /cdecl)
           ELSE : print,'File Lock Error: Operating System not supported: ',ostype
        ENDCASE
        process_id = string(process_id<999999,format='(I6.6)') 

        ;; also get machine name (in the most extreme case, several
        ;; processes start at the same millisecond with the same
        ;; process id on different hosts
        machine_name = strtrim((get_login_info()).(0),2)

        ;; initialize file arrays
        lock_file = strarr(maxfiles)
        temp_file = strarr(maxfiles)

        print,"Done"
     ENDIF
 
     f=0
     WHILE (f LT (maxfiles-1)) AND (lock_file[f] NE '') DO f += 1

     IF f EQ (maxfiles-1) THEN BEGIN
        print,'Too many locked files for this process.'
        stop
     ENDIF

     lock_file[f] = directory+filename+'.lock'     
     temp_file[f] = directory+'temp.' + process_id + '.' + machine_name + '.' + $
                    string((randomu(lock_seed,/double)+1.d0)*1d15,format='(I20.20)')

     IF keyword_set(id) THEN temp_file[f] += '.'+string(id<999999,format='(I6.6)')
          
     ;; create temp file
     spawn,['touch',temp_file[f]],result,error,/noshell
     IF error[0] NE '' THEN BEGIN
        print,'File Lock Error: touching the temp file ',temp_file[f],' was not successful'
        print,result
        print,error
        stop
     ENDIF
     
     IF keyword_set(verbose) THEN print,'Lock File: ',lock_file[f]
     IF keyword_set(verbose) THEN print,'Temp File: ',temp_file[f]
     
     busy = 1B
     count = 0L
     WHILE (busy EQ 1B) AND (count LE maxtries) DO BEGIN
        
        ;; try to link temp file to lock file
        spawn,['ln','-s',temp_file[f],lock_file[f]],result,error,/noshell
        
        ;; here, a error message means that the file could not be
        ;; linked, thus the original lock file is busy.
        IF keyword_set(verbose) THEN print,'File Lock: ln -s reports: ',result,error

        ;; check if link was successful and both files are the same
        IF file_same(temp_file[f],lock_file[f]) AND (error[0] EQ '') THEN BEGIN
           busy=0B
           IF keyword_set(verbose) THEN print,'File Lock: Locked: ',directory+filename
        ENDIF ELSE BEGIN
           IF keyword_set(verbose) THEN print,'File Lock: Busy:  ',directory+filename
           busy=1B
        ENDELSE
        
        ;; wait for a moment before next try
        IF busy THEN wait,0.5

        IF (count GE maxtries) AND keyword_set(force) THEN BEGIN
           ;; remove lock file, assume locking process has left lock
           ;; file sitting there for no reason

           spawn,['rm','-f',lock_file[f]],result,error,/noshell
           IF error[0] NE '' THEN BEGIN
              print,'File Lock Error: cannot forcefully remove lock file ',lock_file[f]
              print,result
              print,error
              stop
           ENDIF

           ;; reset counter, start new try for locking
           count = 0L
        ENDIF ELSE BEGIN
           ;; update counter
           count +=1L
        ENDELSE

     END

     IF count GT maxtries THEN BEGIN
        print,'File Lock Error: cannot lock file (maxtries exceeded) ',directory+filename
        stop
     ENDIF

  ENDIF

  IF keyword_set(unlock) THEN BEGIN
        
     ;; find lock file name
     f=0
     WHILE (f LT (maxfiles-1)) AND (lock_file[f] NE directory+filename+'.lock') DO f += 1

     IF f EQ (maxfiles-1) THEN BEGIN
        print,'File Lock Error: Lock file not defined: ',directory+filename+'.lock'
        stop
     ENDIF

     ;; check if both the lock and the temp file exist
     IF ~FILE_TEST(temp_file[f]) THEN BEGIN
        print,'File Lock Error: The temp_file ',temp_file[f],' does not exist'
        print,'Tried to unlock: ',directory+filename
        stop
     ENDIF
     IF ~FILE_TEST(lock_file[f]) THEN BEGIN
        print,'File Lock Error: The lock_file ',lock_file[f],' does not exist'
        print,'Tried to unlock: ',directory+filename
        stop
     ENDIF

     ;; remove lock and temporary file
     spawn,['rm','-f',lock_file[f]],result,error,/noshell
     IF error[0] NE '' THEN BEGIN
        print,'File Lock Error: cannot remove lock file ',lock_file[f]
        print,result
        print,error
        stop
     ENDIF

     spawn,['rm','-f',temp_file[f]],result,error,/noshell
     IF error[0] NE '' THEN BEGIN
        print,'File Lock Error: cannot remove temporary file ',temp_file[f]
        print,result
        print,error
        stop
     ENDIF

     temp_file[f] = ''
     lock_file[f] = ''

     IF keyword_set(verbose) THEN print,'File Lock: Unlocked: ',directory+filename

  ENDIF

END
