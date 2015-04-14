FUNCTION gauss_fun, x,p

  ;; Evaluates Gauss Function value at x
  
  ;; Input: 
  ;; x : abscissa value
  ;; p[0] : amplitude of gauss function
  ;; p[1] : mean value of gauss function
  ;; p[2] : standard deviation of gauss function

  ;; Returns:
  ;; y : gauss function value
  ;; pder : partial derivates of y to parameters p[0] .. p[2]
  
  y = p[0]*exp(-(x-p[1])^2/(2.*p[2]^2))

  pder = [exp(-(x-p[1])^2/(2.*p[2]^2)),0.,0.]

;  RETURN, [y,pder]
  RETURN, [y]

END

