FUNCTION CORR_TTEST,corr=corr,N=N,sl=sl, ro=ro

; corr:	Pearson's correlation coefficient to be tested
; N: Number of samples
; sl:	Significance level accepted (ex:0.05,0.001,)
; ro:	Correlation value predicted by theory (Null Hypothesis)

; Assumptions
; 1. The N pairs of scores are sampled randomly and independently.
; 2. The distribution of the two variables is bivariate normal.

; NULL hypothesis is ro=ro

corr = double(corr)
sl = double(sl)

IF N_Elements(ro) EQ 0 THEN ro=0.d
IF corr EQ 1.d0 THEN corr=0.9999999d ; avoid Floating divide by 0
IF ro EQ 1.d0 THEN ro=0.9999999d ; avoid Floating divide by 0
IF corr EQ -1.d0 THEN corr=-0.9999999d ; avoid Floating divide by 0
IF ro EQ -1.d0 THEN ro=-0.9999999d ; avoid Floating divide by 0

;----------------------------------------------------------------------------
; COMPUTE CONFIDENCE INTERVAL OF CORRELATION COEFFICIENT
;----------------------------------------------------------------------------
; Conversion of Pearson's correlation to the normally distributed variable zp
; Fisher's transformation

zp=0.5d0*alog((1.d0+corr)/(1.d0-corr))

;Compute zp standard error

sig_zp=1.d0/sqrt(double(N)-3.d0)

; Compute z value from significance level sl
; 99% confidence interval example corresponds sl=0.01 and gives to
z=2.58d0

z=gauss_cvf((sl)/2.d0)

low_zp=zp-z*sig_zp
high_zp=zp+z*sig_zp

r_high=(exp(2.d0*high_zp)-1.d0)/(exp(2.d0*high_zp)+1.d0)
r_low=(exp(2.d0*low_zp)-1.d0)/(exp(2.d0*low_zp)+1.d0)

print,''
print, "High End Case for r: ",r_high," Low End Case for r: ",r_low

; Preliminary result of significance based on Pearson correlation interval
; If the 0 is included in the range between r_low and r_high,
; You can't claim your result is Statistically significant at significance level (sl)
; (or confidence level (1-sl)

IF r_low LT 0. AND r_high GT 0. THEN BEGIN
   print,"This is NOT a statistically significant relationship!" 
   is_significant = 0

ENDIF ELSE BEGIN 
   print,"This is a statistically significant relationship!"
;----------------------------------------------------------------------------

;----------------------------------------------------------------------------
; T test significance
;----------------------------------------------------------------------------
;print,''
;print,'---T TEST result '

; if Null hypothesis is ro=0

IF ro EQ 0. THEN BEGIN

Df=double(N)-2.d0
t=corr*sqrt(Df)/sqrt(1.d0-corr^2)
pt=2.d0*(1.d0-T_PDF(t, Df))

print,df,t,pt

IF pt LT sl THEN $
print,"The correlation with ",string(pt,format='(F12.8)')," is significant repect to significance level ",string(sl,format='(F12.8)') ELSE $
print,"The correlation with ",string(pt,format='(F12.8)')," is NOT significant repect to significance level ",string(sl,format='(F12.8)')

is_significant = pt LT sl

ENDIF ELSE BEGIN

; if Null hypothesis is ro<>0

zpro=0.5d0*alog((1.d0+ro)/(1.d0-ro))
zt=(zp-zpro)/sig_zp
pzt=2.d0*(1.d0-GAUSS_PDF(zt))

IF pzt LT sl THEN $
print,"The null hypothesis that the population correlation is ",string(ro,format='(F12.8)')," can be rejected." ELSE $
print,"The null hypothesis that the population correlation is ",string(ro,format='(F12.8)')," CAN'T be rejected."

is_significant = pt LT sl

ENDELSE

ENDELSE

RETURN, is_significant

END
