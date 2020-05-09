program sqhole
implicit none
!>
!>the spectra of a layer of the one-dimensional periodic array of circular PEC rod
!>
real(8), parameter    :: pi = 3.d0*acos(0.5d0)
complex(8), parameter :: iu=(0.d0,1.d0) 
real(8), parameter :: c = 299792458.d0         !speed of light (m/s)
integer, parameter :: nmode=5000               !maximum of modes of Bloch waves 

!>loop index
integer :: i,j

!>parameters of the incident light 
real(8) :: f,minf,maxf,hf                      !frequency 
real(8) :: minlam,maxlam,hlam,lam              !wavelength 
real(8) :: theta0_deg,theta0_rad               !zenith angle 
real(8) :: alpha0_deg,alpha0_rad               !azimuth angle 

!>parameters of the vacuum medium 
real(8) :: k0                   			   !wavenumber of the vacuum            
real(8) :: kx0,ky0              		       !tangential components of incident wavevector in the vacuum 
real(8) :: kz0                  			   !normal component of the incident wavevector in the vacuum  

!>parameters of the input medium 
real(8) :: epsin       						   !dielectric constant of the input medium 
real(8) :: nin         						   !refractive index of the input medium 
real(8) :: kin         						   !wavenumber of the input medium 
real(8) :: kzin        						   !normal component of the incident wavevector in the input medium 
real(8) :: y0in_p,y0in_s 
real(8) :: y0in

!>parameters of the output medium 
complex(8) :: epsout        !dielectric constant of the output medium 
complex(8) :: nout          !refractive index of the output medium 
complex(8) :: kout          !wavenumber of the output medium 
complex(8) :: kzout         !normal component of the incident wavevector in the output medium 
complex(8) :: y0out         !admittance : y0out = kout/kzout
complex(8) :: y0out_p,y0out_s

!>parameters of a dielectric layer
complex(8) :: epsd        !dielectric constant    
complex(8) :: nd          !refractive index
complex(8) :: kd          !wavenumber
complex(8) :: kzd         !normal component of wavevector 
complex(8) :: y0d_p,y0d_s !0th-order admittance 
real(8)    :: d           !thickness

!>parameters of the square hole array  
real(8) :: l             	  !period 
real(8) :: a             	  !hole width  
real(8) :: h             	  !hole height 
complex(8) :: epsh       	  !dielectric constant of the hole's waveguide 
complex(8) :: muh        	  !magnetic permeability of the hole's waveguide 
complex(8) :: nh         	  !refractive index of the hole's waveguide
complex(8) :: kh         	  !wavenumber of the hole's waveguide
complex(8) :: epsm            !dielectric constant of the metal 
real(8)    :: epsm_re,epsm_im 
complex(8) :: nm              !refractive index od the metal 
real(8)    :: nm_re,nm_im
complex(8) :: zsp             !surface impedance of the metal
complex(8) :: km              !wavenumber in the metal
complex(8) :: kzm             !z-component of wavevector in the metal 
complex(8) :: y0m_p,y0m_s    !admittance of the metal 
complex(8) :: ein,eout 
complex(8) :: gin,gout 
complex(8) :: gv
complex(8) :: gam
complex(8) :: i0 
complex(8) :: disp           

!>total transmission coefficients
complex(8) :: a0_out,an_out 

!>total reflection coefficients 
complex(8) :: b0_in,bn_in 
complex(8) :: b0_in_nd,b0_in_d           !non-diffracted and diffracted terms of b0_in

!>total reflectance and transmittance 
real(8) :: refl0,trans0
real(8) :: refln,transn

!>total absorbance 
real(8) :: absorpn,absorp0  

!>lh = ly-h 
real(8) :: lh 

!>Bloch wavevectors and related parameters 
integer    :: m,mfar_in,mfar_out   !Bloch wavenumber along the x-axis
integer    :: n,nfar_in,nfar_out   !Bloch wavenumber along the y-axis
real(8)    :: betam                !kx0+m*2*pi/l
real(8)    :: betan                !ky0+n*2*pi/l 
real(8)    :: betamn                  
real(8)    :: qmnin                
complex(8) :: qmnout
complex(8) :: qmnd 
complex(8) :: qmnm 
real(8)    :: ymnin_p,ymnin_s      !admittance in the input 
complex(8) :: ymnout_p,ymnout_s    !admittance in the output
complex(8) :: ymnd_p,ymnd_s        !admittance in the dielectric medium
complex(8) :: ymnm_p,ymnm_s        !admittance in the metal medium
integer    :: maxn_in              !maximum order needed for the convergence of ginf_s_in 
integer    :: maxn_out             !maximum order needed for the convergence of ginf_s_out 

!>0th-order reflection coefficients between two regions
complex(8) :: rho_d_in_p
complex(8) :: rho_m_d_p  
  
!>0th-order transmission coefficients between two regions
complex(8) :: tau_d_in_p 
complex(8) :: tau_m_d_p  

!>dot product <betam,betan|TE10>
complex(8) :: dot,dot_p,dot_s

interface 
  
  complex(8) function dot_bloch_te01(l,a,betam,betan)
    implicit none
    real(8), intent(in) :: l
    real(8), intent(in) :: a
    real(8), intent(in) :: betam 
    real(8), intent(in) :: betan 
  end function dot_bloch_te01
     
  subroutine gin_1layer_sub(nmode,k0,kx0,ky0,epsd,kd,d,l,a,epsin,kin,epsm,km,gin,maxn_in)
    implicit none
    integer,    intent(in) :: nmode
    real(8),    intent(in) :: k0 
    real(8),    intent(in) :: kx0
    real(8),    intent(in) :: ky0
    complex(8), intent(in) :: epsd
    complex(8), intent(in) :: kd
    real(8),    intent(in) :: d
    real(8),    intent(in) :: l
    real(8),    intent(in) :: a
    real(8),    intent(in) :: epsin
    real(8),    intent(in) :: kin
    complex(8), intent(in) :: epsm
    complex(8), intent(in) :: km
    complex(8), intent(out) :: gin
    integer,    intent(out) :: maxn_in 
  end subroutine gin_1layer_sub
  
  subroutine gout_0layer_sub(nmode,k0,kx0,ky0,epsout,nout,kout,l,a,epsm,km,gout,maxn_out)
    implicit none
	integer,    intent(in) :: nmode
	real(8),    intent(in) :: k0
	real(8),    intent(in) :: kx0
	real(8),    intent(in) :: ky0
	complex(8), intent(in) :: epsout
	complex(8), intent(in) :: nout
	complex(8), intent(in) :: kout 
	real(8),    intent(in) :: l
	real(8),    intent(in) :: a
	complex(8), intent(in) :: epsm 
	complex(8), intent(in) :: km 
	complex(8), intent(out) :: gout
	integer,    intent(out) :: maxn_out
  end subroutine gout_0layer_sub 
	  
  complex(8) function gv_func(kh,nh,muh,h,a)
    implicit none
    complex(8), intent(in) :: kh
    complex(8), intent(in) :: nh
    complex(8), intent(in) :: muh 
    real(8),    intent(in) :: h 
    real(8),    intent(in) :: a 
  end function gv_func 
  
  complex(8) function gam_func(kh,nh,muh,h,a)
    implicit none
    complex(8), intent(in) :: kh
    complex(8), intent(in) :: nh
    complex(8), intent(in) :: muh 
    real(8),    intent(in) :: h 
    real(8),    intent(in) :: a 
  end function gam_func 

  complex(8) function i0_func(epsin,nin,kin,kx0,ky0,l,a,yin,ym,yd,kzd,d,kzin)
    implicit none
    real(8), intent(in) :: epsin
    real(8), intent(in) :: nin
    real(8), intent(in) :: kin
    real(8), intent(in) :: kx0
    real(8), intent(in) :: ky0 
    real(8), intent(in) :: l
    real(8), intent(in) :: a
    real(8), intent(in)    :: yin 
    complex(8), intent(in) :: ym
    complex(8), intent(in) :: yd  
    complex(8), intent(in) :: kzd
    real(8), intent(in) :: d
    real(8), intent(in) :: kzin
  end function i0_func
    
  complex(8) function eps_Ag(f0)
    real(8), intent(in) :: f0
  end function eps_Ag
  
  complex(8) function eps_Ag_Ordal(f0)
    real(8), intent(in) :: f0
  end function eps_Ag_Ordal
  
  complex(8) function eps_Ni(f0)
    real(8), intent(in) :: f0
  end function eps_Ni
  
end interface

!>ouput files 
open(7,file="spectra-sensor-p300mic-a202mic-h5.1mic-d0.1mic-theta0deg.dat",status="replace")

!>define the geometry of the square hole array (m) 
l = 300.d-6
a = 202d-6
h = 5.1d-6 

!>define optical properties of material filling the hole
!>Polystyrene (PS)
epsh = dcmplx(2.59d0,0.0161d0)
muh  = dcmplx(1.d0,0.d0)
nh   = sqrt(dcmplx(epsh*muh))       

!>define the parameters of the input medium 
!>air
epsin = 1.d0 
nin   = sqrt(epsin)

!>define the parameters of the dielectric layer
!>Polystyrene (PS)
epsd = dcmplx(2.59d0,0.0161d0)
nd   = sqrt(dcmplx(epsd))
d    = 0.1d-6

!>define the wavelength  
hlam   = 0.01d-6
minlam = 301.d-6
maxlam = 350.d-6

!>define the frequency range 
hf   = 0.0000001d+12 
minf = 0.998d+12 
maxf = 1.d+12

!>zenith angle 
theta0_deg = 0.d0
!write(*,*) "Enter the zenith angle (degree) : "
!read (*,*) theta0_deg
theta0_rad = theta0_deg*pi/180.d0

!>azimuth angle 
alpha0_deg = 0.d0
alpha0_rad = alpha0_deg*pi/180.d0 

f = minf
do 
  if (f>maxf) exit  
  
  !>frequency 
  !f = c/lam
  
  !>wavelength 
  lam = c/f 
  
  !>refractive index/dielectric constant of the metal 
  epsm = eps_Ag_Ordal(f) 
  nm   = sqrt(dcmplx(epsm))
  
  !>refractive index/dielectric constant of the output medium
  epsout = epsm
  nout   = sqrt(dcmplx(epsout)) 
  
  !>wavenumber in free space 
  !k0 = 2.d0*pi*f/c  
  k0 = 2.d0*pi/lam
  
  !>kx0 and ky0
  kx0 = k0*sin(theta0_rad)*cos(alpha0_rad)
  ky0 = k0*sin(theta0_rad)*sin(alpha0_rad)
  
  !>check the propagating condition 
  if (k0>abs(sqrt(kx0**2 + ky0**2))) then

    !>wavenumbers 
    kin  = nin*k0 
    kout = nout*k0
    kh   = nh*k0
    kd   = nd*k0
    km   = nm*k0
    
    !>wavevectors  
    kzin  = sqrt(kin**2 - kx0**2 - ky0**2)
    kzout = sqrt(dcmplx(kout**2 - kx0**2 - ky0**2))
    kzd   = sqrt(dcmplx(kd**2 - kx0**2 - ky0**2)) 
    kzm   = sqrt(dcmplx(km**2 - kx0**2 - ky0**2)) 
    
    !>admittance 
    y0in_p  = (epsin*k0+1.d-100)/(kzin+1.d-100)
    y0out_p = (epsout*k0+1.d-100)/(kzout+1.d-100)
    y0d_p   = (epsd*k0+1.d-100)/(kzd+1.d-100)
    y0m_p   = (epsm*k0+1.d-100)/(kzm+1.d-100) 
    
    !>coupled-mode analysis parameters
    call gin_1layer_sub(nmode,k0,kx0,ky0,epsd,kd,d,l,a,epsin,kin,epsm,km,gin,maxn_in) 
    call  gout_0layer_sub(nmode,k0,kx0,ky0,epsout,nout,kout,l,a,epsm,km,gout,maxn_out)
    gv  = gv_func(kh,nh,muh,h,a)
    gam = gam_func(kh,nh,muh,h,a)
    i0  = i0_func(epsin,nin,kin,kx0,ky0,l,a,y0in_p,y0m_p,y0d_p,kzd,d,kzin)  
    
    !>obtain dispersion function 
    disp = (gin-gam)*(gout-gam) - gv**2
    
    !>input and output field coefficients of the hole 
    ein  = i0*(gout-gam)/disp
    eout = i0*gv/disp
    
    !>dot product 
    dot = dot_bloch_te01(l,a,kx0,ky0)
    dot_p=dot
    dot_s=0.d0  
      
    !>0th-order reflection coefficients between two regions
    rho_d_in_p = (y0d_p-y0in_p)/(y0d_p+y0in_p) 
    rho_m_d_p  = (y0m_p-y0d_p)/(y0m_p+y0d_p)
  
    !>0th-order transmission coefficients between two regions
    tau_d_in_p = 2.d0*y0d_p/(y0d_p+y0in_p) 
    tau_m_d_p  = 2.d0*y0m_p/(y0m_p+y0d_p) 
  
    !>***************************************
    !>0th-order total reflection coefficient 
    !>***************************************
    !>non-diffracted term
    b0_in_nd = -((rho_d_in_p+rho_m_d_p*exp(2.d0*iu*d*kzd))/ &
                 (1.d0+rho_m_d_p*rho_d_in_p*exp(2.d0*iu*d*kzd))) 
  
    !>diffracted term   
    b0_in_d  = 0.5d0*tau_m_d_p*tau_d_in_p*exp(iu*d*kzd)*ein*dot_p/(1.d0+rho_m_d_p*rho_d_in_p*exp(2.d0*iu*d*kzd)) 
       
    b0_in = b0_in_nd + b0_in_d
  
    !>0th-order reflectance
    refl0 = abs(b0_in)**2 
  
    !>0th-order ansorbance 
    absorp0 = 1.d0 - refl0 

    !>write file
    write(*,"(4F15.7)") f/1.d+12,lam/1.d-6,refl0,absorp0
    write(7,"(4F15.7)") f/1.d+12,lam/1.d-6,refl0,absorp0

  end if !>end-if of the k0>abs(kx0) condition

  !>update spectral point
  f = f + hf
  
end do

close(7)
end program 

 complex(8) function dot_bloch_te01(l,a,betam,betan)
implicit none
!>
!><betamn,betan|TE01>
!>
real(8), intent(in) :: l
real(8), intent(in) :: a
real(8), intent(in) :: betam 
real(8), intent(in) :: betan 

!>local variables
real(8),    parameter :: pi = 3.d0*acos(0.5d0) 
real(8) :: sincy_minus    !sinc((betan-pi/a)*a/2) 
real(8) :: sincy_plus     !sinc((betan+pi/a)*a/2) 
real(8) :: sincx 

!>obtain the sincy function 
if (abs((betan-pi/a)*a*0.5d0)<1.d-200) then
  sincy_minus=1.d0
else
  sincy_minus=sin((betan-pi/a)*a*0.5d0)/((betan-pi/a)*a*0.5d0)
end if 

if (abs((betan+pi/a)*a*0.5d0)<1.d-200) then
  sincy_plus=1.d0
else
  sincy_plus=sin((betan+pi/a)*a*0.5d0)/((betan+pi/a)*a*0.5d0)
end if 

!>obtain sincx
if (abs(betam*a*0.5d0)<1.d-200) then
  sincx=1.d0
else
  sincx=sin(betam*a*0.5d0)/(betam*a*0.5d0)
end if 

dot_bloch_te01=-(a/l)*(1.d0/sqrt(2.d0))*(sincy_minus+sincy_plus)*sincx
 
return
end function dot_bloch_te01

subroutine gin_1layer_sub(nmode,k0,kx0,ky0,epsd,kd,d,l,a,epsin,kin,epsm,km,gin,maxn_in)
implicit none
integer,    intent(in) :: nmode
real(8),    intent(in) :: k0 
real(8),    intent(in) :: kx0
real(8),    intent(in) :: ky0
complex(8), intent(in) :: epsd
complex(8), intent(in) :: kd
real(8),    intent(in) :: d
real(8),    intent(in) :: l
real(8),    intent(in) :: a
real(8),    intent(in) :: epsin
real(8),    intent(in) :: kin
complex(8), intent(in) :: epsm
complex(8), intent(in) :: km
complex(8), intent(out) :: gin
integer,    intent(out) :: maxn_in 

!>local variables
real(8),    parameter :: pi = 3.d0*acos(0.5d0) 
complex(8), parameter :: iu = (0.d0,1.d0)
integer    :: m,n,mm,nn
real(8)    :: betam           	  !kx0 + m*2*pi/l 
real(8)    :: betan     	  	  !ky0 + n*2*pi/l 
complex(8) :: qmnd      	  	  !sqrt(kd^2 - (betam^2+betan^2))
complex(8) :: qmnin          	  !sqrt(kin^2 - (betam^2+betan^2))
complex(8) :: qmnm            	  !sqrt(km^2 - (betam^2+betan^2))
complex(8) :: ymnd_p,ymnd_s   	  !admittance in dielectric layer
complex(8) :: ymnin_p,ymnin_s     !admittance in input medium 
complex(8) :: ymnm_p,ymnm_s   	  !admittance in metal layer 
complex(8) :: ginf,ginf_old 
complex(8) :: ginf_plus,ginf_minus
real(8)    :: ginf_err,ginf_err_accpt
complex(8) :: dot,dot_p,dot_s     !<betam,betan|TE10>

!>reflection coefficients between two regions
complex(8) :: rho_m_d_p,rho_m_d_s  
complex(8) :: rho_d_in_p,rho_d_in_s

!>multiple reflection parameters  
complex(8) :: fm_p,fm_s  
    
interface
 
  complex(8) function dot_bloch_te01(l,a,betam,betan)
    implicit none
    real(8), intent(in) :: l
    real(8), intent(in) :: a
    real(8), intent(in) :: betam 
    real(8), intent(in) :: betan 
  end function dot_bloch_te01

end interface

!>define accepted error (%) 
ginf_err_accpt = 0.001d0

n=0
do 
  if (n>nmode) then
    write(*,*) "order is too large for gin.."
    stop
  end if
  
  if (n==0) then
    m=0
    !>allowed diffraction wavenumber (1/m) 
    betam  = kx0 + m*(2.d0*pi/l) 
    betan  = ky0 + n*(2.d0*pi/l) 
    qmnd   = sqrt(dcmplx(kd**2 - betam**2 - betan**2))
    qmnin  = sqrt(dcmplx(kin**2 - betam**2 - betan**2))
    qmnm   = sqrt(dcmplx(km**2 - betam**2 - betan**2))
    
    !>admittance 
    ymnd_p = (epsd*k0+1.d-100)/(qmnd+1.d-100)
    ymnd_s = qmnd/k0 
    
    ymnin_p = (epsin*k0+1.d-100)/(qmnin+1.d-100)
    ymnin_s = qmnin/k0 
     
    ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100)
    ymnm_s = qmnm/k0 
                 
    !>obtain <betam,betan|TE10> 
    dot=dot_bloch_te01(l,a,betam,betan)
    if ((m==0).and.(n==0)) then
      dot_p = dot
      dot_s = 0.d0
    else
      dot_p = (betam/sqrt(betam**2+betan**2))*dot
      dot_s = (betan/sqrt(betam**2+betan**2))*dot 
    end if
    
    !>reflection coefficients between two regions
    rho_m_d_p  = (ymnm_p-ymnd_p)/(ymnm_p+ymnd_p)
    rho_m_d_s  = (ymnm_s-ymnd_s)/(ymnm_s+ymnd_s)
    
    rho_d_in_p = (ymnd_p-ymnin_p)/(ymnd_p+ymnin_p)
    rho_d_in_s = (ymnd_s-ymnin_s)/(ymnd_s+ymnin_s) 
    
    !>multiple reflection parameters    
    fm_p = (1.d0-rho_d_in_p*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_p*rho_d_in_p*exp(2.d0*iu*qmnd*d))  
    fm_s = (1.d0-rho_d_in_s*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_s*rho_d_in_s*exp(2.d0*iu*qmnd*d))
    
    !>update the ginf
    ginf_old = iu*((ymnd_p*ymnm_p)/(ymnd_p+ymnm_p))*fm_p*(abs(dot_p)**2) + & 
               iu*((ymnd_s*ymnm_s)/(ymnd_s+ymnm_s))*fm_s*(abs(dot_s)**2) 
    
  else if (n==1) then
  
    !>*****************
    !>positive order n  
    !>***************** 
    ginf_plus = 0.d0
    do m=-n,n,1
      !>allowed diffraction wavenumber (1/m) 
      betam = kx0 + m*(2.d0*pi/l)
      betan = ky0 + n*(2.d0*pi/l)  
      qmnd   = sqrt(dcmplx(kd**2 - betam**2 - betan**2))
      qmnin  = sqrt(dcmplx(kin**2 - betam**2 - betan**2))
      qmnm   = sqrt(dcmplx(km**2 - betam**2 - betan**2))
    
      !>admittance 
      ymnd_p = (epsd*k0+1.d-100)/(qmnd+1.d-100)
      ymnd_s = qmnd/k0 
    
      ymnin_p = (epsin*k0+1.d-100)/(qmnin+1.d-100)
      ymnin_s = qmnin/k0 
     
      ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100)
      ymnm_s = qmnm/k0 
                 
      !>obtain <betam,betan|TE10> 
      dot=dot_bloch_te01(l,a,betam,betan)
      if ((m==0).and.(n==0)) then
        dot_p = dot
        dot_s = 0.d0
      else
        dot_p = (betam/sqrt(betam**2+betan**2))*dot
        dot_s = (betan/sqrt(betam**2+betan**2))*dot 
      end if
    
      !>reflection coefficients between two regions
      rho_m_d_p  = (ymnm_p-ymnd_p)/(ymnm_p+ymnd_p)
      rho_m_d_s  = (ymnm_s-ymnd_s)/(ymnm_s+ymnd_s)
    
      rho_d_in_p = (ymnd_p-ymnin_p)/(ymnd_p+ymnin_p)
      rho_d_in_s = (ymnd_s-ymnin_s)/(ymnd_s+ymnin_s) 
    
      !>multiple reflection parameters    
      fm_p = (1.d0-rho_d_in_p*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_p*rho_d_in_p*exp(2.d0*iu*qmnd*d))  
      fm_s = (1.d0-rho_d_in_s*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_s*rho_d_in_s*exp(2.d0*iu*qmnd*d))
       
      !>update the ginf
      ginf_plus = ginf_plus + iu*((ymnd_p*ymnm_p)/(ymnd_p+ymnm_p))*fm_p*(abs(dot_p)**2) + & 
                              iu*((ymnd_s*ymnm_s)/(ymnd_s+ymnm_s))*fm_s*(abs(dot_s)**2)  
    
    end do
    
    !>*******************
    !>negative order n
    !>*******************
    ginf_minus=0.d0
    do m=-n,n,1
      !>allowed diffraction wavenumber (1/m) 
      betam = kx0 + m*(2.d0*pi/l)
      betan = ky0 - n*(2.d0*pi/l)
      qmnd   = sqrt(dcmplx(kd**2 - betam**2 - betan**2))
      qmnin  = sqrt(dcmplx(kin**2 - betam**2 - betan**2))
      qmnm   = sqrt(dcmplx(km**2 - betam**2 - betan**2))
    
      !>admittance 
      ymnd_p = (epsd*k0+1.d-100)/(qmnd+1.d-100)
      ymnd_s = qmnd/k0 
    
      ymnin_p = (epsin*k0+1.d-100)/(qmnin+1.d-100)
      ymnin_s = qmnin/k0 
     
      ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100)
      ymnm_s = qmnm/k0 
                 
      !>obtain <betam,betan|TE10> 
      dot=dot_bloch_te01(l,a,betam,betan)
      if ((m==0).and.(n==0)) then
        dot_p = dot
        dot_s = 0.d0
      else
        dot_p = (betam/sqrt(betam**2+betan**2))*dot
        dot_s = (betan/sqrt(betam**2+betan**2))*dot 
      end if
    
      !>reflection coefficients between two regions
      rho_m_d_p  = (ymnm_p-ymnd_p)/(ymnm_p+ymnd_p)
      rho_m_d_s  = (ymnm_s-ymnd_s)/(ymnm_s+ymnd_s)
    
      rho_d_in_p = (ymnd_p-ymnin_p)/(ymnd_p+ymnin_p)
      rho_d_in_s = (ymnd_s-ymnin_s)/(ymnd_s+ymnin_s) 
    
      !>multiple reflection parameters    
      fm_p = (1.d0-rho_d_in_p*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_p*rho_d_in_p*exp(2.d0*iu*qmnd*d))  
      fm_s = (1.d0-rho_d_in_s*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_s*rho_d_in_s*exp(2.d0*iu*qmnd*d))
       
      !>update the ginf
      ginf_minus = ginf_minus + iu*((ymnd_p*ymnm_p)/(ymnd_p+ymnm_p))*fm_p*(abs(dot_p)**2) + & 
                                iu*((ymnd_s*ymnm_s)/(ymnd_s+ymnm_s))*fm_s*(abs(dot_s)**2)  
                              
    end do
                                        
    !>***********************
    !>update total ginf_old
    !>***********************
    ginf_old = ginf_old + ginf_plus + ginf_minus
    
    !>*******************
    !>m=n,nn=-(n-1)->(n-1)
    !>*******************
    m=n
    nn=0 
    !>allowed diffraction wavenumber (1/m) 
    betam = kx0 + m*(2.d0*pi/l)
    betan = ky0 + nn*(2.d0*pi/l)
    qmnd   = sqrt(dcmplx(kd**2 - betam**2 - betan**2))
    qmnin  = sqrt(dcmplx(kin**2 - betam**2 - betan**2))
    qmnm   = sqrt(dcmplx(km**2 - betam**2 - betan**2))
    
    !>admittance 
    ymnd_p = (epsd*k0+1.d-100)/(qmnd+1.d-100)
    ymnd_s = qmnd/k0 
    
    ymnin_p = (epsin*k0+1.d-100)/(qmnin+1.d-100)
    ymnin_s = qmnin/k0 
     
    ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100)
    ymnm_s = qmnm/k0 
                 
    !>obtain <betam,betan|TE10> 
    dot=dot_bloch_te01(l,a,betam,betan)
    if ((m==0).and.(n==0)) then
      dot_p = dot
      dot_s = 0.d0
    else
      dot_p = (betam/sqrt(betam**2+betan**2))*dot
      dot_s = (betan/sqrt(betam**2+betan**2))*dot 
    end if
    
    !>reflection coefficients between two regions
    rho_m_d_p  = (ymnm_p-ymnd_p)/(ymnm_p+ymnd_p)
    rho_m_d_s  = (ymnm_s-ymnd_s)/(ymnm_s+ymnd_s)
    
    rho_d_in_p = (ymnd_p-ymnin_p)/(ymnd_p+ymnin_p)
    rho_d_in_s = (ymnd_s-ymnin_s)/(ymnd_s+ymnin_s) 
    
    !>multiple reflection parameters    
    fm_p = (1.d0-rho_d_in_p*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_p*rho_d_in_p*exp(2.d0*iu*qmnd*d))  
    fm_s = (1.d0-rho_d_in_s*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_s*rho_d_in_s*exp(2.d0*iu*qmnd*d))
       
    !>update the ginf
    ginf_old = ginf_old + iu*((ymnd_p*ymnm_p)/(ymnd_p+ymnm_p))*fm_p*(abs(dot_p)**2) + & 
                          iu*((ymnd_s*ymnm_s)/(ymnd_s+ymnm_s))*fm_s*(abs(dot_s)**2)  
      
    !>*******************
    !>m=-n,nn=-(n-1)->(n-1)
    !>*******************
    m=-n
    nn=0
    !>allowed diffraction wavenumber (1/m) 
    betam = kx0 + m*(2.d0*pi/l)
    betan = ky0 + nn*(2.d0*pi/l)
    qmnd   = sqrt(dcmplx(kd**2 - betam**2 - betan**2))
    qmnin  = sqrt(dcmplx(kin**2 - betam**2 - betan**2))
    qmnm   = sqrt(dcmplx(km**2 - betam**2 - betan**2))
    
    !>admittance 
    ymnd_p = (epsd*k0+1.d-100)/(qmnd+1.d-100)
    ymnd_s = qmnd/k0 
    
    ymnin_p = (epsin*k0+1.d-100)/(qmnin+1.d-100)
    ymnin_s = qmnin/k0 
     
    ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100)
    ymnm_s = qmnm/k0 
                 
    !>obtain <betam,betan|TE10> 
    dot=dot_bloch_te01(l,a,betam,betan)
    if ((m==0).and.(n==0)) then
      dot_p = dot
      dot_s = 0.d0
    else
      dot_p = (betam/sqrt(betam**2+betan**2))*dot
      dot_s = (betan/sqrt(betam**2+betan**2))*dot 
    end if
   
    !>reflection coefficients between two regions
    rho_m_d_p  = (ymnm_p-ymnd_p)/(ymnm_p+ymnd_p)
    rho_m_d_s  = (ymnm_s-ymnd_s)/(ymnm_s+ymnd_s)
    
    rho_d_in_p = (ymnd_p-ymnin_p)/(ymnd_p+ymnin_p)
    rho_d_in_s = (ymnd_s-ymnin_s)/(ymnd_s+ymnin_s) 
    
    !>multiple reflection parameters    
    fm_p = (1.d0-rho_d_in_p*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_p*rho_d_in_p*exp(2.d0*iu*qmnd*d))  
    fm_s = (1.d0-rho_d_in_s*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_s*rho_d_in_s*exp(2.d0*iu*qmnd*d))
       
    !>update the ginf
    ginf_old = ginf_old + iu*((ymnd_p*ymnm_p)/(ymnd_p+ymnm_p))*fm_p*(abs(dot_p)**2) + & 
                          iu*((ymnd_s*ymnm_s)/(ymnd_s+ymnm_s))*fm_s*(abs(dot_s)**2)  
    
    
  else if (n>=2) then
    !>*****************
    !>positive order n  
    !>***************** 
    ginf_plus=0.d0 
    do m=-n,n,1 
      !>allowed diffraction wavenumber (1/m) 
      betam = kx0 + m*(2.d0*pi/l)
      betan = ky0 + n*(2.d0*pi/l)  
      qmnd   = sqrt(dcmplx(kd**2 - betam**2 - betan**2))
      qmnin  = sqrt(dcmplx(kin**2 - betam**2 - betan**2))
      qmnm   = sqrt(dcmplx(km**2 - betam**2 - betan**2))
    
      !>admittance 
      ymnd_p = (epsd*k0+1.d-100)/(qmnd+1.d-100)
      ymnd_s = qmnd/k0 
    
      ymnin_p = (epsin*k0+1.d-100)/(qmnin+1.d-100)
      ymnin_s = qmnin/k0 
     
      ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100)
      ymnm_s = qmnm/k0 
                 
      !>obtain <betam,betan|TE10> 
      dot=dot_bloch_te01(l,a,betam,betan)
      if ((m==0).and.(n==0)) then
        dot_p = dot
        dot_s = 0.d0
      else
        dot_p = (betam/sqrt(betam**2+betan**2))*dot
        dot_s = (betan/sqrt(betam**2+betan**2))*dot 
      end if
    
      !>reflection coefficients between two regions
      rho_m_d_p  = (ymnm_p-ymnd_p)/(ymnm_p+ymnd_p)
      rho_m_d_s  = (ymnm_s-ymnd_s)/(ymnm_s+ymnd_s)
    
      rho_d_in_p = (ymnd_p-ymnin_p)/(ymnd_p+ymnin_p)
      rho_d_in_s = (ymnd_s-ymnin_s)/(ymnd_s+ymnin_s) 
    
      !>multiple reflection parameters    
      fm_p = (1.d0-rho_d_in_p*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_p*rho_d_in_p*exp(2.d0*iu*qmnd*d))  
      fm_s = (1.d0-rho_d_in_s*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_s*rho_d_in_s*exp(2.d0*iu*qmnd*d))
       
      !>update the ginf
      ginf_plus = ginf_plus + iu*((ymnd_p*ymnm_p)/(ymnd_p+ymnm_p))*fm_p*(abs(dot_p)**2) + & 
                              iu*((ymnd_s*ymnm_s)/(ymnd_s+ymnm_s))*fm_s*(abs(dot_s)**2)  
         
    end do !>end-of n loop
    
    !>*******************
    !>negative order n
    !>*******************
    ginf_minus=0.d0
    do m=-n,n,1
      !>allowed diffraction wavenumber (1/m) 
      betam = kx0 + m*(2.d0*pi/l)
      betan = ky0 - n*(2.d0*pi/l)
      qmnd   = sqrt(dcmplx(kd**2 - betam**2 - betan**2))
      qmnin  = sqrt(dcmplx(kin**2 - betam**2 - betan**2))
      qmnm   = sqrt(dcmplx(km**2 - betam**2 - betan**2))
    
      !>admittance 
      ymnd_p = (epsd*k0+1.d-100)/(qmnd+1.d-100)
      ymnd_s = qmnd/k0 
    
      ymnin_p = (epsin*k0+1.d-100)/(qmnin+1.d-100)
      ymnin_s = qmnin/k0 
     
      ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100)
      ymnm_s = qmnm/k0 
                 
      !>obtain <betam,betan|TE10> 
      dot=dot_bloch_te01(l,a,betam,betan)
      if ((m==0).and.(n==0)) then
        dot_p = dot
        dot_s = 0.d0
      else
        dot_p = (betam/sqrt(betam**2+betan**2))*dot
        dot_s = (betan/sqrt(betam**2+betan**2))*dot 
      end if
    
      !>reflection coefficients between two regions
      rho_m_d_p  = (ymnm_p-ymnd_p)/(ymnm_p+ymnd_p)
      rho_m_d_s  = (ymnm_s-ymnd_s)/(ymnm_s+ymnd_s)
    
      rho_d_in_p = (ymnd_p-ymnin_p)/(ymnd_p+ymnin_p)
      rho_d_in_s = (ymnd_s-ymnin_s)/(ymnd_s+ymnin_s) 
    
      !>multiple reflection parameters    
      fm_p = (1.d0-rho_d_in_p*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_p*rho_d_in_p*exp(2.d0*iu*qmnd*d))  
      fm_s = (1.d0-rho_d_in_s*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_s*rho_d_in_s*exp(2.d0*iu*qmnd*d))
       
      !>update the ginf
      ginf_minus = ginf_minus + iu*((ymnd_p*ymnm_p)/(ymnd_p+ymnm_p))*fm_p*(abs(dot_p)**2) + & 
                                iu*((ymnd_s*ymnm_s)/(ymnd_s+ymnm_s))*fm_s*(abs(dot_s)**2)  
                              
    end do !>end-do of n loop
    
    !>***********************
    !>update total ginf_old
    !>***********************
    ginf = ginf_old + ginf_plus + ginf_minus
    
    !>*******************
    !>m=n,nn=-(n-1)->(n-1)
    !>*******************
    m=n
    do nn=-(n-1),n-1,1 
      !>allowed diffraction wavenumber (1/m) 
      betam = kx0 + m*(2.d0*pi/l)
      betan = ky0 + nn*(2.d0*pi/l)
      qmnd   = sqrt(dcmplx(kd**2 - betam**2 - betan**2))
      qmnin  = sqrt(dcmplx(kin**2 - betam**2 - betan**2))
      qmnm   = sqrt(dcmplx(km**2 - betam**2 - betan**2))
    
      !>admittance 
      ymnd_p = (epsd*k0+1.d-100)/(qmnd+1.d-100)
      ymnd_s = qmnd/k0 
    
      ymnin_p = (epsin*k0+1.d-100)/(qmnin+1.d-100)
      ymnin_s = qmnin/k0 
     
      ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100)
      ymnm_s = qmnm/k0 
                 
      !>obtain <betam,betan|TE10> 
      dot=dot_bloch_te01(l,a,betam,betan)
      if ((m==0).and.(n==0)) then
        dot_p = dot
        dot_s = 0.d0
      else
        dot_p = (betam/sqrt(betam**2+betan**2))*dot
        dot_s = (betan/sqrt(betam**2+betan**2))*dot 
      end if
    
      !>reflection coefficients between two regions
      rho_m_d_p  = (ymnm_p-ymnd_p)/(ymnm_p+ymnd_p)
      rho_m_d_s  = (ymnm_s-ymnd_s)/(ymnm_s+ymnd_s)
    
      rho_d_in_p = (ymnd_p-ymnin_p)/(ymnd_p+ymnin_p)
      rho_d_in_s = (ymnd_s-ymnin_s)/(ymnd_s+ymnin_s) 
    
      !>multiple reflection parameters    
      fm_p = (1.d0-rho_d_in_p*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_p*rho_d_in_p*exp(2.d0*iu*qmnd*d))  
      fm_s = (1.d0-rho_d_in_s*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_s*rho_d_in_s*exp(2.d0*iu*qmnd*d))
       
      !>update the ginf
      ginf = ginf + iu*((ymnd_p*ymnm_p)/(ymnd_p+ymnm_p))*fm_p*(abs(dot_p)**2) + & 
                    iu*((ymnd_s*ymnm_s)/(ymnd_s+ymnm_s))*fm_s*(abs(dot_s)**2)  
                              
    end do
      
    !>*******************
    !>m=-n,nn=-(n-1)->(n-1)
    !>*******************
    m=-n
    do nn=-(n-1),n-1,1
      !>allowed diffraction wavenumber (1/m) 
      betam = kx0 + m*(2.d0*pi/l)
      betan = ky0 + nn*(2.d0*pi/l)
      qmnd   = sqrt(dcmplx(kd**2 - betam**2 - betan**2))
      qmnin  = sqrt(dcmplx(kin**2 - betam**2 - betan**2))
      qmnm   = sqrt(dcmplx(km**2 - betam**2 - betan**2))
    
      !>admittance 
      ymnd_p = (epsd*k0+1.d-100)/(qmnd+1.d-100)
      ymnd_s = qmnd/k0 
    
      ymnin_p = (epsin*k0+1.d-100)/(qmnin+1.d-100)
      ymnin_s = qmnin/k0 
     
      ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100)
      ymnm_s = qmnm/k0 
                 
      !>obtain <betam,betan|TE10> 
      dot=dot_bloch_te01(l,a,betam,betan)
      if ((m==0).and.(n==0)) then
        dot_p = dot
        dot_s = 0.d0
      else
        dot_p = (betam/sqrt(betam**2+betan**2))*dot
        dot_s = (betan/sqrt(betam**2+betan**2))*dot 
      end if
    
      !>reflection coefficients between two regions
      rho_m_d_p  = (ymnm_p-ymnd_p)/(ymnm_p+ymnd_p)
      rho_m_d_s  = (ymnm_s-ymnd_s)/(ymnm_s+ymnd_s)
    
      rho_d_in_p = (ymnd_p-ymnin_p)/(ymnd_p+ymnin_p)
      rho_d_in_s = (ymnd_s-ymnin_s)/(ymnd_s+ymnin_s) 
    
      !>multiple reflection parameters    
      fm_p = (1.d0-rho_d_in_p*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_p*rho_d_in_p*exp(2.d0*iu*qmnd*d))  
      fm_s = (1.d0-rho_d_in_s*exp(2.d0*iu*qmnd*d))/(1.d0+rho_m_d_s*rho_d_in_s*exp(2.d0*iu*qmnd*d))
       
      !>update the ginf
      ginf = ginf + iu*((ymnd_p*ymnm_p)/(ymnd_p+ymnm_p))*fm_p*(abs(dot_p)**2) + & 
                    iu*((ymnd_s*ymnm_s)/(ymnd_s+ymnm_s))*fm_s*(abs(dot_s)**2)  
                              
    end do

    !>*************************
    !>check the convergence of ginf
    !>*************************
    ginf_err = ((abs(ginf-ginf_old)+1.d-200)/(abs(ginf_old)+1.d-200))*100.d0
      
    if (ginf_err<ginf_err_accpt) then
      gin = ginf
      maxn_in = n  
      return
    else
      ginf_old = ginf
    end if
  end if !>end-if of n>=1 loop
 
  !>update n
  n = n+1
end do !>end-do of n loop

end subroutine gin_1layer_sub 

subroutine gout_0layer_sub(nmode,k0,kx0,ky0,epsout,nout,kout,l,a,epsm,km,gout,maxn_out)
implicit none
integer,    intent(in) :: nmode
real(8),    intent(in) :: k0
real(8),    intent(in) :: kx0
real(8),    intent(in) :: ky0
complex(8), intent(in) :: epsout
complex(8), intent(in) :: nout
complex(8), intent(in) :: kout 
real(8),    intent(in) :: l
real(8),    intent(in) :: a
complex(8), intent(in) :: epsm 
complex(8), intent(in) :: km 
complex(8), intent(out) :: gout
integer,    intent(out) :: maxn_out

!>local variables
real(8),    parameter :: pi = 3.d0*acos(0.5d0) 
complex(8), parameter :: iu = (0.d0,1.d0)
integer    :: m,n,mm,nn
real(8)    :: betam     		!kx0 + m*2*pi/l 
real(8)    :: betan     		!ky0 + n*2*pi/l 
complex(8) :: qmnout     		!sqrt(kout^2 - (betam^2+betan^2))
complex(8) :: qmnm              !sqrt(km^2 - (betam^2+betan^2))
complex(8) :: ymnout_p,ymnout_s !admittance in output medium  
complex(8) :: ymnm_p,ymnm_s     !admittance in metal  
complex(8) :: ginf,ginf_old 
complex(8) :: ginf_plus,ginf_minus
real(8)    :: ginf_err,ginf_err_accpt
complex(8) :: dot,dot_p,dot_s          !<betam,betan|TE10>

interface
 
  complex(8) function dot_bloch_te01(l,a,betam,betan)
    implicit none
    real(8), intent(in) :: l
    real(8), intent(in) :: a
    real(8), intent(in) :: betam 
    real(8), intent(in) :: betan 
  end function dot_bloch_te01

end interface

!>define accepted error (%) 
ginf_err_accpt = 0.001d0

n=0
do 
  if (n>nmode) then
    write(*,*) "order is too large for ginf_s_in.."
    stop
  end if
  
  if (n==0) then
    m=0
    !>allowed diffraction wavenumber (1/m) 
    betam = kx0 + m*(2.d0*pi/l) 
    betan = ky0 + n*(2.d0*pi/l) 
    qmnout= sqrt(dcmplx(kout**2 - betam**2 -betan**2))
    qmnm  = sqrt(dcmplx(km**2 - betam**2 -betan**2))
    
    !>admittance
    ymnout_p = (epsout*k0+1.d-100)/(qmnout+1.d-100) 
    ymnout_s = qmnout/k0
    
    ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100) 
    ymnm_s = qmnm/k0 
    
    !>obtain <betam,betan|TE10> 
    dot=dot_bloch_te01(l,a,betam,betan)
    if ((m==0).and.(n==0)) then
      dot_p = dot
      dot_s = 0.d0
    else
      dot_p = (betam/sqrt(betam**2+betan**2))*dot
      dot_s = (betan/sqrt(betam**2+betan**2))*dot 
    end if
      
    !>update the ginf
    ginf_old = iu*((ymnout_p*ymnm_p+1.d-100)/(ymnm_p+ymnout_p+1.d-100))*(abs(dot_p)**2) + &
               iu*((ymnout_s*ymnm_s+1.d-100)/(ymnm_s+ymnout_s+1.d-100))*(abs(dot_s)**2)
     
  else if (n==1) then 
    !>*****************
    !>positive order n  
    !>*****************  
    ginf_plus = 0.d0
    do m=-n,n,1
      !>allowed diffraction wavenumber (1/m) 
      betam = kx0 + m*(2.d0*pi/l) 
      betan = ky0 + n*(2.d0*pi/l) 
      qmnout= sqrt(dcmplx(kout**2 - betam**2 -betan**2))
      qmnm  = sqrt(dcmplx(km**2 - betam**2 -betan**2))
    
      !>admittance
      ymnout_p = (epsout*k0+1.d-100)/(qmnout+1.d-100) 
      ymnout_s = qmnout/k0
    
      ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100) 
      ymnm_s = qmnm/k0 
     
      !>obtain <betam,betan|TE10> 
      dot=dot_bloch_te01(l,a,betam,betan)
      if ((m==0).and.(n==0)) then
        dot_p = dot
        dot_s = 0.d0
      else
        dot_p = (betam/sqrt(betam**2+betan**2))*dot
        dot_s = (betan/sqrt(betam**2+betan**2))*dot 
      end if
      
      !>update the ginf
      ginf_plus = ginf_plus + iu*((ymnout_p*ymnm_p+1.d-100)/(ymnm_p+ymnout_p+1.d-100))*(abs(dot_p)**2) + &
                              iu*((ymnout_s*ymnm_s+1.d-100)/(ymnm_s+ymnout_s+1.d-100))*(abs(dot_s)**2)   
    end do
      
    !>*******************
    !>negative order n
    !>*******************
    ginf_minus=0.d0
    do m=-n,n,1
      !>allowed diffraction wavenumber (1/m) 
      betam = kx0 + m*(2.d0*pi/l) 
      betan = ky0 - n*(2.d0*pi/l) 
      qmnout= sqrt(dcmplx(kout**2 - betam**2 -betan**2))
      qmnm  = sqrt(dcmplx(km**2 - betam**2 -betan**2))
    
      !>admittance
      ymnout_p = (epsout*k0+1.d-100)/(qmnout+1.d-100) 
      ymnout_s = qmnout/k0
    
      ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100) 
      ymnm_s = qmnm/k0 
     
      !>obtain <betam,betan|TE10> 
      dot=dot_bloch_te01(l,a,betam,betan)
      if ((m==0).and.(n==0)) then
        dot_p = dot
        dot_s = 0.d0
      else
        dot_p = (betam/sqrt(betam**2+betan**2))*dot
        dot_s = (betan/sqrt(betam**2+betan**2))*dot 
      end if
      
      !>update the ginf
      ginf_minus = ginf_minus + iu*((ymnout_p*ymnm_p+1.d-100)/(ymnm_p+ymnout_p+1.d-100))*(abs(dot_p)**2) + &
                                iu*((ymnout_s*ymnm_s+1.d-100)/(ymnm_s+ymnout_s+1.d-100))*(abs(dot_s)**2) 
    end do
    
    !>***********************
    !>update total ginf_old
    !>***********************
    ginf_old = ginf_old + ginf_plus + ginf_minus
    
    !>*******************
    !>m=n,nn=-(n-1)->(n-1)
    !>*******************
    m=n
    nn=0 
    !>allowed diffraction wavenumber (1/m) 
    betam = kx0 + m*(2.d0*pi/l)
    betan = ky0 + nn*(2.d0*pi/l)
    qmnout= sqrt(dcmplx(kout**2 - betam**2 -betan**2))
    qmnm  = sqrt(dcmplx(km**2 - betam**2 -betan**2))
    
    !>admittance
    ymnout_p = (epsout*k0+1.d-100)/(qmnout+1.d-100) 
    ymnout_s = qmnout/k0
    
    ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100) 
    ymnm_s = qmnm/k0 
     
    !>obtain <betam,betan|TE10> 
    dot=dot_bloch_te01(l,a,betam,betan)
    if ((m==0).and.(n==0)) then
      dot_p = dot
      dot_s = 0.d0
    else
      dot_p = (betam/sqrt(betam**2+betan**2))*dot
      dot_s = (betan/sqrt(betam**2+betan**2))*dot 
    end if
      
    !>update the ginf
    ginf_old = ginf_old + iu*((ymnout_p*ymnm_p+1.d-100)/(ymnm_p+ymnout_p+1.d-100))*(abs(dot_p)**2) + &
                          iu*((ymnout_s*ymnm_s+1.d-100)/(ymnm_s+ymnout_s+1.d-100))*(abs(dot_s)**2) 
      
    !>*******************
    !>m=-n,nn=-(n-1)->(n-1)
    !>*******************
    m=-n
    nn=0
    !>allowed diffraction wavenumber (1/m) 
    betam = kx0 + m*(2.d0*pi/l)
    betan = ky0 + nn*(2.d0*pi/l)
    qmnout= sqrt(dcmplx(kout**2 - betam**2 -betan**2))
    qmnm  = sqrt(dcmplx(km**2 - betam**2 -betan**2))
    
    !>admittance
    ymnout_p = (epsout*k0+1.d-100)/(qmnout+1.d-100) 
    ymnout_s = qmnout/k0
    
    ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100) 
    ymnm_s = qmnm/k0 
     
    !>obtain <betam,betan|TE10> 
    dot=dot_bloch_te01(l,a,betam,betan)
    if ((m==0).and.(n==0)) then
      dot_p = dot
      dot_s = 0.d0
    else
      dot_p = (betam/sqrt(betam**2+betan**2))*dot
      dot_s = (betan/sqrt(betam**2+betan**2))*dot 
    end if
      
    !>update the ginf
    ginf_old = ginf_old + iu*((ymnout_p*ymnm_p+1.d-100)/(ymnm_p+ymnout_p+1.d-100))*(abs(dot_p)**2) + &
                          iu*((ymnout_s*ymnm_s+1.d-100)/(ymnm_s+ymnout_s+1.d-100))*(abs(dot_s)**2) 
    
  else if (n>=2) then
    !>*****************
    !>positive order n
    !>***************** 
    ginf_plus=0.d0 
    do m=-n,n,1 
      !>allowed diffraction wavenumber (1/m) 
      betam = kx0 + m*(2.d0*pi/l)
      betan = ky0 + n*(2.d0*pi/l)  
      qmnout= sqrt(dcmplx(kout**2 - betam**2 -betan**2))
      qmnm  = sqrt(dcmplx(km**2 - betam**2 -betan**2))
    
      !>admittance
      ymnout_p = (epsout*k0+1.d-100)/(qmnout+1.d-100) 
      ymnout_s = qmnout/k0
    
      ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100) 
      ymnm_s = qmnm/k0 
     
      !>obtain <betam,betan|TE10> 
      dot=dot_bloch_te01(l,a,betam,betan)
      if ((m==0).and.(n==0)) then
        dot_p = dot
        dot_s = 0.d0
      else
        dot_p = (betam/sqrt(betam**2+betan**2))*dot
        dot_s = (betan/sqrt(betam**2+betan**2))*dot 
      end if
      
      !>update the ginf
      ginf_plus = ginf_plus + iu*((ymnout_p*ymnm_p+1.d-100)/(ymnm_p+ymnout_p+1.d-100))*(abs(dot_p)**2) + &
                              iu*((ymnout_s*ymnm_s+1.d-100)/(ymnm_s+ymnout_s+1.d-100))*(abs(dot_s)**2) 
    end do !>end-of m loop
    
    !>*******************
    !>negative order n
    !>*******************
    ginf_minus=0.d0
    do m=-n,n,1
      !>allowed diffraction wavenumber (1/m) 
      betam = kx0 + m*(2.d0*pi/l)
      betan = ky0 - n*(2.d0*pi/l)
      qmnout= sqrt(dcmplx(kout**2 - betam**2 -betan**2))
      qmnm  = sqrt(dcmplx(km**2 - betam**2 -betan**2))
    
      !>admittance
      ymnout_p = (epsout*k0+1.d-100)/(qmnout+1.d-100) 
      ymnout_s = qmnout/k0
    
      ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100) 
      ymnm_s = qmnm/k0 
     
      !>obtain <betam,betan|TE10> 
      dot=dot_bloch_te01(l,a,betam,betan)
      if ((m==0).and.(n==0)) then
        dot_p = dot
        dot_s = 0.d0
      else
        dot_p = (betam/sqrt(betam**2+betan**2))*dot
        dot_s = (betan/sqrt(betam**2+betan**2))*dot 
      end if
      
      !>update the ginf
      ginf_minus = ginf_minus + iu*((ymnout_p*ymnm_p+1.d-100)/(ymnm_p+ymnout_p+1.d-100))*(abs(dot_p)**2) + &
                                iu*((ymnout_s*ymnm_s+1.d-100)/(ymnm_s+ymnout_s+1.d-100))*(abs(dot_s)**2) 
    end do !>end-do of n loop
    
    !>***********************
    !>update total ginf_old
    !>***********************
    ginf = ginf_old + ginf_plus + ginf_minus
    
    !>*******************
    !>m=n,nn=-(n-1)->(n-1)
    !>*******************
    m=n
    do nn=-(n-1),n-1,1 
      !>allowed diffraction wavenumber (1/m) 
      betam = kx0 + m*(2.d0*pi/l)
      betan = ky0 + nn*(2.d0*pi/l)
      qmnout= sqrt(dcmplx(kout**2 - betam**2 -betan**2))
      qmnm  = sqrt(dcmplx(km**2 - betam**2 -betan**2))
    
      !>admittance
      ymnout_p = (epsout*k0+1.d-100)/(qmnout+1.d-100) 
      ymnout_s = qmnout/k0
    
      ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100) 
      ymnm_s = qmnm/k0 
     
      !>obtain <betam,betan|TE10> 
      dot=dot_bloch_te01(l,a,betam,betan)
      if ((m==0).and.(n==0)) then
        dot_p = dot
        dot_s = 0.d0
      else
        dot_p = (betam/sqrt(betam**2+betan**2))*dot
        dot_s = (betan/sqrt(betam**2+betan**2))*dot 
      end if
      
      !>update the ginf
      ginf = ginf + iu*((ymnout_p*ymnm_p+1.d-100)/(ymnm_p+ymnout_p+1.d-100))*(abs(dot_p)**2) + &
                    iu*((ymnout_s*ymnm_s+1.d-100)/(ymnm_s+ymnout_s+1.d-100))*(abs(dot_s)**2) 
    end do
      
    !>*******************
    !>m=-n,nn=-(n-1)->(n-1)
    !>*******************
    m=-n
    do nn=-(n-1),n-1,1
      !>allowed diffraction wavenumber (1/m) 
      betam = kx0 + m*(2.d0*pi/l)
      betan = ky0 + nn*(2.d0*pi/l)
      qmnout= sqrt(dcmplx(kout**2 - betam**2 -betan**2))
      qmnm  = sqrt(dcmplx(km**2 - betam**2 -betan**2))
    
      !>admittance
      ymnout_p = (epsout*k0+1.d-100)/(qmnout+1.d-100) 
      ymnout_s = qmnout/k0
    
      ymnm_p = (epsm*k0+1.d-100)/(qmnm+1.d-100) 
      ymnm_s = qmnm/k0 
     
      !>obtain <betam,betan|TE10> 
      dot=dot_bloch_te01(l,a,betam,betan)
      if ((m==0).and.(n==0)) then
        dot_p = dot
        dot_s = 0.d0
      else
        dot_p = (betam/sqrt(betam**2+betan**2))*dot
        dot_s = (betan/sqrt(betam**2+betan**2))*dot 
      end if
      
      !>update the ginf
      ginf = ginf + iu*((ymnout_p*ymnm_p+1.d-100)/(ymnm_p+ymnout_p+1.d-100))*(abs(dot_p)**2) + &
                    iu*((ymnout_s*ymnm_s+1.d-100)/(ymnm_s+ymnout_s+1.d-100))*(abs(dot_s)**2) 
    end do

    !>*************************
    !>check the convergence of ginf
    !>*************************
    ginf_err = ((abs(ginf-ginf_old)+1.d-200)/(abs(ginf_old)+1.d-200))*100.d0
      
    if (ginf_err<ginf_err_accpt) then
      gout = ginf
      maxn_out = n   
      return
    else
      ginf_old = ginf
    end if
    
  end if !>end-if of n>=1 loop
 
  !>update n
  n = n+1
end do !>end-do of n loop

end subroutine gout_0layer_sub


complex(8) function gv_func(kh,nh,muh,h,a)
implicit none
complex(8), intent(in) :: kh
complex(8), intent(in) :: nh
complex(8), intent(in) :: muh 
real(8),    intent(in) :: h 
real(8),    intent(in) :: a 

!>local variables
real(8), parameter :: pi = 3.d0*acos(0.5d0) 
complex(8) :: gv
real(8)    :: gx10
complex(8) :: gz10 

!>gx10
gx10 = pi/a 

!>gz10 
gz10 = sqrt(dcmplx(kh**2 - gx10**2))

!>gv
gv = (gz10/(kh/nh))*((1.d0+1.d-200)/(sin(gz10*h)+1.d-200))

gv_func = gv 

return 
end function gv_func

complex(8) function gam_func(kh,nh,muh,h,a)
implicit none
complex(8), intent(in) :: kh
complex(8), intent(in) :: nh
complex(8), intent(in) :: muh 
real(8),    intent(in) :: h 
real(8),    intent(in) :: a 

!>local variables
real(8),    parameter :: pi = 3.d0*acos(0.5d0) 
complex(8), parameter :: iu=dcmplx(0.d0,1.d0)
complex(8) :: gam
real(8)    :: gx10
complex(8) :: gz10 
complex(8) :: tan_x 

!>gx10
gx10 = pi/a 

!>gz10 
gz10 = sqrt(dcmplx(kh**2 - gx10**2))

!>obtain tan_x
tan_x = -iu*(exp(iu*gz10*h) - exp(-iu*gz10*h))/(exp(iu*gz10*h) + exp(-iu*gz10*h))

!>gv
gam = (gz10/(kh/nh))*((1.d0+1.d-200)/(tan_x+1.d-200))

gam_func = gam
return 
end function gam_func

complex(8) function i0_func(epsin,nin,kin,kx0,ky0,l,a,yin,ym,yd,kzd,d,kzin)
implicit none
real(8), intent(in) :: epsin
real(8), intent(in) :: nin
real(8), intent(in) :: kin
real(8), intent(in) :: kx0
real(8), intent(in) :: ky0 
real(8), intent(in) :: l
real(8), intent(in) :: a
real(8), intent(in)    :: yin 
complex(8), intent(in) :: ym
complex(8), intent(in) :: yd  
complex(8), intent(in) :: kzd
real(8), intent(in) :: d
real(8), intent(in) :: kzin

!>local variables
complex(8), parameter :: iu = (0.d0,1.d0)
real(8)    :: sincx      !sin(0.5*kx0*w)/(0.5*kx0*w)
real(8)    :: dot        !<waveguide|kx0> 
complex(8) :: rho_m_d,rho_d_in        !reflection coefficients
complex(8) :: tau_m_d,tau_in_d        !transmission coefficients

interface
 
  complex(8) function dot_bloch_te01(l,a,betam,betan)
    implicit none
    real(8), intent(in) :: l
    real(8), intent(in) :: a
    real(8), intent(in) :: betam 
    real(8), intent(in) :: betan 
  end function dot_bloch_te01

end interface

!>obtain dot product 
dot = dot_bloch_te01(l,a,kx0,ky0) 

!>reflection cofficients
rho_m_d  = (ym-yd)/(ym+yd)
rho_d_in = (yd-yin)/(yd+yin)

!>transmission coefficients
tau_m_d  = 2.d0*ym/(ym+yd)
tau_in_d = 2.d0*yin/(yin+yd)

!>obtain i0
!i0_func = iu*yd*((tau_m_d*tau_in_d*exp(iu*kzd*d))/(1.d0+rho_m_d*rho_d_in*exp(2.d0*iu*kzd*d))) &
!          *dot*exp(-iu*kzin*d)
         
i0_func = iu*yd*((tau_m_d*tau_in_d*exp(iu*kzd*d))/(1.d0+rho_m_d*rho_d_in*exp(2.d0*iu*kzd*d))) &
          *dot

return
end function i0_func

complex(8) function eps_Ag(f0)
implicit none
!>
!>dielectric constant of silver from the Drude model 
!>Reference : Optical properties of metallic films for vertical-cavity optoeletronic devices 
!>
real(8),intent(in) :: f0           !frequency (1/s)

!>local variables
real(8),parameter    :: pi = 3.d0*acos(0.5d0)
complex(8),parameter :: iu=(0.d0,1.d0)
real(8),parameter    :: hbar = 6.58211928d-16            !planck constant eV*s

!>angular frequency 
real(8) :: w0 

!>parameters of Drude's model 
real(8) :: wp,fp                !plasma frequency     
real(8) :: sig0                 !oscillator strength 
real(8) :: gam0                 !damping constant
complex(8) :: epsf              !free dielectric constant 

!>parameters of Lorentz's model
real(8) :: sigj(5)          !oscillator strength
real(8) :: wj(5)            !resonance frequency
real(8) :: gamj(5)          !damping constants 
complex(8) :: epsb          !bound dielectric constant
integer :: j

!>angular frequency in eV unit 
w0 = hbar*2.d0*pi*f0          

!>************
!>obtain epsf
!>************
wp = 9.01d0             !eV                  
sig0 = 0.845d0
gam0 = 0.048d0          !eV
epsf = 1.d0 - sig0*(wp**2)/(w0*w0 + iu*w0*gam0) 

!>************
!>obtain epsb
!>************
epsb=0.d0
do j=1,5
  if (j==1) then
    sigj(j) = 0.065d0
    gamj(j) = 3.886d0   !(eV)
    wj(j)   = 0.816d0   !(eV)
  else if (j==2) then
    sigj(j) = 0.124d0 
    gamj(j) = 0.452d0   !(eV)
    wj(j)   = 4.481d0   !(eV)
  else if (j==3) then
    sigj(j) = 0.011d0 
    gamj(j) = 0.065d0   !(eV)
    wj(j)   = 8.185d0   !(eV)
  else if (j==4) then
    sigj(j) = 0.840d0
    gamj(j) = 0.916d0   !(eV)
    wj(j)   = 9.083d0   !(eV)
  else if (j==5) then
    sigj(j) = 5.646d0 
    gamj(j) = 2.419d0   !(eV)
    wj(j)   = 20.29d0   !(eV)
  end if
  
  epsb = epsb + sigj(j)*(wp**2)/(wj(j)**2 - w0**2 - iu*w0*gamj(j))
end do

eps_Ag = epsf+epsb
return
end function eps_Ag

complex(8) function eps_Ag_ordal(f0)
implicit none
!>
!>dielectric constant of silver from the Drude model 
!>Reference : M. A. Ordal et al., 
!>            "Optical properties of the metals Al, Co, Cu, Au, Fe, Pb, Ni, Pd,
!>             Pt, Ag, Ti, and W in the infrared and far infrared," Appl. Opt. 22(7), 1099-1120 (1983). 
!>
real(8),intent(in) :: f0           !frequency (1/s)

!>local variable 
complex(8), parameter :: iu=(0.d0,1.d0) 
real(8), parameter :: c = 299792458.d0         !speed of light (m/s)
real(8), parameter :: kp=7.25d+4   !plasma wavenumber (1/cm)
real(8), parameter :: kt=1.45d+2   !relaxation constant (1/cm)\
complex(8) :: epsc                 !complex dielectric constant 
real(8) :: epsc_1,epsc_2           !real/imaginary part of epsc 
real(8) :: k0                      !wavenumber of f0 
real(8) :: lam0                    !wavelength of f0 

!>convert f0 to k0 
lam0 = c/f0            !m 
lam0 = lam0*1.d+2      !cm 
k0   = 1.d0/lam0       !1/cm 

!>obtain epsc 
epsc_1 = 1.d0 - kp**2/(k0**2 + kt**2) 
epsc_2 = (kp**2)*kt/(k0*(k0**2+kt**2))  
epsc = dcmplx(epsc_1,epsc_2)

!>return 
eps_Ag_ordal = epsc 
return 

end function eps_Ag_ordal

complex(8) function eps_Ni(f0)
implicit none
!>
!>dielectric constant of silver from the Drude model 
!>Reference : Optical properties of metallic films for vertical-cavity optoeletronic devices 
!>
real(8),intent(in) :: f0           !frequency (1/s)

!>local variables
real(8),parameter    :: pi = 3.d0*acos(0.5d0)
complex(8),parameter :: iu=(0.d0,1.d0)
real(8),parameter    :: hbar = 6.58211928d-16            !planck constant eV*s

!>angular frequency 
real(8) :: w0 

!>parameters of Drude's model 
real(8) :: wp,fp                !plasma frequency     
real(8) :: sig0                 !oscillator strength 
real(8) :: gam0                 !damping constant
complex(8) :: epsf              !free dielectric constant 

!>parameters of Lorentz's model
real(8) :: sigj(4)          !oscillator strength
real(8) :: wj(4)            !resonance frequency
real(8) :: gamj(4)          !damping constants 
complex(8) :: epsb          !bound dielectric constant
integer :: j

!>angular frequency in eV unit 
w0 = hbar*2.d0*pi*f0          

!>************
!>obtain epsf
!>************
wp = 15.92d0             !eV                  
sig0 = 0.096d0
gam0 = 0.048d0          !eV
epsf = 1.d0 - sig0*(wp**2)/(w0*w0 + iu*w0*gam0) 

!>************
!>obtain epsb
!>************
epsb=0.d0
do j=1,4
  if (j==1) then
    sigj(j) = 0.100d0
    gamj(j) = 4.511d0   !(eV)
    wj(j)   = 0.174d0   !(eV)
  else if (j==2) then
    sigj(j) = 0.135d0 
    gamj(j) = 1.334d0   !(eV)
    wj(j)   = 0.582d0   !(eV)
  else if (j==3) then
    sigj(j) = 0.106d0 
    gamj(j) = 2.178d0   !(eV)
    wj(j)   = 1.597d0   !(eV)
  else if (j==4) then
    sigj(j) = 0.729d0
    gamj(j) = 6.292d0   !(eV)
    wj(j)   = 6.089d0   !(eV)
  end if
  
  epsb = epsb + sigj(j)*(wp**2)/(wj(j)**2 - w0**2 - iu*w0*gamj(j))
end do

eps_Ni = epsf+epsb
return
end function eps_Ni
