c==============================================================================

      subroutine fcn()

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      parameter (nprmax=1000, nprdatmax=3000, npxdtmax=3000)
      parameter (nmaxdata=50000)
      integer npxdt1,npxdt2,nprdt1,nprdt2,npr,npr2
      integer pixel(2), ccd_id(2), poly
      double precision pxdt(npxdtmax),pxflx(npxdtmax),l2,sfcpm,
     & pxsg(npxdtmax),prflx(nprmax*nprdatmax),tprpr(nprmax*nprmax),
     & cpmflx(npxdtmax),sgcpmflx(npxdtmax),cpmrs(npxdtmax)
      character*1024 infile2, tmpfile3
      character*80 fitname
      double precision t0, u0, tE, chi2_cpm
      double precision yfit(nmaxdata)

      infile2 = 'lcOB160241i'
      fitname = 'test'
      t0 = 7492.77
      u0 = 0.07
      tE = 70.8

      call initcpm(npxdt1,npxdt2,nprdt1,nprdt2,npr,npr2,prflx,
     &             pxdt,pxflx,pxsg,l2,ccd_id,pixel,infile2,tprpr,
     &             nprmax,nprdatmax,npxdtmax)


      do 100 i=1,nprdt2
          yfit(i) = 2.0*i
c          write(*,*) yfit(i)
 100  continue


      call run_cpm_cpp(yfit,maxdata,sfcpm,chi2_cpm,
     &           npxdt2,nprdt2,npr,prflx,pxdt,pxflx,pxsg,l2,npr2,
     &           tprpr,nprmax,nprdatmax,npxdtmax,fitname,cpmflx,
     &           sgcpmflx,cpmrs)

       END

c===============================================================================
c This subroutine runs the CPM+LCM external library. 
c Following subroutines added by C. RANC
c===============================================================================
      subroutine run_cpm_cpp(yfit,maxdata,sfcpm,chi2_cpm,
     &           npxdt2,nprdt2,npr,prflx,pxdt,pxflx,pxsg,l2,npr2,
     &           tprpr,nprmax,nprdatmax,npxdtmax,fitname,cpmflx,
     &           sgcpmflx,cpmrs)

      integer npxdt2,nprdt2,npr,npr2,fllc,i,k,nprmax,nprdatmax,npxdtmax
      double precision pxdt(npxdtmax),pxflx(npxdtmax),l2,sfcpm,
     & pxsg(npxdtmax),prflx(nprmax*nprdatmax),tprpr(nprmax*nprmax),
     & yfit(maxdata),chi2_cpm,cpmflx(npxdtmax),
     & sgcpmflx(npxdtmax),cpmrs(npxdtmax)
      character*80 fitname 

c     Add the microlensing model
c     --------------------------
      do 100 i=1,nprdt2
          prflx(i*npr2) = yfit(i)
 100  continue

c     Run the C++ routine
c     -------------------
      fllc=0
      sfcpm=0
      call run_cpm(pxdt,pxflx,pxsg,prflx,tprpr,npxdt2,nprdt2,npr,npr2,
     &             l2,sfcpm,chi2_cpm,fllc,fitname,cpmflx,sgcpmflx,cpmrs)

      write(*,*) chi2_cpm

      end
c===============================================================================
      subroutine run_cpmlc_cpp(yfit,maxdata,mxclr,k,sfcpm,chi2_cpm,
     &           npxdt2,nprdt2,npr,prflx,pxdt,pxflx,pxsg,l2,npr2,
     &           tprpr,nprmax,nprdatmax,npxdtmax,fitname,cpmflx,
     &           sgcpmflx,cpmrs)

      integer npxdt2,nprdt2,npr,npr2,fllc,i,k,nprmax,nprdatmax,npxdtmax
      double precision pxdt(npxdtmax),pxflx(npxdtmax),l2,sfcpm,
     & pxsg(npxdtmax),prflx(nprmax*nprdatmax),tprpr(nprmax*nprmax),
     & yfit(maxdata,0:mxclr),chi2_cpm,cpmflx(npxdtmax),
     & sgcpmflx(npxdtmax),cpmrs(npxdtmax)
      character*80 fitname 

c     Add the microlensing model
c     --------------------------
      do 100 i=1,nprdt2
          prflx(i*npr2) = yfit(i,k)
 100  continue

c     Run the C++ routine
c     -------------------
      fllc=1
      sfcpm=0
      call run_cpm(pxdt,pxflx,pxsg,prflx,tprpr,npxdt2,nprdt2,npr,npr2,
     &             l2,sfcpm,chi2_cpm,fllc,fitname,cpmflx,sgcpmflx,cpmrs)

      end
c===============================================================================
      subroutine initcpm(npxdt1,npxdt2,nprdt1,nprdt2,npr,npr2,
     & prflx,pxdt,pxflx,pxsg,l2,ccd_id,pixel,infile2,tprpr,nprmax,
     & nprdatmax,npxdtmax)

      integer nprmax,nprdatmax,npxdtmax
      integer npxdt1,npxdt2,nprdt1,nprdt2,npr,npr2,poly
      integer i,j,k,pixel(2),ccd_id(2)
      double precision l2,x,y,z
      double precision pxdt(npxdtmax),pxflx(npxdtmax),
     & pxsg(npxdtmax),prflx(nprmax*nprdatmax),tprpr(nprmax*nprmax)
      logical pxmsk(npxdtmax)
      character*128 prefix 
      character*4 pixel_char(2), ccd_char(2) 
      character*1024 fname, infile2
      character tf

c     Load the paraters for cpm
c     -------------------------
c     We take an example.
      pixel(1) = 1021
      pixel(2) = 118
      ccd_id(1) = 91
      ccd_id(2) = 49
      poly = 0
      l2 = 1000.0
      write(6,*) 'CPM epoch ', ccd_id(1) 
      write(6,*) 'CPM CCD ', ccd_id(2)
      write(6,*) 'CPM Pixel ', pixel(1), pixel(2)
      write(6,*) 'Tuning params (l2, poly):', l2, poly
      
      write(ccd_char(1), '(I4)') ccd_id(1) 
      write(ccd_char(2), '(I4)') ccd_id(2) 
      write(pixel_char(1), '(I4)') pixel(1) 
      write(pixel_char(2), '(I4)') pixel(2) 
      prefix = trim(adjustl(ccd_char(1)))//'_' //
     & trim(adjustl(ccd_char(2)))//'_'//
     & trim(adjustl(pixel_char(1)))//'_'//
     & trim(adjustl(pixel_char(2)))//'_'

c     Load the timeserie of the pixel
c     -------------------------------
      npxdt1 = 0
      npxdt2 = 0
      fname = trim(adjustl(prefix))//'epoch_mask.dat'
      open(unit=1,file=fname,status='old',err=101)
      do j=1,9999
            read(1,*,end=101) tf
            if (tf .eq. 'T') then
               pxmsk(j)=.TRUE.
               npxdt2 = npxdt2 + 1
            else 
               pxmsk(j)=.FALSE.
            endif
            npxdt1 = npxdt1 + 1
      end do  
      close(1)

 101  fname = trim(adjustl(prefix))//'pixel_flux.dat'
      i = 0
      open(unit=1,file=fname,status='old',err=103)
      do j=1,npxdt1
            read(1,*,end=103) x, y, z 
            if (pxmsk(j).eqv..TRUE.) then 
               i = i + 1
               pxdt(i)=x
               pxflx(i)=y
               pxsg(i)=z
            endif
      end do  
      close(1)

c     Load the predictors flux
c     ------------------------
 103  nprdt1 = 0
      nprdt2 = 0
      fname = trim(adjustl(prefix))//'predictor_epoch_mask.dat'
      open(unit=1,file=fname,status='old',err=105)
      do j=1,9999
            read(1,*,end=105) tf
            if (tf .eq. 'T') then
               pxmsk(j)=.TRUE.
               nprdt2 = nprdt2 + 1
            else 
               pxmsk(j)=.FALSE.
            endif
            nprdt1 = nprdt1 + 1
      end do  
      close(1)

 105  npr = 0
      fname = trim(adjustl(prefix))//'pre_matrix_xy.dat'
      open(unit=1,file=fname,status='old',err=108)
      do j=1,9999
          do i=1,nprdt2
             read(1,*,end=108) x, y, z
          end do  
          npr=npr+1
      end do  
      close(1)

 108  npr2 = npr + poly + 2

      rewind(unit=1,err=111)
      do i=1,nprdt2
          do j=1,npr
             read(1,*,end=111) x, y, z
             prflx((i-1)*npr2+j) = z
          end do  
      end do  
      close(1)

c     Add the Vandermonde matrix
c     --------------------------
 111  do i=1,nprdt2
          do j=npr+1,npr2-1
              prflx((i-1)*npr2+j) = i**(j-npr-1)
          end do  
      end do  

c     Compute an auxiliary matrix 
c     ---------------------------
      do i=1,npr2-1
          do j=1,npr2-1
              x=0
              do k=1,nprdt2 
                  x=x+prflx((k-1)*npr2+i)*prflx((k-1)*npr2+j)/pxsg(k)**2
              end do
              tprpr((i-1)*(npr2-1)+j)=x
          end do  
      end do  

c     Write an input file for the K2 data
c     -----------------------------------
      fname = infile2(1:11)//'.k2c9'
      x = 3000
      y = 300
      open(unit=7,file=fname,err=114)
      do j=1,npxdt2
            write(7,200) pxdt(j)-2450000, x, y 
      end do  
 200  format(f11.6,f7.1,f6.1)
      close(7)

 114  end

