












      MODULE Ocean_Output

      IMPLICIT NONE

      include 'netcdf.inc'

      include 'comclio.h'
!      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
!      include 'ice.com'
      include 'dynami.com'
!      include 'thermo.com'
!      include 'densit.com'
!      include 'isoslope.com'
      include 'reper.com'
      include 'comunit.h'

C     ----------------------------------------------------------------------------
C     Module Interface Specification

      PRIVATE
      PUBLIC :: open
      PUBLIC :: close
      PUBLIC :: incretime
      PUBLIC :: writet
      PUBLIC :: writed
      PUBLIC :: write_2d
      PUBLIC :: write_3d
      PUBLIC :: ncid
      PUBLIC :: status
      PUBLIC :: output

C     Codes for output file types
      PUBLIC :: Monthly_Means
      PUBLIC :: Yearly_Means

C     Codes for one dimensional fields
      PUBLIC :: ptlon
      PUBLIC :: pulon
      PUBLIC :: ptlat
      PUBLIC :: pulat
      PUBLIC :: tdepth
      PUBLIC :: wdepth
      PUBLIC :: wedges
      PUBLIC :: corners
      PUBLIC :: sflat
      PUBLIC :: sfdepth
      PUBLIC :: sfedges
      PUBLIC :: basidx
      PUBLIC :: ptlonp
      PUBLIC :: pulonp
      PUBLIC :: ptlatp
      PUBLIC :: pulatp
      PUBLIC :: time

C     fields on 3 d tracer grid
      PUBLIC :: Potential_temperature
      PUBLIC :: Salinity

C     fields on 3 d momentum grid
      PUBLIC :: Zonal_velocity
      PUBLIC :: Meridional_velocity
      PUBLIC :: Vertical_velocity

C     fields on horizontal momentum grid
      PUBLIC :: Zonal_barotropic_momentum
      PUBLIC :: Meridional_barotropic_momentum

C     fields on horizontal tracer grid
      PUBLIC :: Sea_surface_height
      PUBLIC :: Surface_temperature
      PUBLIC :: Surface_salinity
      PUBLIC :: Surface_heat_flux
      PUBLIC :: Surface_freshwater_flux
      PUBLIC :: Depth_surface_mixed_layer
      PUBLIC :: Depth_convection
      PUBLIC :: G_M_slope
      PUBLIC :: Ice_thickness
      PUBLIC :: Ice_production
      PUBLIC :: Lead_fraction
      PUBLIC :: Snow_thickness
      PUBLIC :: Snow_precipitation
      PUBLIC :: Ice_temperature
      PUBLIC :: Heat_flux_ice_base
      PUBLIC :: Zonal_ice_velocity
      PUBLIC :: Meridional_ice_velocity
      PUBLIC :: Zonal_wind_stress
      PUBLIC :: Meridional_wind_stress

C     meridional streamfunctions
      PUBLIC :: meridional_overturning_streamfunction
      PUBLIC :: meridional_heat_transport
      PUBLIC :: meridional_salt_transport

      PUBLIC :: isotope_oxygen
      PUBLIC :: isotope_deuterium
      PUBLIC :: averaged_CFC11
      PUBLIC :: averaged_CFC12

C     ----------------------------------------------------------------------------
C     Type Declarations

      TYPE Info
      CHARACTER(len=64) :: unit
      CHARACTER(len=56) :: long
      CHARACTER(len=56) :: modps
      CHARACTER(len=64) :: postopo
      CHARACTER(len=8)  :: axised
      END TYPE Info

      TYPE Infovar
      CHARACTER(len=64) :: unit
      CHARACTER(len=56) :: long
      CHARACTER(len=56) :: coord
      CHARACTER(len=56) :: miss
      END TYPE Infovar

      TYPE Infodimvar
      CHARACTER(len=64) :: unit
      CHARACTER(len=56) :: long
      CHARACTER(len=56) :: ps
      CHARACTER(len=64) :: bnds
      CHARACTER(len=56) :: coord
      CHARACTER(len=56) :: miss
      END TYPE Infodimvar

C     ----------------------------------------------------------------------------
C     Module Parameters

C     Output file types
      INTEGER, PARAMETER :: Monthly_Means =  1
      INTEGER, PARAMETER :: Yearly_Means  =  2

C     Grid coordinates
      TYPE(Info), PARAMETER ::
     & I_info = Info("units","long_name","modulo","topology","axis"),
     & I_info2 = Info("units","long_name","point_spacing","positive","edges")
C     Grid coordinates
      TYPE(Infovar), PARAMETER ::
     & I_infovar = Infovar("units","long_name","coordinates","missing_value"),
     & I_infovar2 = Infovar("units","long_name","comment","missing_value")
C     Grid coordinates
      TYPE(Infodimvar), PARAMETER ::
     & I_infodimvar   = Infodimvar("units",  "long_name", "point_spacing",
     &          "bounds",                 "coordinates", "missing_value")

C     One dimensional fields
      INTEGER, PARAMETER :: ptlon   =  1
      INTEGER, PARAMETER :: pulon   =  2
      INTEGER, PARAMETER :: ptlat   =  3
      INTEGER, PARAMETER :: pulat   =  4
      INTEGER, PARAMETER :: tdepth  =  5
      INTEGER, PARAMETER :: wdepth  =  6
      INTEGER, PARAMETER :: wedges  =  7
      INTEGER, PARAMETER :: corners =  8
      INTEGER, PARAMETER :: sflat   =  9
      INTEGER, PARAMETER :: sfdepth = 10
      INTEGER, PARAMETER :: sfedges = 11
      INTEGER, PARAMETER :: basidx  = 12
      INTEGER, PARAMETER :: ptlonp  = 13
      INTEGER, PARAMETER :: pulonp  = 14
      INTEGER, PARAMETER :: ptlatp  = 15
      INTEGER, PARAMETER :: pulatp  = 16
      INTEGER, PARAMETER :: time    = 17

C     fields on 3 d tracer grid
      INTEGER, PARAMETER :: Potential_temperature                 = 1
      INTEGER, PARAMETER :: Salinity                              = 2

C     fields on 3 d momentum grid
      INTEGER, PARAMETER :: Zonal_velocity                        = 3
      INTEGER, PARAMETER :: Meridional_velocity                   = 4
      INTEGER, PARAMETER :: Vertical_velocity                     = 5

C     fields on horizontal momentum grid
      INTEGER, PARAMETER :: Zonal_barotropic_momentum             = 6
      INTEGER, PARAMETER :: Meridional_barotropic_momentum        = 7

C     fields on horizontal tracer grid
      INTEGER, PARAMETER :: Sea_surface_height                    = 8
      INTEGER, PARAMETER :: Surface_temperature                   = 9
      INTEGER, PARAMETER :: Surface_salinity                      = 10
      INTEGER, PARAMETER :: Surface_heat_flux                     = 11
      INTEGER, PARAMETER :: Surface_freshwater_flux               = 12
      INTEGER, PARAMETER :: Depth_surface_mixed_layer             = 13
      INTEGER, PARAMETER :: Depth_convection                      = 14
      INTEGER, PARAMETER :: G_M_slope                             = 15
      INTEGER, PARAMETER :: Ice_thickness                         = 16
      INTEGER, PARAMETER :: Ice_production                        = 17
      INTEGER, PARAMETER :: Lead_fraction                         = 18
      INTEGER, PARAMETER :: Snow_thickness                        = 19
      INTEGER, PARAMETER :: Snow_precipitation                    = 20
      INTEGER, PARAMETER :: Ice_temperature                       = 21
      INTEGER, PARAMETER :: Heat_flux_ice_base                    = 22
      INTEGER, PARAMETER :: Zonal_ice_velocity                    = 23
      INTEGER, PARAMETER :: Meridional_ice_velocity               = 24
      INTEGER, PARAMETER :: Zonal_wind_stress                     = 25
      INTEGER, PARAMETER :: Meridional_wind_stress                = 26

C     meridional streamfunctions
      INTEGER, PARAMETER :: meridional_overturning_streamfunction = 27
      INTEGER, PARAMETER :: meridional_heat_transport             = 28
      INTEGER, PARAMETER :: meridional_salt_transport             = 29

      INTEGER, PARAMETER :: isotope_oxygen                        = 30
      INTEGER, PARAMETER :: isotope_deuterium                     = 31
      INTEGER, PARAMETER :: averaged_CFC11                        = 32
      INTEGER, PARAMETER :: averaged_CFC12                        = 33

      INTEGER, PARAMETER :: Compute_Time_in_Months  = 1
      INTEGER, PARAMETER :: Compute_Time_in_Years   = 2

C     ----------------------------------------------------------------------------
C     Module Variables

      LOGICAL, SAVE :: initialize_module = .TRUE.
      INTEGER, SAVE :: how_to_compute_time, ncid
      INTEGER       :: status

      CONTAINS

C     ----------------------------------------------------------------------------
C     PUBLIC FUNCTIONS AND SUBROUTINES
C     ----------------------------------------------------------------------------

C     ----------------------------------------------------------------------------
C     open

      SUBROUTINE open(file,refexp)

      INTEGER,     INTENT(in) :: file
      character*6, INTENT(in) :: refexp

      CHARACTER(len=256) :: output_filename
      LOGICAL            :: existe

      SELECT CASE (file)
      CASE (Monthly_Means)
         output_filename = "outputdata/ocean/CLIO3_"//trim(refexp)//
     &                                             "_mon_avrg.nc"
         how_to_compute_time = Compute_Time_in_Months
      CASE (Yearly_Means)
         output_filename = "outputdata/ocean/CLIO3_"//trim(refexp)//
     &                                             "_ann_avrg.nc"
         how_to_compute_time = Compute_Time_in_Years
      CASE DEFAULT
         write(iuo+66,*) 'ERROR: unable to create netcdf file ',file
      END SELECT

      INQUIRE(file=output_filename,exist=existe)
      IF (existe) THEN
         initialize_module = .FALSE.
      ELSE
         initialize_module = .TRUE.
      END IF

      IF (existe) THEN
         status=nf_open(output_filename,NF_WRITE,ncid)
      ELSE
         status=nf_create(output_filename,nf_clobber,ncid)
      END IF
      if (status.ne.nf_noerr) then
        write(iuo+66,*) 'ERROR: unable to create netcdf file ',output_filename
        stop
      endif

      IF (initialize_module) CALL initializeModule

      RETURN
      END SUBROUTINE open

C     ----------------------------------------------------------------------------
C     close

      SUBROUTINE close()

      IF (initialize_module) CALL initializeModule

      status=NF_CLOSE(ncid)
      IF (status/=nf_noerr) THEN
         write(iuo+66,*) nf_strerror(status)
      END IF

      RETURN
      END SUBROUTINE close

      SUBROUTINE incretime(ja,ijour,timrec)

      INTEGER, INTENT(in) :: ja,ijour,timrec
      INTEGER             :: dimid,varid,l
      REAL*8              :: rtim

      SELECT CASE (how_to_compute_time)
      CASE (Compute_Time_in_Months)
         rtim=dfloat((ja-1)*360+ijour)/30.0
      CASE (Compute_Time_in_Years)
         rtim=dfloat(ja)
      END SELECT

      rtim=dfloat((ja-1)*360+ijour)/30.0
c      status=NF_INQ_DIMID(ncid,'time',dimid)
c      status=NF_INQ_DIMLEN(ncid,dimid,l)
      status=NF_INQ_VARID(ncid,'time',varid)
c      l=l+1
      status=nf_put_var1_double(ncid,varid,timrec,rtim)

      IF (status/=nf_noerr) THEN
         write(iuo+66,*) nf_strerror(status)
      END IF

      RETURN
      END SUBROUTINE incretime

      SUBROUTINE writet(numid,start,tab_dim,D_data)

      INTEGER, INTENT(in) :: numid
      INTEGER, INTENT(in) :: tab_dim(4)
      INTEGER, INTENT(in) :: start(4)
      INTEGER             :: varid,dimid
      LOGICAL             :: write_test

      REAL*8, DIMENSION(:,:,:), INTENT(in)    :: D_data

      write_test = .FALSE.

      SELECT CASE (how_to_compute_time)
      CASE (Compute_Time_in_Months)
         if ( newtotvaro(numid,2)==1 ) write_test = .TRUE.
      CASE (Compute_Time_in_Years)
         if ( newtotvaro(numid,3)==1 ) write_test = .TRUE.
      END SELECT

      if ( write_test ) then

      status=NF_INQ_VARID(ncid,nametotvaro(numid,2),varid)
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,tab_dim,D_data)

      IF (status/=nf_noerr) THEN
         write(iuo+66,*) nf_strerror(status)
         call ec_error(123)
      END IF

      end if

      RETURN
      END SUBROUTINE writet

      SUBROUTINE writed(numid,start,tab_dim,D_data)

      INTEGER, INTENT(in) :: numid
      INTEGER, INTENT(in) :: tab_dim(4)
      INTEGER, INTENT(in) :: start(4)
      INTEGER             :: varid,dimid
      LOGICAL             :: write_test

      REAL*8, DIMENSION(:,:), INTENT(in)    :: D_data

      write_test = .FALSE.

      SELECT CASE (how_to_compute_time)
      CASE (Compute_Time_in_Months)
         if ( newtotvaro(numid,2)==1 ) write_test = .TRUE.
      CASE (Compute_Time_in_Years)
         if ( newtotvaro(numid,3)==1 ) write_test = .TRUE.
      END SELECT

      if ( write_test ) then

      status=NF_INQ_VARID(ncid,nametotvaro(numid,2),varid)
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,tab_dim,D_data)

      IF (status/=nf_noerr) THEN
         write(iuo+66,*) nf_strerror(status)
         call ec_error(123)
      END IF

      end if

      RETURN
      END SUBROUTINE writed

      SUBROUTINE write_3d(numid,f3d,kdim,itim)
c
c no masking needs to be applied here, this is all done when
c averaging data!
c
      INTEGER, INTENT(in) :: numid,kdim,itim
      INTEGER             :: varid
      INTEGER             :: count(4),start(4),j,k
      REAL*8              :: f3d(imax,jmax,kdim)
c
      do k=1,kdim
        start(1)=1
        count(1)=imax-2
        count(2)=1
        start(3)=k
        count(3)=1
        start(4)=itim
        count(4)=1
        status=NF_INQ_VARID(ncid,nametotvaro(numid,2),varid)
        do j=1,jmax
          start(2)=j
          status=nf_put_vara_double(ncid,varid,start,count,f3d(2,j,k))
          if (status.ne.nf_noerr) then
            status=nf_close(ncid)
            write(iuo+66,*) 'ERROR3D ',status,' writing netcdf file!'
            write(iuo+66,*) ' start:',start
            write(iuo+66,*) ' count:',count
            stop
          endif
        end do
      end do
c
      RETURN
      END SUBROUTINE write_3d

      SUBROUTINE write_2d(numid,varidin,field,msk_s,msk_u,itim)

      INTEGER, INTENT(in) :: numid,varidin,itim
      logical, INTENT(in) :: msk_s,msk_u
      REAL*8,  INTENT(in) :: field(imax,jmax)

      integer i,j,start(3),count(3),varid
      real*8 xxAA,xxWW,yyAA,yyWW,spval
      integer inAA,inWW
      REAL*8 f2d(imax,jmax)

      varid=varidin
      spval=missing_valo

c fill work array
c
      do j=1,jmax
        do i=1,imax
          f2d(i,j)=field(i,j)
        end do
      end do
c
c Perform combined WW/AA-grid masking. Indices and lon/lat ranges for
c the two distinct subgrid have been taken from netdatL.f. We undefine
c gridpoints (by setting the spval flag) that separate the two partial
c grids, but we do apply the tracer or momentum grid mask here only
c upon request.
c
      do j=1,jmax
        xxAA = xaj1 + dxaj * DFLOAT(j-1)
        yyWW = ywj1 + dywj * DFLOAT(j-1)
        do i=1,imax
          yyAA = yai1 + dyai * DFLOAT(i-1)
          inAA = 1
          if ( yyAA.le.-47. .or. yyAA.ge.68. ) inAA = 0
          if ( xxAA.lt.(2.*yyAA-12.) ) inAA = 0
          if ( xxAA.gt.(284.-2*yyAA) ) inAA = 0
          if ( xxAA.gt.190. .and. xxAA.gt.(220.-yyAA) ) inAA = 0
          if ( xxAA.gt.(2.*yyAA+237.) ) inAA = 0
          xxWW = xwi1 + dxwi * DFLOAT(i-1)
          xxWW = xxWW - 20. + untour
          xxWW = mod(xxWW,untour) + 20.
          inWW = 1
          if ( yyWW.ge.70. ) inWW = 0
          if ( (yyWW.ge.39.) .and. (xxWW.lt.42.) ) inWW = 0
          if ( yyWW.ge.30. .and. yyWW.lt.39. .and. xxWW.lt.36. ) inWW=0
          if ( (yyWW.ge.20.) .and. (xxWW.gt.263.) ) inWW = 0
cPEDRO!!!!!!!          if ( (yyWW.ge.dyy) .and. (yyWW.lt.20.) .and.
          if ( (yyWW.ge.0) .and. (yyWW.lt.20.) .and.
     &            xxWW.gt.(305.-2.*yyWW) ) inWW = 0
c
          if ( yyWW.ge.30. .and. yyWW.le.70. .and. xxWW.lt.264.
     &              .and. xxWW.gt.(219.+(70.-yyWW)*45./40.) ) inWW=0
c
          if ((inWW.eq.0).and.(inAA.eq.0)) f2d(i,j) = spval
c
c upon request: tracer-/momentum-grid masking
c
          if (msk_s.and.(msks(i,j).eq.0)) f2d(i,j) = spval
          if (msk_u.and.(msku(i,j).eq.0)) f2d(i,j) = spval
        end do
      end do
c
c dump to netcdf file
c
      start(1)=1
      count(1)=imax-2
      start(2)=1
      count(2)=1
      start(3)=itim
      count(3)=0
      if (itim.gt.0) count(3)=1
      if (varid.eq.999) status=NF_INQ_VARID(ncid,nametotvaro(numid,2),varid)
        if (status.ne.nf_noerr) then
          status=nf_close(ncid)
          write(iuo+66,*) 'ERROR2D ',status,' writing netcdf file!'
          write(iuo+66,*) numid,nametotvaro(numid,2)
          stop
        endif

      do j=1,jmax
        start(2)=j
        status=nf_put_vara_double(ncid,varid,start,count,f2d(2,j))
        if (status.ne.nf_noerr) then
          status=nf_close(ncid)
          write(iuo+66,*) 'ERROR2D ',status,' writing netcdf file!'
          write(iuo+66,*) ' start:',start
          write(iuo+66,*) ' count:',count
          stop
        endif
      end do
c
      RETURN
      END SUBROUTINE write_2d

C     ----------------------------------------------------------------------------
C     output
C
C     This function checks whether the output flags of a field indicate that
C     the field should be written to any of the output files.

      FUNCTION output(flag) RESULT(result)

      REAL*8, INTENT(in) :: flag
      LOGICAL            :: result

      result = ( flag /= 0 )

      RETURN
      END FUNCTION output


C     ----------------------------------------------------------------------------
C     PRIVATE FUNCTIONS AND SUBROUTINES
C     ----------------------------------------------------------------------------

C     ----------------------------------------------------------------------------
C     initialize

      SUBROUTINE initializeModule

      REAL*8, PARAMETER :: pi               = 2*acos(0.)
      REAL*8, PARAMETER :: radian_to_degree = 180.0/pi

      REAL*8, ALLOCATABLE, DIMENSION(:)   :: templon, templat
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: temp

      integer tab_dimid(5)
      integer tab_dim(5)

      INTEGER :: i,j,l,dimid,varid
      LOGICAL :: write_test = .FALSE.

      character(len=120):: longdata
      integer,dimension(8) :: cvalues
      REAL*8 :: realmonth

      REAL*8 :: tabdimvar(22,imax+1,jmax+1,4)

      CHARACTER*60, ALLOCATABLE, DIMENSION(:,:) :: namedimvar

      allocate(namedimvar(9,22))
      namedimvar = reshape( (/ character*60 ::  ! GFORTRAN
!      namedimvar = reshape( (/
     &    'tlon','ptlon','ptlat','','degrees_east','longitude','uneven','tlon_bounds','',
     &    'tlon_bounds','ptlon','ptlat','corners','degrees_east','longitude edges','uneven','','',
     &    'tlonp','ptlonp','ptlatp','','degrees_east','longitude bounds','uneven','','',
     &    'tlat','ptlon','ptlat','','degrees_north','latitude','uneven','tlat_bounds','',
     &    'tlat_bounds','ptlon','ptlat','corners','degrees_north','latitude edges','uneven','','',
     &    'tlatp','ptlonp','ptlatp','','degrees_north','latitude bounds','uneven','','',
     &    'ulon','pulon','pulat','','degrees_east','longitude','uneven','ulon_bounds','',
     &    'ulon_bounds','pulon','pulat','corners','degrees_east','longitude edges','uneven','','',
     &    'ulonp','pulonp','pulatp','','degrees_east','longitude bounds','uneven','','',
     &    'ulat','pulon','pulat','','degrees_north','latitude','uneven','ulat_bounds','',
     &    'ulat_bounds','pulon','pulat','corners','degrees_north','latitude edges','uneven','','',
     &    'ulatp','pulonp','pulatp','','degrees_north','latitude bounds','uneven','','',
     &    'angle','pulon','pulat','','radians','clockwise rotation of x-axis','uneven','','ulon ulat',
     &    'dxs1','ptlon','ptlat','','meters','zonal length of grid box sides','uneven','','tlon tlat',
     &    'dxs2','ptlon','ptlat','','meters','meridional length of grid box sides','uneven','','tlon tlat',
     &    'dxc1','ptlon','ptlat','','meters','zonal width at grid box centres','uneven','','tlon tlat',
     &    'dxc2','ptlon','ptlat','','meters','meridional width at grid box centres','uneven','','tlon tlat',
     &    'area','ptlon','ptlat','','meters**2','surface area of tracer grid box','uneven','','tlon tlat',
     &    'tmask','ptlon','ptlat','','none','horizontal tracer grid mask','','','tlon tlat',
     &    'umask','pulon','pulat','','none','horizontal momentum grid mask','','','ulon ulat',
     &    'h','ptlon','ptlat','','meters','bathymetry on tracer gridpoints','','','tlon tlat',
     &    'fcor','pulon','pulat','','1/seconds','Coriolis parameter on momentum grid','','','ulon ulat' /),
     & shape(namedimvar) )

      IF (initialize_module) THEN

C     -----------------------------------------------------------------------

         CALL initdim()

!         status=nf_sync(ncid)
!         IF (status/=nf_noerr) THEN
!           write(iuo+66,*) nf_strerror(status)
!         END IF

         tabdimvar( :,:,:,:)           = 0.0
         tabdimvar( 1,1,1:jmax,1)      = xslon(2,:)
         tabdimvar( 2,1,1:jmax,:)      = xsedg(2,:,:)
         tabdimvar( 3,1,:,1)           = xslonp(2,:)
         tabdimvar( 4,1,1:jmax,1)      = yslat(2,:)
         tabdimvar( 5,1,1:jmax,:)      = ysedg(2,:,:)
         tabdimvar( 6,1,:,1)           = yslatp(2,:)
         tabdimvar( 7,1,1:jmax,1)      = xulon(2,:)
         tabdimvar( 8,1,1:jmax,:)      = xuedg(2,:,:)
         tabdimvar( 9,1,:,1)           = xulonp(2,:)
         tabdimvar(10,1,1:jmax,1)      = yulat(2,:)
         tabdimvar(11,1,1:jmax,:)      = yuedg(2,:,:)
         tabdimvar(12,1,:,1)           = yulatp(2,:)
         tabdimvar(13,1,1:jmax,1)      = angle(2,:)
         tabdimvar(14,1:imax,1:jmax,1) = dxs1(:,:)
         tabdimvar(15,1:imax,1:jmax,1) = dxs2(:,:)
         tabdimvar(16,1:imax,1:jmax,1) = dxc1(:,:)
         tabdimvar(17,1:imax,1:jmax,1) = dxc2(:,:)
         tabdimvar(18,1:imax,1:jmax,1) = area(:,:)
         tabdimvar(19,1,1:jmax,1)      = real(msks(2,:))
         tabdimvar(20,1,1:jmax,1)      = real(msku(2,:))
         tabdimvar(21,1:imax,1:jmax,1) = hs(:,:)
         do j=1,jmax
           do i=1,imax
             tabdimvar(22,i,j,1)=2.0*fs2cor(i,j)
           end do
         end do

         do i=1,22
            CALL decdimvar(namedimvar(:,i))
         enddo
         status=NF_ENDDEF(ncid)
         do i=1,22
            CALL initdimvar(namedimvar(:,i),tabdimvar(i,:,:,:))
         enddo

C     -----------------------------------------------------------------------

         do i=1,numtotvaro

            SELECT CASE (how_to_compute_time)
            CASE (Compute_Time_in_Months)
               if ( newtotvaro(i,2)==1 ) write_test = .TRUE.
            CASE (Compute_Time_in_Years)
               if ( newtotvaro(i,3)==1 ) write_test = .TRUE.
            END SELECT

            if ( write_test ) then
               status=NF_REDEF(ncid)
               CALL initvar(nametotvaro(i,:),dimtotvaro(i,:))
               status=NF_INQ_VARID(ncid,nametotvaro(i,2),varid)
               status=NF_ENDDEF(ncid)
            end if

            write_test = .FALSE.

         enddo

C     -----------------------------------------------------------------------
C     Write the global attributes

         status=NF_REDEF(ncid)
c         write(*,*) globalatto(:,1)
         globalatto(6,1)="table_id"
         globalatto(6,2)="Table Omon (7 April 2004)"
         globalatto(22,1)="tracking_id"
         call system('uuidgen > uuidgen_file')
         open(84,file = 'uuidgen_file',status='old',form='formatted')
         read(84,*) longdata
         close(84)
         call system('rm -rf uuidgen_file')
         globalatto(22,2)=trim(longdata)
         globalatto(24,1)="frequency"
         SELECT CASE (how_to_compute_time)
             CASE (Compute_Time_in_Months)
                globalatto(24,2)="mon"
             CASE (Compute_Time_in_Years)
                globalatto(24,2)="yr"
         END SELECT
         call date_and_time(VALUES=cvalues)
         longdata=''
         write(longdata,999) cvalues(1), cvalues(2), cvalues(3), cvalues(5), cvalues(6), cvalues(7)
 999     format(i4.4,'-',i2.2,'-',i2.2,'T',i2.2,':',i2.2,':',i2.2,'Z')
         globalatto(25,1)="creation_date"
         globalatto(25,2)=trim(longdata)
         globalatto(26,1)="modeling_realm"
         globalatto(26,2)="ocean"
         do i=1,26
           status=NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,globalatto(i,1),len_trim(globalatto(i,2)),trim(globalatto(i,2)))
         enddo
         status=NF_ENDDEF(ncid)

C     -----------------------------------------------------------------------
C     6. Indicate that initialization has been done.

         initialize_module = .FALSE.

      END IF
      RETURN

      CONTAINS


      SUBROUTINE initdim()

      INTEGER                 :: dimid,varid
c      CHARACTER(len=256)      :: string,tmpstri,tmpstr
c      INTEGER                 :: realmonth
c
c      character*(*) avgname
c
c      integer ldim

c      parameter (ldim=26)

c-----------------------------------------------------------------------
c
      integer i,j,k,recdim,idim,jdim,kdim,idimp1,jdimp1,kdimp1,kdimp2
      integer DimIDs(17), dims(4), start(3), count(3)
c
      real*8 rdelta,rmisc(max(imax,jmax,kmax))
c
c-----------------------------------------------------------------------

      idim=imax-2
      idimp1=idim+1
      jdim=jmax
      jdimp1=jdim+1
      kdim=kmax
      kdimp1=kdim+1
      kdimp2=kdim+2
c
c-----------------------------------------------------------------------
c Define the dimensions of staggered fields. Please note that the data
c are stored in the netcdf dataset on their native gridis, unlike those
c generated by the postprocessing routine cresuint, which performs
c horizontal interpolation onto a regular grid with a spacing of 2.5
c degrees longitude. CLIO is discretized on the Arakawa B-grid, and thus
c we deal with both tracer and momentum point longitudes and latitudes
c in the horizontal, and tracer (and horizontal momentum) depths at the
c centers of the cells and vertical momentum at their bottom in the
c vertical. This adds up to 6 dimensions, plus one for the number of
c tracers, and, finally, one for time. CLIO uses a grid combined of a
c regular longitude latitude grid all over the globe except the North
c Atlantic and Arctic Ocean, where a rotated subgrid is used. The two
c meshes are being merged in the equatorial Atlantic. We thus define
c PSEUDO LONgitudes and LATitudes here...
c-----------------------------------------------------------------------
c
      status=nf_def_dim(ncid,'ptlon',idim,dimid)
      start(1)=dimid
c   pseudo tracer grid longitudes:
      status=nf_def_var(ncid,'ptlon',NF_DOUBLE,1,start,varid)
      status=nf_put_att_text(ncid,varid,'units',12,'degrees_east')
      status=nf_put_att_text(ncid,varid,'long_name',28,
     &                       'pseudo tracer grid longitude')
      status=nf_put_att_text(ncid,varid,'modulo',1,'m')
      status=nf_put_att_text(ncid,varid,'topology',8,'circular')
      status=NF_ENDDEF(ncid)
C     Output data for Variable
      start(1) = 1
      count(1) = idim
      rmisc(1)=xslon(2,1)
      rdelta=xslon(3,1)-xslon(2,1)
      do i=2,idim
        rmisc(i)=rmisc(i-1)+rdelta
      end do
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,count,rmisc)
      status=NF_REDEF(ncid)


      status=nf_def_dim(ncid,'pulon' ,idim,dimid)
      start(1)=dimid
c   pseudo momentum grid longitudes:
      status=nf_def_var(ncid,'pulon',NF_DOUBLE,1,start,varid)
      status=nf_put_att_text(ncid,varid,'units',12,'degrees_east')
      status=nf_put_att_text(ncid,varid,'long_name',30,
     &                       'pseudo momentum grid longitude')
      status=nf_put_att_text(ncid,varid,'modulo',1,'m')
      status=nf_put_att_text(ncid,varid,'topology',8,'circular')
C     Output data for Variable
      status=NF_ENDDEF(ncid)
      start(1) = 1
      count(1) = idim
      rmisc(1)=xulon(2,1)
      rdelta=xulon(3,1)-xulon(2,1)
      do i=2,idim
        rmisc(i)=rmisc(i-1)+rdelta
      end do
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,count,rmisc)
      status=NF_REDEF(ncid)


      status=nf_def_dim(ncid,'ptlat' ,jdim,dimid)
      start(1)=dimid
c   pseudo tracer grid latitudes:
      status=nf_def_var(ncid,'ptlat',NF_DOUBLE,1,start,varid)
      status=nf_put_att_text(ncid,varid,'units',13,'degrees_north')
      status=nf_put_att_text(ncid,varid,'long_name',27,
     &                       'pseudo tracer grid latitude')
C     Output data for Variable
      status=NF_ENDDEF(ncid)
      start(1) = 1
      count(1) = jdim
      rmisc(1)=yslat(1,1)
      rdelta=yslat(1,2)-yslat(1,1)
      do j=2,jdim
        rmisc(j)=rmisc(j-1)+rdelta
      end do
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,count,rmisc)
      status=NF_REDEF(ncid)


      status=nf_def_dim(ncid,'pulat' ,jdim,dimid)
      start(1)=dimid
c   pseudo momentum grid latitudes:
      status=nf_def_var(ncid,'pulat',NF_DOUBLE,1,start,varid)
      status=nf_put_att_text(ncid,varid,'units',13,'degrees_north')
      status=nf_put_att_text(ncid,varid,'long_name',29,
     &                       'pseudo momentum grid latitude')
C     Output data for Variable
      status=NF_ENDDEF(ncid)
      start(1) = 1
      count(1) = jdim
      rmisc(1)=yulat(1,1)
      rdelta=yulat(1,2)-yulat(1,1)
      do j=2,jdim
        rmisc(j)=rmisc(j-1)+rdelta
      end do
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,count,rmisc)
      status=NF_REDEF(ncid)


      status=nf_def_dim(ncid,'tdepth',kdim,dimid)
      start(1)=dimid
c   tracer (and horizontal momentum) grid depths:
      status=nf_def_var(ncid,'tdepth',NF_DOUBLE,1,start,varid)
      status=nf_put_att_text(ncid,varid,'units',6,'meters')
      status=nf_put_att_text(ncid,varid,'long_name',22,
     &                       'depth of tracer points')
      status=nf_put_att_text(ncid,varid,'point_spacing',6,'uneven')
      status=nf_put_att_text(ncid,varid,'positive',2,'up')
      status=nf_put_att_text(ncid,varid,'edges',6,'wdepth')
C     Output data for Variable
      status=NF_ENDDEF(ncid)
      status=nf_put_var_double(ncid,varid,z)
      status=NF_REDEF(ncid)


      status=nf_def_dim(ncid,'wdepth',kdimp1,dimid)
      start(1)=dimid
c   vertical momentum grid depths (also edges of tracer grid depth:
      status=nf_def_var(ncid,'wdepth',NF_DOUBLE,1,start,varid)
      status=nf_put_att_text(ncid,varid,'units',6,'meters')
      status=nf_put_att_text(ncid,varid,'long_name',33,
     &                       'depth of vertical momentum points')
      status=nf_put_att_text(ncid,varid,'point_spacing',6,'uneven')
      status=nf_put_att_text(ncid,varid,'positive',2,'up')
      status=nf_put_att_text(ncid,varid,'edges',6,'wedges')
C     Output data for Variable
      status=NF_ENDDEF(ncid)
      status=nf_put_var_double(ncid,varid,zw)
      status=NF_REDEF(ncid)


      status=nf_def_dim(ncid,'wedges',kdimp2,dimid)
      start(1)=dimid
c   edges of vertical momentum grid depths:
      status=nf_def_var(ncid,'wedges',NF_DOUBLE,1,start,varid)
      status=nf_put_att_text(ncid,varid,'units',6,'meters')
      status=nf_put_att_text(ncid,varid,'long_name',32,
     &                       'edges of vertical momentum boxes')
      status=nf_put_att_text(ncid,varid,'point_spacing',6,'uneven')
      status=nf_put_att_text(ncid,varid,'positive',2,'up')
C     Output data for Variable
      status=NF_ENDDEF(ncid)
      start(1) = 1
      count(1) = kdimp2
      rmisc(1)=zw(1)
      do k=2,kdimp1
        rmisc(k)=z(k-1)
      end do
      rmisc(kdimp2)=zw(kdimp1)
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,count,rmisc)
      status=NF_REDEF(ncid)


      status=nf_def_dim(ncid,'corners',4,dimid)


      status=nf_def_dim(ncid,'sflat',57,dimid)
      start(1)=dimid
c   latitudes at which overturning streamfunctions are evaluated:
      status=nf_def_var(ncid,'sflat',NF_DOUBLE,1,start,varid)
      status=nf_put_att_text(ncid,varid,'units',13,'degrees_north')
      status=nf_put_att_text(ncid,varid,'long_name',8,'latitude')
C     Output data for Variable
      status=NF_ENDDEF(ncid)
      start(1) = 1
      count(1) = 57
      rmisc(1)=yulat(1,1)
      rdelta=yulat(1,2)-yulat(1,1)
      do j=2,jdim
        rmisc(j)=rmisc(j-1)+rdelta
      end do
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,count,rmisc)
      status=NF_REDEF(ncid)


      status=nf_def_dim(ncid,'sfdepth',kdimp2,dimid)
      start(1)=dimid
c   depths for streamfunction computation:
      status=nf_def_var(ncid,'sfdepth',NF_DOUBLE,1,start,varid)
      status=nf_put_att_text(ncid,varid,'units',6,'meters')
      status=nf_put_att_text(ncid,varid,'long_name',20,
     &                       'streamfunction depth')
      status=nf_put_att_text(ncid,varid,'point_spacing',6,'uneven')
      status=nf_put_att_text(ncid,varid,'positive',2,'up')
      status=nf_put_att_text(ncid,varid,'edges',7,'sfedges')
C     Output data for Variable
      status=NF_ENDDEF(ncid)
      start(1) = 1
      count(1) = kdimp2
      rmisc(1)=zw(1)
      do k=1,kmax
        rmisc(k+1)=z(k)
      end do
      rmisc(kmax+2)=0.0
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,count,rmisc)
      status=NF_REDEF(ncid)


      status=nf_def_dim(ncid,'sfedges',kdimp2+1,dimid)
      start(1)=dimid
c   corresponding edges
      status=nf_def_var(ncid,'sfedges',NF_DOUBLE,1,start,varid)
      status=nf_put_att_text(ncid,varid,'units',6,'meters')
      status=nf_put_att_text(ncid,varid,'long_name',30,
     &                       'edges of streamfunction depths')
      status=nf_put_att_text(ncid,varid,'point_spacing',6,'uneven')
      status=nf_put_att_text(ncid,varid,'positive',2,'up')
C     Output data for Variable
      status=NF_ENDDEF(ncid)
      start(1) = 1
      count(1) = kdimp2+1
      rmisc(1)=zw(1)
      rmisc(2)=zw(1)+1.0d-6
      do k=2,kmax
        rmisc(k+1)=zw(k)
      end do
      rmisc(kmax+2)=-1.0d-6
      rmisc(kmax+3)=0.0
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,count,rmisc)
      status=NF_REDEF(ncid)


      status=nf_def_dim(ncid,'basidx',4,dimid)
      start(1)=dimid
c   basin index axis
      status=nf_def_var(ncid,'basidx',NF_INT,1,start,varid)
      status=nf_put_att_text(ncid,varid,'long_name',9,'longitude')
C     Output data for Variable
      status=NF_ENDDEF(ncid)
      do k=1,4
        status=nf_put_var1_int(ncid,varid,k,k)
      end do
      status=NF_REDEF(ncid)


      status=nf_def_dim(ncid,'ptlonp',idimp1,dimid)
      status=nf_def_dim(ncid,'pulonp',idimp1,dimid)
      status=nf_def_dim(ncid,'ptlatp',jdimp1,dimid)
      status=nf_def_dim(ncid,'pulatp',jdimp1,dimid)


      status=nf_def_dim(ncid,'time',nf_unlimited,dimid)
      start(1)=dimid
c   time
      status=nf_def_var(ncid,'time',NF_DOUBLE,1,start,varid)
      status=nf_put_att_text(ncid,varid,'calendar',7,'360_day')

      SELECT CASE (how_to_compute_time)
      CASE (Compute_Time_in_Months)
         status=nf_put_att_text(ncid,varid,'units',32,
     &                         'months since 0000-00-00 00:00:00')
         status=nf_put_att_text(ncid,varid,'delta_t',19,
     &                            '0000-01-00 00:00:00')
         status=nf_put_att_text(ncid,varid,'avg_period',19,
     &                            '0000-01-00 00:00:00')
      CASE (Compute_Time_in_Years)
         status=nf_put_att_text(ncid,varid,'units',31,
     &                         'years since 0000-00-00 00:00:00')
         status=nf_put_att_text(ncid,varid,'delta_t',19,
     &                            '0001-00-00 00:00:00')
         status=nf_put_att_text(ncid,varid,'avg_period',19,
     &                            '0001-00-00 00:00:00')
      END SELECT

!      IF (I_num%std=="time") THEN
!         write(tmpstr,*) realmonth
!         write(tmpstri,*) irunlabel
!         SELECT CASE (I_num%unit)
!         CASE ("years")
!           string="years since 1-"//trim(adjustl(tmpstr))//"-"//trim(adjustl(tmpstri))
!         CASE ("months")
!           string="months since 1-"//trim(adjustl(tmpstr))//"-"//trim(adjustl(tmpstri))
!         END SELECT
!      END IF
C     Output data for Variable
      status=NF_ENDDEF(ncid)
      start(1) = 1
      count(1) = 1
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,count,(/1.0d0/))
      status=NF_REDEF(ncid)

      RETURN
      END SUBROUTINE initdim



      SUBROUTINE decdimvar(tab_c)

      character*60, INTENT(in)  :: tab_c(9)
      INTEGER                   :: varid,N_num
      INTEGER                   :: tab_dimid(3),tmpint(1)
      REAL*8                    :: tmp(1)

      tab_dimid(:)=0
      N_num=2
      status=NF_INQ_DIMID(ncid,tab_c(2),dimid)
      tab_dimid(1)=dimid
      status=NF_INQ_DIMID(ncid,tab_c(3),dimid)
      tab_dimid(2)=dimid
      IF ( len_trim(tab_c(4))/=0 ) THEN
         status=NF_INQ_DIMID(ncid,tab_c(4),dimid)
         tab_dimid(3)=dimid
         N_num=3
      END IF

C     Definition
      IF ( tab_c(1).eq."tmask".or.tab_c(1).eq."umask" ) THEN
         tmpint(1)=99
         status=NF_DEF_VAR(ncid,tab_c(1),NF_INT,N_num,tab_dimid,varid)
         status=NF_PUT_ATT_INT(ncid,varid,
     &      I_infodimvar%miss,NF_INT,1,tmpint)
      ELSE
         tmp(1)=missing_valo
         status=NF_DEF_VAR(ncid,tab_c(1),NF_DOUBLE,N_num,tab_dimid,varid)
         status=NF_PUT_ATT_DOUBLE(ncid,varid,
     &      I_infodimvar%miss,NF_DOUBLE,1,tmp)
      END IF

C     Define the attributes
      status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_infodimvar%unit,len_trim(tab_c(5)),trim(tab_c(5)))
      status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_infodimvar%long,len_trim(tab_c(6)),trim(tab_c(6)))
      IF ( len_trim(tab_c(7))/=0 ) status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_infodimvar%ps,len_trim(tab_c(7)),trim(tab_c(7)))

      IF ( len_trim(tab_c(8))/=0 ) status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_infodimvar%bnds,len_trim(tab_c(8)),trim(tab_c(8)))
      IF ( len_trim(tab_c(9))/=0 ) status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_infodimvar%coord,len_trim(tab_c(9)),trim(tab_c(9)))

      RETURN
      END SUBROUTINE decdimvar


      SUBROUTINE initdimvar(tab_c,tabdimvar)

      character*60, INTENT(in)  :: tab_c(9)
      REAL*8, INTENT(in)        :: tabdimvar(imax+1,jmax+1,4)
      REAL*8                    :: tabdimvarint(imax+1,jmax+1)
      INTEGER                   :: j,n,varid1,varid2,varid3,varid4,varid5
      INTEGER                   :: varid,start(3),count(3)

      status=NF_INQ_VARID(ncid,tab_c(1),varid)

      start(1)=1
      count(1)=imax-2
      start(2)=1
      count(2)=1
      start(3)=0
      count(3)=0
      SELECT CASE (tab_c(1))
         CASE ("tlon")
            status=NF_INQ_VARID(ncid,"tlon",varid1)
            status=NF_INQ_VARID(ncid,"tlat",varid2)
            status=NF_INQ_VARID(ncid,"ulon",varid3)
            status=NF_INQ_VARID(ncid,"ulat",varid4)
            status=NF_INQ_VARID(ncid,"angle",varid5)
            do j=1,jmax
              start(2)=j
              status=nf_put_vara_double(ncid,varid1,start,count,xslon(2,j))
              status=nf_put_vara_double(ncid,varid2,start,count,yslat(2,j))
              status=nf_put_vara_double(ncid,varid3,start,count,xulon(2,j))
              status=nf_put_vara_double(ncid,varid4,start,count,yulat(2,j))
              status=nf_put_vara_double(ncid,varid5,start,count,angle(2,j))
            end do
         CASE ("tlon_bounds")
            status=NF_INQ_VARID(ncid,"tlon_bounds",varid1)
            status=NF_INQ_VARID(ncid,"tlat_bounds",varid2)
            status=NF_INQ_VARID(ncid,"ulon_bounds",varid3)
            status=NF_INQ_VARID(ncid,"ulat_bounds",varid4)
            do n=1,4
              start(3)=n
              count(3)=1
              do j=1,jmax
                start(2)=j
                status=nf_put_vara_double(ncid,varid1,start,count,xsedg(2,j,n))
                status=nf_put_vara_double(ncid,varid2,start,count,ysedg(2,j,n))
                status=nf_put_vara_double(ncid,varid3,start,count,xuedg(2,j,n))
                status=nf_put_vara_double(ncid,varid4,start,count,yuedg(2,j,n))
              end do
            end do
         CASE ("tlonp")
            status=NF_INQ_VARID(ncid,"tlonp",varid1)
            status=NF_INQ_VARID(ncid,"tlatp",varid2)
            status=NF_INQ_VARID(ncid,"ulonp",varid3)
            status=NF_INQ_VARID(ncid,"ulatp",varid4)
            count(1)=imax-1
            do j=1,jmax+1
              start(2)=j
              status=nf_put_vara_double(ncid,varid1,start,count,xslonp(2,j))
              status=nf_put_vara_double(ncid,varid2,start,count,yslatp(2,j))
              status=nf_put_vara_double(ncid,varid3,start,count,xulonp(2,j))
              status=nf_put_vara_double(ncid,varid4,start,count,yulatp(2,j))
            end do
         CASE ("tmask")
            status=NF_INQ_VARID(ncid,"tmask",varid1)
            status=NF_INQ_VARID(ncid,"umask",varid2)
            do j=1,jmax
              start(2)=j
              status=nf_put_vara_int(ncid,varid1,start,count,msks(2,j))
              status=nf_put_vara_int(ncid,varid2,start,count,msku(2,j))
            end do
         CASE ("tlat","ulon","ulat","angle","tlatp","ulonp","ulatp",
     &                "tlat_bounds","ulon_bounds","ulat_bounds","umask")
         CASE DEFAULT
            call write_2d(0,varid,tabdimvar(:,:,1),.true.,.false.,0)
      END SELECT

      IF (status/=nf_noerr) THEN
         write(iuo+66,*) nf_strerror(status)
         !call ec_error(123)
      END IF

      RETURN
      END SUBROUTINE initdimvar


      SUBROUTINE initvar(tab_c,tab_d)

      character*80, INTENT(in)  :: tab_c(4)
      character*60, INTENT(in)  :: tab_d(4)
      character*80              :: tempstr
      INTEGER                   :: j,varid,N_num
      INTEGER                   :: tab_dimid(4)
      REAL*8                    :: tmp(1)

      N_num=3
      IF (tab_d(4)/="_") N_num=4

      do j=1,4
        tab_dimid(j)=convertDim(tab_d(j))
      end do
C     Definition
      status=NF_DEF_VAR(ncid,tab_c(2),NF_DOUBLE,N_num,tab_dimid,varid)

      tempstr=tab_c(4)
      j=index(tempstr,"_")
      do while(j.ne.0)
        tempstr(j:j)=" "
        j=index(tempstr,"_")
      enddo

C     Define the attributes
      status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_infovar%unit,len_trim(tab_c(3)),trim(tab_c(3)))
      status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_infovar%long,len_trim(tab_c(1)),trim(tab_c(1)))
      IF ( tab_d(1)=="basidx" ) THEN
         status=NF_PUT_ATT_TEXT(ncid,varid,
     &         I_infovar2%coord,len_trim(tempstr),trim(tempstr))
      ELSE
         status=NF_PUT_ATT_TEXT(ncid,varid,
     &         I_infovar%coord,len_trim(tempstr),trim(tempstr))
      END IF

      tmp(1)=missing_valo
      status=NF_PUT_ATT_DOUBLE(ncid,varid,
     &      I_infovar%miss,NF_DOUBLE,1,tmp)

      IF (status/=nf_noerr) THEN
         write(iuo+66,*) nf_strerror(status)
         call ec_error(123)
      END IF

      RETURN
      END SUBROUTINE initvar

      END SUBROUTINE initializeModule



      FUNCTION convertDim(striid) RESULT(result)

      character*60, INTENT(in) :: striid
      INTEGER                  :: result

      SELECT CASE (striid)
      CASE ("ptlon")
        result=ptlon
      CASE ("pulon")
        result=pulon
      CASE ("ptlat")
        result=ptlat
      CASE ("pulat")
        result=pulat
      CASE ("tdepth")
        result=tdepth
      CASE ("wdepth")
        result=wdepth
      CASE ("wedges")
        result=wedges
      CASE ("corners")
        result=corners
      CASE ("sflat")
        result=sflat
      CASE ("sfdepth")
        result=sfdepth
      CASE ("sfedges")
        result=sfedges
      CASE ("basidx")
        result=basidx
      CASE ("ptlonp")
        result=ptlonp
      CASE ("pulonp")
        result=pulonp
      CASE ("ptlatp")
        result=ptlatp
      CASE ("pulatp")
        result=pulatp
      CASE ("time")
        result=time
      END SELECT

      RETURN
      END FUNCTION convertDim



      END MODULE Ocean_Output
