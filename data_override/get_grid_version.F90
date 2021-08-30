!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!> @defgroup get_grid_version_fms2_io_mod get_grid_version_fms2_io_mod
!> @ingroup data_override
!> @brief fms2_io implementations of grid routines for @ref data_override_mod

!> @file
!> @brief File for @ref get_grid_version_mod

!> @addtogroup get_grid_version_mod
!> @{
module get_grid_version_mod
use constants_mod, only: PI
use mpp_mod, only : mpp_error,FATAL,NOTE, mpp_min, mpp_max
use mpp_domains_mod, only : domain2d, operator(.NE.),operator(.EQ.)
use mpp_domains_mod, only : mpp_get_global_domain, mpp_get_data_domain
use fms2_io_mod,     only : FmsNetcdfDomainFile_t, FmsNetcdfFile_t, open_file, close_file, &
                            variable_exists, read_data, get_variable_size, get_variable_num_dimensions
use mosaic2_mod,     only : get_mosaic_tile_grid
use platform_mod

implicit none

real(r4_kind), parameter    :: deg_to_radian_4=PI/180.
real(r8_kind), parameter    :: deg_to_radian_8=PI/180.
contains

!> Get lon and lat of three model (target) grids from grid_spec.nc
subroutine check_grid_sizes(domain_name, Domain, nlon, nlat)
character(len=12), intent(in) :: domain_name
type (domain2d),   intent(in) :: Domain
integer,           intent(in) :: nlon, nlat

character(len=184) :: error_message
integer            :: xsize, ysize

call mpp_get_global_domain(Domain, xsize=xsize, ysize=ysize)
if(nlon .NE. xsize .OR. nlat .NE. ysize) then
  error_message = 'Error in data_override_init. Size of grid as specified by '// &
                  '             does not conform to that specified by grid_spec.nc.'// &
                  '  From             :     by      From grid_spec.nc:     by    '
  error_message( 59: 70) = domain_name
  error_message(130:141) = domain_name
  write(error_message(143:146),'(i4)') xsize
  write(error_message(150:153),'(i4)') ysize
  write(error_message(174:177),'(i4)') nlon
  write(error_message(181:184),'(i4)') nlat
  call mpp_error(FATAL,error_message)
endif
end subroutine check_grid_sizes

!> Get global lon and lat of three model (target) grids, with a given file name
subroutine get_grid_version_1(grid_file, mod_name, domain, isc, iec, jsc, jec, lon, lat, min_lon, max_lon, grid_center_bug)
  character(len=*),            intent(in) :: grid_file !< name of grid file
  character(len=*),            intent(in) :: mod_name !< module name
  type(domain2d),              intent(in) :: domain !< 2D domain
  integer,                     intent(in) :: isc, iec, jsc, jec
  class(*), dimension(isc:,jsc:), intent(out) :: lon, lat
  class(*),                       intent(out) :: min_lon, max_lon
  logical,           intent(in), optional :: grid_center_bug !< Enables legacy behaviour

  integer                                      :: i, j, siz(4)
  integer                                      :: nlon, nlat !< size of global lon and lat
  real(r4_kind), dimension(:,:,:), allocatable :: lon_vert_4, lat_vert_4 !< of OCN grid vertices
  real(r8_kind), dimension(:,:,:), allocatable :: lon_vert_8, lat_vert_8 !< of OCN grid vertices
  real(r4_kind), dimension(:),     allocatable :: glon_4, glat_4  !< lon and lat of 1-D grid of atm/lnd
  real(r8_kind), dimension(:),     allocatable :: glon_8, glat_8  !< lon and lat of 1-D grid of atm/lnd
  logical                                      :: is_new_grid
  integer                                      :: is, ie, js, je
  integer                                      :: isd, ied, jsd, jed
  integer                                      :: isg, ieg, jsg, jeg
  character(len=3)                             :: xname, yname
  integer                                      :: start(2), nread(2)
  type(FmsNetcdfDomainFile_t)                  :: fileobj
  integer                                      :: ndims  !< Number of dimensions
  logical                                      :: gc_bug !< local grid_center_bug variable, default is .false.

  if(.not. open_file(fileobj, grid_file, 'read', domain )) then
     call mpp_error(FATAL, 'data_override_mod(get_grid_version_1): Error in opening file '//trim(grid_file))
  endif

  call mpp_get_data_domain(domain, isd, ied, jsd, jed)
  call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)

  select case(mod_name)
  case('ocn', 'ice')
    is_new_grid = .FALSE.
    if(variable_exists(fileobj, 'x_T')) then
       is_new_grid = .true.
    else if(variable_exists(fileobj, 'geolon_t')) then
       is_new_grid = .FALSE.
    else
       call mpp_error(FATAL,'data_override: both x_T and geolon_t is not in the grid file '//trim(grid_file) )
    endif

    if(is_new_grid) then
      ndims = get_variable_num_dimensions(fileobj, 'x_T')
      call get_variable_size(fileobj, 'x_T', siz(1:ndims))
      nlon = siz(1); nlat = siz(2)
      call check_grid_sizes(trim(mod_name)//'_domain  ', domain, nlon, nlat)
      select type (lon)
      type is (real(r4_kind))
         select type (lat)
         type is (real(r4_kind))
            allocate(lon_vert_4(isc:iec,jsc:jec,4), lat_vert_4(isc:iec,jsc:jec,4) )

            call read_data(fileobj, 'x_vert_T', lon_vert_4)
            call read_data(fileobj, 'y_vert_T', lat_vert_4)

!2 Global lon and lat of ocean grid cell centers are determined from adjacent vertices
            lon(:,:) = (lon_vert_4(:,:,1) + lon_vert_4(:,:,2) + lon_vert_4(:,:,3) + lon_vert_4(:,:,4))*0.25
            lat(:,:) = (lat_vert_4(:,:,1) + lat_vert_4(:,:,2) + lat_vert_4(:,:,3) + lat_vert_4(:,:,4))*0.25
            deallocate(lon_vert_4)
            deallocate(lat_vert_4)
         end select
      type is (real(r8_kind))
         select type (lat)
         type is (real(r8_kind))
            allocate(lon_vert_8(isc:iec,jsc:jec,4), lat_vert_8(isc:iec,jsc:jec,4) )

            call read_data(fileobj, 'x_vert_T', lon_vert_8)
            call read_data(fileobj, 'y_vert_T', lat_vert_8)

!2 Global lon and lat of ocean grid cell centers are determined from adjacent vertices
            lon(:,:) = (lon_vert_8(:,:,1) + lon_vert_8(:,:,2) + lon_vert_8(:,:,3) + lon_vert_8(:,:,4))*0.25
            lat(:,:) = (lat_vert_8(:,:,1) + lat_vert_8(:,:,2) + lat_vert_8(:,:,3) + lat_vert_8(:,:,4))*0.25
            deallocate(lon_vert_8)
            deallocate(lat_vert_8)
         end select
      end select
    else

      if (present(grid_center_bug)) then
          gc_bug = grid_center_bug
      else
          gc_bug = .false.
      endif

      if(gc_bug) call mpp_error(NOTE, &
           'data_override: grid_center_bug is set to true, the grid center location may be incorrect')

      ndims = get_variable_num_dimensions(fileobj, 'geolon_vert_t')
      call get_variable_size(fileobj, 'geolon_vert_t', siz(1:ndims))
      nlon = siz(1) - 1; nlat = siz(2) - 1;
      call check_grid_sizes(trim(mod_name)//'_domain  ', domain, nlon, nlat)

      start(1) = isc; nread(1) = iec-isc+2
      start(2) = jsc; nread(2) = jec-jsc+2

      select type (lon)
      type is (real(r4_kind))
         select type (lat)
         type is (real(r4_kind))
            allocate(lon_vert_4(isc:iec+1,jsc:jec+1,1))
            allocate(lat_vert_4(isc:iec+1,jsc:jec+1,1))

            call read_data(fileobj, 'geolon_vert_t', lon_vert_4(:,:,1), corner=start, edge_lengths=nread)
            call read_data(fileobj, 'geolat_vert_t', lat_vert_4(:,:,1), corner=start, edge_lengths=nread)

            if(gc_bug) then
               do j = jsc, jec
                  do i = isc, iec
                     lon(i,j) = (lon_vert_4(i,j,1) + lon_vert_4(i+1,j,1))/2.
                     lat(i,j) = (lat_vert_4(i,j,1) + lat_vert_4(i,j+1,1))/2.
                  enddo
               enddo
            else
               do j = jsc, jec
                  do i = isc, iec
                     lon(i,j) = (lon_vert_4(i,j,1) + lon_vert_4(i+1,j,1) + &
                          lon_vert_4(i+1,j+1,1) + lon_vert_4(i,j+1,1))*0.25
                     lat(i,j) = (lat_vert_4(i,j,1) + lat_vert_4(i+1,j,1) + &
                          lat_vert_4(i+1,j+1,1) + lat_vert_4(i,j+1,1))*0.25
                  enddo
               enddo
            end if
            deallocate(lon_vert_4)
            deallocate(lat_vert_4)
         end select
      type is (real(r8_kind))
         select type (lat)
         type is (real(r8_kind))
            allocate(lon_vert_8(isc:iec+1,jsc:jec+1,1))
            allocate(lat_vert_8(isc:iec+1,jsc:jec+1,1))

            call read_data(fileobj, 'geolon_vert_t', lon_vert_8(:,:,1), corner=start, edge_lengths=nread)
            call read_data(fileobj, 'geolat_vert_t', lat_vert_8(:,:,1), corner=start, edge_lengths=nread)

            if(gc_bug) then
               do j = jsc, jec
                  do i = isc, iec
                     lon(i,j) = (lon_vert_8(i,j,1) + lon_vert_8(i+1,j,1))/2.
                     lat(i,j) = (lat_vert_8(i,j,1) + lat_vert_8(i,j+1,1))/2.
                  enddo
               enddo
            else
               do j = jsc, jec
                  do i = isc, iec
                     lon(i,j) = (lon_vert_8(i,j,1) + lon_vert_8(i+1,j,1) + &
                          lon_vert_8(i+1,j+1,1) + lon_vert_8(i,j+1,1))*0.25
                     lat(i,j) = (lat_vert_8(i,j,1) + lat_vert_8(i+1,j,1) + &
                          lat_vert_8(i+1,j+1,1) + lat_vert_8(i,j+1,1))*0.25
                  enddo
               enddo
            end if
            deallocate(lon_vert_8)
            deallocate(lat_vert_8)
         end select
      end select
    endif
  case('atm', 'lnd')
     if(trim(mod_name) == 'atm') then
        xname = 'xta'; yname = 'yta'
     else
        xname = 'xtl'; yname = 'ytl'
     endif
     select type (lon)
     type is (real(r4_kind))
        select type (lat)
        type is (real(r4_kind))
           ndims = get_variable_num_dimensions(fileobj, xname)
           call get_variable_size(fileobj, xname, siz(1:ndims))
           nlon = siz(1); allocate(glon_4(nlon))
           call read_data(fileobj, xname, glon_4)

           ndims = get_variable_num_dimensions(fileobj, xname)
           call get_variable_size(fileobj, yname, siz(1:ndims))
           nlat = siz(1); allocate(glat_4(nlat))
           call read_data(fileobj, yname, glat_4)
           call check_grid_sizes(trim(mod_name)//'_domain  ', domain, nlon, nlat)

           is = isc - isg + 1; ie = iec - isg + 1
           js = jsc - jsg + 1; je = jec - jsg + 1
           do j = js, jec
              do i = is, ie
                 lon(i,j) = glon_4(i)
                 lat(i,j) = glat_4(j)
              enddo
           enddo
           deallocate(glon_4)
           deallocate(glat_4)
        end select
     type is (real(r8_kind))
        select type (lat)
        type is (real(r8_kind))
           ndims = get_variable_num_dimensions(fileobj, xname)
           call get_variable_size(fileobj, xname, siz(1:ndims))
           nlon = siz(1); allocate(glon_8(nlon))
           call read_data(fileobj, xname, glon_8)

           ndims = get_variable_num_dimensions(fileobj, xname)
           call get_variable_size(fileobj, yname, siz(1:ndims))
           nlat = siz(1); allocate(glat_8(nlat))
           call read_data(fileobj, yname, glat_8)
           call check_grid_sizes(trim(mod_name)//'_domain  ', domain, nlon, nlat)

           is = isc - isg + 1; ie = iec - isg + 1
           js = jsc - jsg + 1; je = jec - jsg + 1
           do j = js, jec
              do i = is, ie
                 lon(i,j) = glon_8(i)
                 lat(i,j) = glat_8(j)
              enddo
           enddo
           deallocate(glon_8)
           deallocate(glat_8)
        end select
     end select
  case default
     call mpp_error(FATAL, "data_override_mod: mod_name should be 'atm', 'ocn', 'ice' or 'lnd' ")
  end select

  call close_file(fileobj)

  ! convert from degree to radian
  select type (lon)
  type is (real(r4_kind))
     select type (lat)
     type is (real(r4_kind))
        select type (min_lon)
        type is (real(r4_kind))
           select type (max_lon)
           type is (real(r4_kind))
              lon = lon * deg_to_radian_4
              lat = lat * deg_to_radian_4
              min_lon = minval(lon)
              max_lon = maxval(lon)
              call mpp_min(min_lon)
              call mpp_max(max_lon)
           end select
        end select
     end select
  type is (real(r8_kind))
     select type (lat)
     type is (real(r8_kind))
        select type (min_lon)
        type is (real(r8_kind))
           select type (max_lon)
           type is (real(r8_kind))
              lon = lon * deg_to_radian_8
              lat = lat * deg_to_radian_8
              min_lon = minval(lon)
              max_lon = maxval(lon)
              call mpp_min(min_lon)
              call mpp_max(max_lon)
           end select
        end select
     end select
  end select


end subroutine get_grid_version_1

!> Get global lon and lat of three model (target) grids from mosaic.nc.
!! Currently we assume the refinement ratio is 2 and there is one tile on each pe.
subroutine get_grid_version_2(fileobj, mod_name, domain, isc, iec, jsc, jec, lon, lat, min_lon, max_lon)
  type(FmsNetcdfFile_t),       intent(in) :: fileobj !< file object for grid file
  character(len=*),            intent(in) :: mod_name !< module name
  type(domain2d),              intent(in) :: domain !< 2D domain
  integer,                     intent(in) :: isc, iec, jsc, jec
  class(*), dimension(isc:,jsc:), intent(out) :: lon, lat
  class(*),                       intent(out) :: min_lon, max_lon

  integer            :: i, j, siz(2)
  integer            :: nlon, nlat             ! size of global grid
  integer            :: nlon_super, nlat_super ! size of global supergrid.
  integer            :: isd, ied, jsd, jed
  integer            :: isg, ieg, jsg, jeg
  integer            :: isc2, iec2, jsc2, jec2
  character(len=256) :: solo_mosaic_file, grid_file
  real(r4_type), allocatable  :: tmpx_4(:,:), tmpy_4(:,:)
  real(r8_type), allocatable  :: tmpx_8(:,:), tmpy_8(:,:)
  type(domain2d)     :: domain2
  logical            :: open_solo_mosaic
  type(FmsNetcdfFile_t) :: mosaicfileobj, tilefileobj
  integer            :: start(2), nread(2)

  if(trim(mod_name) .NE. 'atm' .AND. trim(mod_name) .NE. 'ocn' .AND. &
     trim(mod_name) .NE. 'ice' .AND. trim(mod_name) .NE. 'lnd' ) call mpp_error(FATAL, &
        "data_override_mod: mod_name should be 'atm', 'ocn', 'ice' or 'lnd' ")

  call mpp_get_data_domain(domain, isd, ied, jsd, jed)
  call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)

  ! get the grid file to read

  if(variable_exists(fileobj, trim(mod_name)//'_mosaic_file' )) then
     call read_data(fileobj, trim(mod_name)//'_mosaic_file', solo_mosaic_file)

     solo_mosaic_file = 'INPUT/'//trim(solo_mosaic_file)
     if(.not. open_file(mosaicfileobj, solo_mosaic_file, 'read')) then
        call mpp_error(FATAL, 'data_override_mod(get_grid_version_2: Error in opening solo mosaic file '//trim(solo_mosaic_file))
     endif
     open_solo_mosaic=.true.
  else
     mosaicfileobj = fileobj
     open_solo_mosaic = .false.
  end if

  call get_mosaic_tile_grid(grid_file, mosaicfileobj, domain)

  if(.not. open_file(tilefileobj, grid_file, 'read')) then
     call mpp_error(FATAL, 'data_override_mod(get_grid_version_2: Error in opening tile file '//trim(grid_file))
  endif

  call get_variable_size(tilefileobj, 'area', siz)
  nlon_super = siz(1); nlat_super = siz(2)
  if( mod(nlon_super,2) .NE. 0) call mpp_error(FATAL,  &
       'data_override_mod: '//trim(mod_name)//' supergrid longitude size can not be divided by 2')
  if( mod(nlat_super,2) .NE. 0) call mpp_error(FATAL,  &
       'data_override_mod: '//trim(mod_name)//' supergrid latitude size can not be divided by 2')
  nlon = nlon_super/2;
  nlat = nlat_super/2;
  call check_grid_sizes(trim(mod_name)//'_domain  ', domain, nlon, nlat)
  isc2 = 2*isc-1; iec2 = 2*iec+1
  jsc2 = 2*jsc-1; jec2 = 2*jec+1

  start(1) = isc2; nread(1) = iec2-isc2+1
  start(2) = jsc2; nread(2) = jec2-jsc2+1

  select type (lon)
  type is (real(r4_kind))
     select type (lat)
     type is (real(r4_kind))
        allocate(tmpx_4(isc2:iec2, jsc2:jec2), tmpy_4(isc2:iec2, jsc2:jec2) )

        call read_data( tilefileobj, 'x', tmpx_4, corner=start,edge_lengths=nread)
        call read_data( tilefileobj, 'y', tmpy_4, corner=start,edge_lengths=nread)

        ! copy data onto model grid
        if(trim(mod_name) == 'ocn' .OR. trim(mod_name) == 'ice') then
           do j = jsc, jec
              do i = isc, iec
                 lon(i,j) = (tmpx_4(i*2-1,j*2-1)+tmpx_4(i*2+1,j*2-1)+tmpx_4(i*2+1,j*2+1)+tmpx_4(i*2-1,j*2+1))*0.25
                 lat(i,j) = (tmpy_4(i*2-1,j*2-1)+tmpy_4(i*2+1,j*2-1)+tmpy_4(i*2+1,j*2+1)+tmpy_4(i*2-1,j*2+1))*0.25
              end do
           end do
        else
           do j = jsc, jec
              do i = isc, iec
                 lon(i,j) = tmpx_4(i*2,j*2)
                 lat(i,j) = tmpy_4(i*2,j*2)
              end do
           end do
        endif

        ! convert to radian
        lon = lon * deg_to_radian_4
        lat = lat * deg_to_radian_4

        deallocate(tmpx_4, tmpy_4)
        select type (min_lon)
        type is (real(r4_kind))
           select type (max_lon)
           type is (real(r4_kind))
              min_lon = minval(lon)
              max_lon = maxval(lon)
              call mpp_min(min_lon)
              call mpp_max(max_lon)
           end type
        end type
     end type
  type is (real(r8_kind))
     select type (lat)
     type is (real(r8_kind))
        allocate(tmpx_8(isc2:iec2, jsc2:jec2), tmpy_8(isc2:iec2, jsc2:jec2) )

        call read_data( tilefileobj, 'x', tmpx_8, corner=start,edge_lengths=nread)
        call read_data( tilefileobj, 'y', tmpy_8, corner=start,edge_lengths=nread)

        ! copy data onto model grid
        if(trim(mod_name) == 'ocn' .OR. trim(mod_name) == 'ice') then
           do j = jsc, jec
              do i = isc, iec
                 lon(i,j) = (tmpx_8(i*2-1,j*2-1)+tmpx_8(i*2+1,j*2-1)+tmpx_8(i*2+1,j*2+1)+tmpx_8(i*2-1,j*2+1))*0.25
                 lat(i,j) = (tmpy_8(i*2-1,j*2-1)+tmpy_8(i*2+1,j*2-1)+tmpy_8(i*2+1,j*2+1)+tmpy_8(i*2-1,j*2+1))*0.25
              end do
           end do
        else
           do j = jsc, jec
              do i = isc, iec
                 lon(i,j) = tmpx_8(i*2,j*2)
                 lat(i,j) = tmpy_8(i*2,j*2)
              end do
           end do
        endif

        ! convert to radian
        lon = lon * deg_to_radian_8
        lat = lat * deg_to_radian_8

        deallocate(tmpx_8, tmpy_8)
        select type (min_lon)
        type is (real(r8_kind))
           select type (max_lon)
           type is (real(r8_kind))
              min_lon = minval(lon)
              max_lon = maxval(lon)
              call mpp_min(min_lon)
              call mpp_max(max_lon)
           end type
        end type
     end type

  call close_file(tilefileobj)
  if(open_solo_mosaic)  call close_file(mosaicfileobj)

end subroutine get_grid_version_2

end module get_grid_version_mod
!> @}
! close documentation grouping
