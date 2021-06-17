module freesurface

  use decomp_2d, only : mytype, xsize

  implicit none

  real(mytype), save, allocatable, dimension(:,:,:) :: nx1, ny1, nz1
  real(mytype), save, allocatable, dimension(:,:,:) :: nx1_old, ny1_old, nz1_old
  real(mytype), save, allocatable, dimension(:,:,:) :: norm_mag_s, norm_mag

  real(mytype), save, allocatable, dimension(:,:,:) :: nx1_2, ny1_2, nz1_2
  real(mytype), save, allocatable, dimension(:,:,:) :: nx1_old_2, ny1_old_2, nz1_old_2

  real(mytype), save, allocatable, dimension(:,:,:) :: pot, bulk
  real(mytype), save, allocatable, dimension(:,:,:) :: wx1, wy1, wz1

  private
  public :: update_fluid_properties, &
       surface_tension_force, source_pf, speed_max

contains

  subroutine update_fluid_properties(rho1, mu1, levelset1)

    use decomp_2d, only : mytype, xsize, ysize, zsize
    use var, only : half, one, two, three, pi, ten, sixteen
    use var, only : dx, dy, dz
    use var, only : numscalar
    use param, only : nrhotime, ilevelset, xlx
    use param, only : dens1, dens2, visc1, visc2
    use var, only : nclx1, ncly1, nclz1, nclxn, nclyn, nclzn
    use var, only : nclxS1, nclyS1, nclzS1, nclxSn, nclySn, nclzSn

    use var, only : lsp1 => pgy1, lspv21 => pp1
    use var, only : lsp12 => uyp2, lsp2 => upi2, lspv32 => pp2, lspv2 => ppi2
    use var, only : lsp23 => uzp3, lsp3 => po3, lspv3 => ppi3, prop3 => dv3

    use var, only : di1, sx, cifxp6, cisxp6, ciwxp6, nxmsize
    use var, only : dipp2, sy, cifyp6, cisyp6, ciwyp6, nymsize
    use var, only : dipp3, sz, cifzp6, ciszp6, ciwzp6, nzmsize
    use var, only : dip3, cifip6z, cisip6z, ciwip6z, cifz6, cisz6, ciwz6
    use var, only : dip2, cifip6y, cisip6y, ciwip6y, cify6, cisy6, ciwy6
    use var, only : cifip6, cisip6, ciwip6, cifx6, cisx6, ciwx6
    use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x

    use var, only : ph1, ph2, ph3, ph4

    implicit none

    !! Input
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: levelset1

    !! InOut
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: mu1

    !! Local
    integer :: i, j, k



       rho1(:,:,:,1) = levelset1(:, :, :)*dens1 + ((one-levelset1(:, :, :))*dens2)

       mu1(:,:,:) = levelset1(:, :, :)*visc1 + ((one-levelset1(:, :, :))*visc2)


  endsubroutine update_fluid_properties

  subroutine surface_tension_force(dux1, duy1, duz1, rho1, levelset1)


    use decomp_2d, only : mytype
    use decomp_2d, only : xsize, ysize, zsize
    use decomp_2d, only : xstart, xend
    use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x
    use decomp_2d, only : alloc_x


    use variables, only : derxS, deryS, derzS
    use variables, only : sx, sy, sz

    use var, only : ep1, divu3, drho1

    use var, only : zero, half, one, two, three, pi, four

    use var, only : di1
    use var, only : ny, nz

    use variables, only : ffz, fsz, fwz, ffy, fsy, fwy, ffx, fsx, fwx

    use var, only : di2
    use var, only : di3

    use param, only : nclx1, ncly1, nclz1
    use param, only : dx, dy, dz, xlx, yly, zlz
    use param, only : nrhotime
    use param, only : dens1, dens2
    use param, only : itime, ifirst, ntime
    use param, only : six

    use var, only : nxmsize, nymsize, nzmsize
    use var, only : ph1, ph2, ph3, ph4
    use var, only : nclx1, ncly1, nclz1, nclxn, nclyn, nclzn
    use var, only : ux1, uy1, uz1
    use var, only : ux2, uy2, uz2
    use var, only : ux3, uy3, uz3

    USE var, ONLY : tc1,td1, te1, tf1, tg1, di1
    USE var, ONLY : tc2,td2, te2, tf2, tg2, di2
    USE var, ONLY : tc3,td3, te3, tf3, tg3, di3
    use var, only : nclxS1, nclyS1, nclzS1, nclxSn, nclySn, nclzSn
    use variables, only : ffxpS, fsxpS, fwxpS
    use variables, only : ffypS, fsypS, fwypS, ppy
    use variables, only : ffzpS, fszpS, fwzpS

    use weno, only : weno5

    implicit none

    !! IN
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: levelset1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime), intent(in) :: rho1

    !! INOUT
    !real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(inout) :: sx1, sy1, sz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1

    !! LOCAL
    real(mytype) :: normmag
    integer :: nlock
    integer :: i, j, k
    real(mytype) :: sigma ! Surface tension coefficient
    real(mytype) :: x, y, z
    real(mytype) :: ddelta, eta, alpha_pf, rhomean

    sigma = 0.00001_mytype
    eta = xlx/(xsize(1) - one)
    alpha_pf = six*sqrt(two)
    rhomean = (dens1 + dens2) / two


    if (itime.eq.ifirst) then
       call alloc_x(nx1); call alloc_x(ny1); call alloc_x(nz1)
       call alloc_x(nx1_old); call alloc_x(ny1_old); call alloc_x(nz1_old)
       call alloc_x(norm_mag_s)
    endif

    call transpose_x_to_y(levelset1(:,:,:),td2(:,:,:))
    call transpose_y_to_z(td2(:,:,:),td3(:,:,:))

    call weno5(te1, levelset1(:,:,:), ux1, 1, nclxS1, nclxSn, xsize(1), xsize(2), xsize(3), 1)

    call weno5(tf2,td2(:,:,:), uy2, 2, nclyS1, nclySn, ysize(1), ysize(2), ysize(3), 1)

    call transpose_y_to_x(tf2(:,:,:),tf1(:,:,:))

    call weno5(tg3,td3(:,:,:), uz3, 3, nclzS1, nclzSn, zsize(1), zsize(2), zsize(3), 1)

    call transpose_z_to_y(tg3(:,:,:),tg2(:,:,:))
    call transpose_y_to_x(tg2(:,:,:),tg1(:,:,:))

    norm_mag_s(:, :, :) = sqrt(te1(:,:,:)**2 + tf1(:,:,:)**2 + tg1(:,:,:)**2)

    nx1_old(:,:,:) = te1(:,:,:)
    ny1_old(:,:,:) = tf1(:,:,:)
    nz1_old(:,:,:) = tg1(:,:,:)

    do k = 1, xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)

             if (norm_mag_s(i,j,k).gt.zero) then
               nx1(i,j,k) = te1(i,j,k) / norm_mag_s(i,j,k)
               ny1(i,j,k) = tf1(i,j,k) / norm_mag_s(i,j,k)
               nz1(i,j,k) = tg1(i,j,k) / norm_mag_s(i,j,k)
             else
                nx1(i,j,k) = zero
                ny1(i,j,k) = zero
                nz1(i,j,k) = zero
             endif
          enddo
       enddo
    enddo

    call derxS (te1, nx1(:,:,:),di1,sx,ffxpS,fsxpS,fwxpS,xsize(1),xsize(2),xsize(3),1)


    call transpose_x_to_y(ny1(:,:,:),tc2(:,:,:))
    call deryS (tf2,tc2(:,:,:),di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1)
    call transpose_y_to_x(tf2(:,:,:),tf1(:,:,:))


    call transpose_x_to_y(nz1(:,:,:),tc2(:,:,:))
    call transpose_y_to_z(tc2(:,:,:),tc3(:,:,:))
    call derzS (tg3,tc3(:,:,:),di3,sz,ffzpS,fszpS,fwzpS,zsize(1),zsize(2),zsize(3),1)
    call transpose_z_to_y(tg3(:,:,:),tg2(:,:,:))
    call transpose_y_to_x(tg2(:,:,:),tg1(:,:,:))


    !! Add contribution to forcing terms
    dux1(:,:,:,1) = dux1(:,:,:,1) - (sigma * (te1(:,:,:) + tf1(:,:,:) + tg1(:,:,:)) * nx1_old(:,:,:))
    duy1(:,:,:,1) = duy1(:,:,:,1) - (sigma * (te1(:,:,:) + tf1(:,:,:) + tg1(:,:,:)) * ny1_old(:,:,:))
    duz1(:,:,:,1) = duz1(:,:,:,1) - (sigma * (te1(:,:,:) + tf1(:,:,:) + tg1(:,:,:)) * nz1_old(:,:,:))

  end subroutine surface_tension_force

  subroutine source_pf(sx1, rho1, levelset1)

    use decomp_2d, only : mytype
    use decomp_2d, only : xsize, ysize, zsize
    use decomp_2d, only : xstart, xend
    use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x
    use decomp_2d, only : alloc_x


    use variables, only : derxS, deryS, derzS
    use variables, only : sx, sy, sz

    use var, only : ep1, divu3, drho1

    use var, only : zero, half, one, two, three, pi, four

    use var, only : di1
    use var, only : ny, nz

    use variables, only : ffz, fsz, fwz, ffy, fsy, fwy, ffx, fsx, fwx

    use var, only : di2
    use var, only : di3

    use param, only : nclx1, ncly1, nclz1
    use param, only : dx, dy, dz, xlx, yly, zlz
    use param, only : nrhotime
    use param, only : dens1, dens2
    use param, only : itime, ifirst
    use param, only : six

    use var, only : nxmsize, nymsize, nzmsize
    use var, only : ph1, ph2, ph3, ph4
    use var, only : nclx1, ncly1, nclz1, nclxn, nclyn, nclzn
    use var, only : ux1, uy1, uz1
    use var, only : ux2, uy2, uz2
    use var, only : ux3, uy3, uz3

    USE var, ONLY : tc1,td1, te1, tf1, tg1, di1
    USE var, ONLY : tc2,td2, te2, tf2, tg2, di2
    USE var, ONLY : tc3,td3, te3, tf3, tg3, di3
    use var, only : nclxS1, nclyS1, nclzS1, nclxSn, nclySn, nclzSn
    use variables, only : ffxpS, fsxpS, fwxpS
    use variables, only : ffypS, fsypS, fwypS, ppy
    use variables, only : ffzpS, fszpS, fwzpS

    use weno, only : weno5

    implicit none

    !! IN
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: levelset1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime), intent(in) :: rho1

    !! INOUT
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(inout) :: sx1

    !! LOCAL
    real(mytype) :: normmag
    integer :: nlock
    integer :: i, j, k
    real(mytype) :: x, y, z
    real(mytype) :: ddelta, eta, rhomean, ux_max,uy_max,uz_max, ux_max0,uy_max0,uz_max0
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: vel_mag


    eta = xlx/(xsize(1) - one)

    rhomean = (dens1 + dens2) / two


    if (itime.eq.ifirst) then
       call alloc_x(nx1_2); call alloc_x(ny1_2); call alloc_x(nz1_2)
       call alloc_x(nx1_old_2); call alloc_x(ny1_old_2); call alloc_x(nz1_old_2)
       call alloc_x(norm_mag)

    endif

    call transpose_x_to_y(levelset1(:,:,:),td2(:,:,:))
    call transpose_y_to_z(td2(:,:,:),td3(:,:,:))

    call derxS (te1, levelset1(:,:,:),di1,sx,ffxpS,fsxpS,fwxpS,xsize(1),xsize(2),xsize(3),1)

    call deryS (tf2,td2(:,:,:),di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1)
    call transpose_y_to_x(tf2(:,:,:),tf1(:,:,:))

    call derzS (tg3,td3(:,:,:),di3,sz,ffzpS,fszpS,fwzpS,zsize(1),zsize(2),zsize(3),1)
    call transpose_z_to_y(tg3(:,:,:),tg2(:,:,:))
    call transpose_y_to_x(tg2(:,:,:),tg1(:,:,:))

    norm_mag(:, :, :) = sqrt(te1(:,:,:)**2 + tf1(:,:,:)**2 + tg1(:,:,:)**2)

    nx1_old_2(:,:,:) = te1(:,:,:)
    ny1_old_2(:,:,:) = tf1(:,:,:)
    nz1_old_2(:,:,:) = tg1(:,:,:)

    do k = 1, xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)

             if (norm_mag(i,j,k).gt.zero) then
               nx1_2(i,j,k) = te1(i,j,k) / norm_mag(i,j,k)
               ny1_2(i,j,k) = tf1(i,j,k) / norm_mag(i,j,k)
               nz1_2(i,j,k) = tg1(i,j,k) / norm_mag(i,j,k)
             else
                nx1_2(i,j,k) = zero
                ny1_2(i,j,k) = zero
                nz1_2(i,j,k) = zero
             endif
          enddo
       enddo
    enddo

    call speed_max(ux1,uy1,uz1,ux_max0,uy_max0,uz_max0)

    ux_max = eta*ux_max0
    uy_max = eta*uy_max0
    uz_max = eta*uz_max0

    nx1_old_2(:,:,:) = ux_max * (eta*nx1_old_2(:,:,:))
    ny1_old_2(:,:,:) = uy_max * (eta*ny1_old_2(:,:,:))
    nz1_old_2(:,:,:) = uz_max * (eta*nz1_old_2(:,:,:))

    nx1_2(:,:,:) = ux_max * (((levelset1(:,:,:))*(one - levelset1(:,:,:))*(nx1_2(:,:,:))))
    ny1_2(:,:,:) = uy_max * (((levelset1(:,:,:))*(one - levelset1(:,:,:))*(ny1_2(:,:,:))))
    nz1_2(:,:,:) = uz_max * (((levelset1(:,:,:))*(one - levelset1(:,:,:))*(nz1_2(:,:,:))))

    call derxS (te1, nx1_old_2(:,:,:),di1,sx,ffxpS,fsxpS,fwxpS,xsize(1),xsize(2),xsize(3),1)


    call transpose_x_to_y(ny1_old_2(:,:,:),tc2(:,:,:))
    call deryS (tf2,tc2(:,:,:),di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1)
    call transpose_y_to_x(tf2(:,:,:),tf1(:,:,:))


    call transpose_x_to_y(nz1_old_2(:,:,:),tc2(:,:,:))
    call transpose_y_to_z(tc2(:,:,:),tc3(:,:,:))
    call derzS (tg3,tc3(:,:,:),di3,sz,ffzpS,fszpS,fwzpS,zsize(1),zsize(2),zsize(3),1)
    call transpose_z_to_y(tg3(:,:,:),tg2(:,:,:))
    call transpose_y_to_x(tg2(:,:,:),tg1(:,:,:))


    sx1(:,:,:) =  sx1(:,:,:) - (te1(:,:,:) + tf1(:,:,:) + tg1(:,:,:))

    call derxS (te1, nx1_2(:,:,:),di1,sx,ffxpS,fsxpS,fwxpS,xsize(1),xsize(2),xsize(3),1)


    call transpose_x_to_y(ny1_2(:,:,:),tc2(:,:,:))
    call deryS (tf2,tc2(:,:,:),di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1)
    call transpose_y_to_x(tf2(:,:,:),tf1(:,:,:))


    call transpose_x_to_y(nz1_2(:,:,:),tc2(:,:,:))
    call transpose_y_to_z(tc2(:,:,:),tc3(:,:,:))
    call derzS (tg3,tc3(:,:,:),di3,sz,ffzpS,fszpS,fwzpS,zsize(1),zsize(2),zsize(3),1)
    call transpose_z_to_y(tg3(:,:,:),tg2(:,:,:))
    call transpose_y_to_x(tg2(:,:,:),tg1(:,:,:))

    sx1(:,:,:) =  sx1(:,:,:) + (te1(:,:,:) + tf1(:,:,:) + tg1(:,:,:))

  end subroutine source_pf

  subroutine speed_max(ux,uy,uz,ux_max,uy_max,uz_max)

    USE decomp_2d
    USE decomp_2d_poisson
    USE variables
    USE param
    USE var
    USE MPI

    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ux,uy,uz
    real(mytype), intent(inout) :: ux_max,uy_max,uz_max

    integer :: code,ierror,i,j,k
    real(mytype) :: uxmax,uymax,uzmax,uxmin,uymin,uzmin,u_mag
    real(mytype) :: uxmax1,uymax1,uzmax1,uxmin1,uymin1,uzmin1


    uxmax=-1609.;uymax=-1609.;uzmax=-1609.;uxmin=1609.;uymin=1609.;uzmin=1609.

    do k = 1, xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)
            u_mag = sqrt((ux(i,j,k)**2)+(uy(i,j,k)**2)+(uz(i,j,k)**2))
                if (u_mag.gt.uxmax) then
                   uxmax=u_mag
                endif
          enddo
       enddo
    enddo

    call MPI_ALLREDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)


    ux_max = uxmax1
    uy_max = uxmax1
    uz_max = uxmax1


  end subroutine speed_max

endmodule freesurface
