subroutine velocityfinder1D(E,Mx,elab,nx,uxnew,vnew)
    implicit none
    integer :: nx,ny
    integer, parameter:: DP=selected_real_kind(8)
    real(kind=DP),intent(in)::E(nx),Mx(nx)
    real(kind=DP),intent(out)::elab(nx),uxnew(nx)
    real(kind=DP),intent(out)::vnew
    real(kind=DP), allocatable :: ecopy(:),mxcopy(:),unewabs(:),mabs(:)
    real(kind=DP), parameter :: toler = 1.0_DP-10.0_DP**(-8)
    real(kind=DP), parameter :: cut = 10.0_DP**(-12)
    real(kind=DP)            :: cx=0.0_DP,cy=0.0_DP
    integer :: i
    integer :: j

    allocate(ecopy(nx))
    allocate(mxcopy(nx))
    !allocate(mycopy(nx,ny))
    allocate(unewabs(nx))
    allocate(mabs(nx))         
   
    do i=1,nx
            mabs(i)=dabs((Mx(i)**2)**0.5)
            if (mabs(i) >= dabs(E(i))) then
                cx=Mx(i)/mabs(i)
                !cy=My(i,j)/mabs(i,j)
                mabs(i)=dabs(E(i))*toler
                mxcopy(i)=cx*mabs(i)
                !mycopy(i,j)=cy*mabs(i,j)
            else
            mabs(i)=mabs(i)
            mxcopy(i)=Mx(i)
            !mycopy(i,j)=My(i,j)
            endif
        end do
            
    do i=1,nx
            mabs(i)=dabs((mxcopy(i)**2)**0.5)
            if (dabs(E(i)) <= cut) then
                ecopy(i)=0.0_DP
                uxnew(i)=0.0_DP
                !uynew(i)=0.0_DP
                elab(i)=0.0_DP
            else if (mabs(i) <= cut) then
                ecopy(i)=E(i)
                uxnew(i)=0.0_DP
                !uynew(i)=0.0_DP
                elab(i)=dabs(E(i))
            else
                call Root(dabs(E(i)),dabs(mabs(i)),vnew)
                uxnew(i)=(vnew*mxcopy(i))/dabs(mabs(i))
                !uynew(i)=(vnew*mycopy(i))/dabs(mabs(i))               
                elab(i)=dabs(E(i))-dabs(mabs(i))*vnew
            endif
   enddo
    deallocate(ecopy)
    deallocate(mxcopy)
    !deallocate(mycopy)
    deallocate(unewabs)
    deallocate(mabs) 
end subroutine velocityfinder1D


subroutine Root(En,M,vnew)
    implicit none
    integer, parameter:: DP=selected_real_kind(8)
    real(kind=DP),intent(in)::En,M
    real(kind=DP),intent(out)::vnew
    integer :: n
    integer :: it
    real(kind=DP):: vtemp=0.0_DP
    real(kind=DP):: v=0.0_DP
    real(kind=DP), parameter :: tol = 10.0_DP**(-20)
    n=1
    it=10000
    call func(v,En,M,vnew)
    vtemp=vnew
    do while ((n<it).AND.(dabs(v-vtemp)>tol))
        call func(v,En,M,vnew)
        v=vnew
        call func(v,En,M,vnew)
        vtemp=vnew
        n=n+1
    end do

end subroutine Root

subroutine func(vold,En,M,vnew)
    implicit none
    integer, parameter:: DP=selected_real_kind(8)
    real(kind=DP),intent(in)::En,M,vold
    real(kind=DP),intent(out)::vnew
    real(kind=DP) ::val=0.0_DP
!    call eos(En - M*vold,val)
    val=1/3.0_DP*(En - M*vold)
    vnew= M/(En + val)
    
end subroutine
