subroutine FCT(A,Bx,C,ux,dsx,nx,dt,Ajnew)
    implicit none
    integer :: nx
    integer :: minusOne,minusTwo
    integer, parameter :: DP=selected_real_kind(8)
    real(kind=DP),intent(in)    ::A(nx),Bx(nx),C(nx)
    real(kind=DP),intent(in)    ::ux(nx)
!   Define Variables A-Energy/Momentum, B-Source term,     
!   ux,uy-velocity, dsx,dsy-CFL number
    real(kind=DP),intent(in)    ::dsx,dt
    real(kind=DP),intent(out)    ::Ajnew(nx)
    real(kind=DP),dimension(0:nx-1)  ::   AjnewPrimecopy
    real(kind=DP),dimension(0:nx-1)  ::     Acopy,Bxcopy,Ccopy,uxcopy,Ajnewcopy
    real(kind=DP),dimension(0:nx-1)  ::     Qplx,  Qmnx, Diff
    real(kind=DP),dimension(0:nx-1)  ::     Abarjx,Sx
    real(kind=DP),dimension(0:nx-1)  ::     Abarj,Deltajx
    real(kind=DP),dimension(0:nx-1)  ::    Ajexx,AjTildex
    real(kind=DP),dimension(0:nx-1)  ::    ATildeMin,ATildeMax, Ajin, Ajout
    real(kind=DP),dimension(0:nx-1)  ::     Fin,Fout,Ahatx
    real(kind=DP),dimension(0:nx-1)  ::    epsilonx
    real(kind=DP),dimension(0:nx-1)  ::     ZeroArray, OneArray
    
    
    real(kind=DP)  ::   cut=10.0_DP**(-10)
    real(kind=DP)  ::   Difsn=1.0_DP/8.0_DP
!    #Difsn=1/8                      #Diffusion coeffiecient

    real(kind=DP)  ::   Amskx=1.0_DP/4.0_DP, Amsky=0.25_DP
    ZeroArray=0.0_DP
    OneArray=1.0_DP

!   Assign all elements to zero   

    
    
    Ajnewcopy=0.0_DP
    AjnewPrimecopy=0.0_DP
    Acopy=0.0_DP
    Bxcopy=0.0_DP
    Ccopy=0.0_DP
    uxcopy=0.0_DP
 
    Qplx=1.0_DP
    Qmnx=1.0_DP
    Abarjx=1.0_DP
    Sx=0.0_DP
    Abarj=1.0_DP
    Deltajx=1.0_DP
    Ajexx=0.0_DP
    AjTildex=1.0_DP
    ATildeMin=0.0_DP
    ATildeMax=0.0_DP
    Ajin=0.0_DP
    Ajout=0.0_DP
    Fin=0.0_DP
    Fout=0.0_DP
    Ahatx=0.0_DP
    epsilonx=1.0_DP
    minusOne=nx-1
    minusTwo=nx-2
    
    
!    Acopy(:nx-1,:ny-1)=A 
    Acopy=A 
    Bxcopy=Bx
    Ccopy=C
    uxcopy=ux
    
    
    epsilonx(:)=uxcopy(:)*dsx  
    Diff(:minusOne)=Acopy(1:)-Acopy(:minusOne)       
 !  Calculating epsilon
    
   
    
    Qplx(:minusOne)=(0.5_DP - epsilonx(:minusOne))/(1+(epsilonx(1:)-epsilonx(:minusOne)))
    
    Qmnx(1:)=(0.5_DP + epsilonx(1:))/(1-(epsilonx(:minusOne)-epsilonx(1:)))
    
!    #Calculating transported and diffused quantitites, with source term treated using
!    #original explicit SHASTA proposed by Boris and Book
    
!    #Abarj(1:minusOne)=1/2*Qpl(1:minusOne)**2*(A(2:)-A(1:minusOne))-1/2*Qmn(1:minusOne)**2*(A(1:minusOne)-A(:minusTwo)) \
!    #+Qpl(1:minusOne)*(A(1:minusOne)+(B(2:)-B(1:minusOne))*ds)+Qmn(1:minusOne)*(A(1:minusOne)+(B(1:minusOne)-B(:minusTwo))*ds)
    
!    #Calculating transported and diffused quantitites, with Rischke's method
  
    Abarj(1:minusOne)= 0.5_DP*Qplx(1:minusOne)**2*(Acopy(2:)-Acopy(1:minusOne))-0.5_DP*&
                   Qmnx(1:minusOne)**2*(Acopy(1:minusOne)-Acopy(:minusTwo)) +&
                   (Qplx(1:minusOne)+Qmnx(1:minusOne))*Acopy(1:minusOne)-(Bxcopy(2:)-Bxcopy(:minusTwo))*dsx/2.0_DP-(Ccopy(:))*dt
    
!    #Calculating the amount of diffusion that leads to the antidiffusion fluxes using either
!    #explicit or phonenical SHASTA
    Deltajx(:minusOne)=(Abarj(1:)-Abarj(:minusOne))
        
!    #Explicit SHASTA with mask coefficient
!    Ajexx(:)=Amskx*Deltajx(:)
    
!    #Phoenical SHASTA
    Ajexx(1:minusOne)=Difsn*1.0_DP/8.0_DP*(Deltajx(1:minusOne)-1.0_DP/8.0_DP*(Diff(2:)-2*Diff(1:minusOne)+Diff(:minusTwo)))
    
!    #Calculating the “flux corrected” antidiffusion fluxes
    AjTildex(1:minusOne)=sign(OneArray(1:minusOne), Deltajx(1:minusOne))*max(ZeroArray(1:minusOne),min(Deltajx(2:)*&
                     sign(OneArray(1:minusOne),Ajexx(1:minusOne)), dabs(Ajexx(1:minusOne)),Deltajx(:minusTwo)*&
                     sign(OneArray(1:minusOne),Ajexx(1:minusOne))) )
    
                  
    
!    #Calculating the final time-advanced quantities
    Ajnewcopy(1:)=Abarj(1:)-AjTildex(1:)+AjTildex(:minusOne)

                        
!    #Handling boundary values

    AjnewPrimeCopy=Ajnewcopy(:minusOne)



     AjnewPrimecopy(1)=AjnewPrimecopy(2); AjnewPrimecopy(0)=AjnewPrimecopy(2)
     AjnewPrimecopy(minusOne)=AjnewPrimecopy(minusTwo-1); AjnewPrimecopy(minusTwo)=AjnewPrimecopy(minusTwo-1)
!   Ajnewcopy(:,minusTwo)=Ajnewcopy(:,minusTwo-2); Ajnewcopy(:,minusOne)=Ajnewcopy(:,minusTwo-2)
!   Ajnewcopy(:,minusOne-1)=Ajnewcopy(:,minusTwo-2)
!    Ajnewcopy(nx,ny)=Ajnewcopy(nx-2,ny-2)
    Ajnew=AjnewPrimecopy


!    #print("Ajnew",Ajnew(:,3:4))

end subroutine FCT

