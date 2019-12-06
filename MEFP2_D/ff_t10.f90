module t10
    use VAR
    implicit none
    
    contains
    
    !-----------------------------------------------------t10----------------------------------------------------!
    subroutine et10()
        if (allocated(coord_el)==.false.) then
            allocate(coord_el(10,2))
        else
        end if
        
        coord_el(:,1)=(/3.0d0,2.0d0,1.0d0,0.0d0,2.0d0,1.0d0,0.0d0,1.0d0,0.0d0,0.0d0/) !xsi1
        coord_el(:,2)=(/0.0d0,1.0d0,2.0d0,3.0d0,0.0d0,1.0d0,2.0d0,0.0d0,1.0d0,0.0d0/) !xsi2
        coord_el=coord_el/3.0d0
    
    end subroutine et10
    
    !-----------------------------------------------------fi-----------------------------------------------------!

    subroutine fi()
        if (allocated(fis)==.false.) then
            allocate(fis(10))
        else
        end if
        
        !xsi3 = 1. - xsi1 - xsi2 
        
        fis(1) = xsi1+4.5d0*(-1.d0+xsi1)*xsi1**2
        fis(2) = 4.5d0*xsi1*(-1.d0+3.d0*xsi1)*xsi2
        fis(3) = 4.5d0*xsi1*xsi2*(-1.d0+3.d0*xsi2)
        fis(4) = xsi2+4.5d0*(-1.d0+xsi2)*xsi2**2
        fis(5) = -4.5d0*xsi1*(-1.d0+3.d0*xsi1)*(-1.d0+xsi1+xsi2)
        fis(6) = -27.d0*xsi1*xsi2*(-1.d0+xsi1+xsi2)
        fis(7) = -4.5d0*xsi2*(-1.d0+xsi1+xsi2)*(-1.d0+3.d0*xsi2)
        fis(8) = 4.5d0*xsi1*(-1.d0+xsi1+xsi2)*(-2.d0+3.d0*xsi1+3.d0*xsi2)
        fis(9) = 4.5d0*xsi2*(-1.d0+xsi1+xsi2)*(-2.d0+3.d0*xsi1+3.d0*xsi2)
        fis(10) = -0.5d0*(-1.d0+xsi1+xsi2)*(-2.d0+3.d0*xsi1+3.d0*xsi2)*(-1.d0+3.d0*xsi1+3.d0*xsi2)
    
    end subroutine fi
    
    !-----------------------------------------------------dfi1----------------------------------------------------!
    
    subroutine dfi1()
        if (allocated(dfis1)==.false.) then
            allocate(dfis1(10))
        else
        end if
        
        dfis1(1) = 1.d0+4.5d0*xsi1*(-2.d0+3.d0*xsi1)
        dfis1(2) = 4.5d0*(-1.d0+6.d0*xsi1)*xsi2
        dfis1(3) = 4.5d0*xsi2*(-1.d0+3.d0*xsi2)
        dfis1(4) = 00.d0
        dfis1(5) = -4.5d0*(1.d0-1.d0*xsi2+xsi1*(-8.d0+9.d0*xsi1+6.d0*xsi2))
        dfis1(6) = -27.d0*xsi2*(-1.d0+2.d0*xsi1+xsi2)
        dfis1(7) = 4.5d0*(1.d0 -3.d0*xsi2)*xsi2
        dfis1(8) = 4.5d0*(2.d0+9.d0*xsi1**2-5.d0*xsi2+3.d0*xsi2**2+2.d0*xsi1*(-5.d0+6.d0*xsi2))
        dfis1(9) = 4.5d0*xsi2*(-5.d0+6.d0*xsi1+6.d0*xsi2)
        dfis1(10) = 0.5d0*(-11.d0+36.d0*xsi2-9.d0*(xsi1*(-4.d0+3.d0*xsi1)+6.d0*xsi1*xsi2+3.d0*xsi2**2))
    
    end subroutine dfi1
    
    !-----------------------------------------------------dfi2----------------------------------------------------!
    
    subroutine dfi2()
        if (allocated(dfis2)==.false.) then
            allocate(dfis2(10))
        else
        end if
        
        dfis2(1) = 0.d0
        dfis2(2) = 4.5d0*xsi1*(-1.d0+3.d0*xsi1)
        dfis2(3) = 4.5d0*xsi1*(-1.d0+6.d0*xsi2)
        dfis2(4) = 1.d0+4.5d0*xsi2*(-2.d0+3.d0*xsi2)
        dfis2(5) = 4.5d0*(1.d0-3.d0*xsi1)*xsi1
        dfis2(6) = -27.d0*xsi1*(-1.d0+xsi1+2.d0*xsi2)
        dfis2(7) = -4.5d0*(1.d0+xsi1*(-1.d0+6.d0*xsi2)+xsi2*(-8.d0+9.d0*xsi2))
        dfis2(8) = 4.5d0*xsi1*(-5.d0+6.d0*xsi1+6.d0*xsi2)
        dfis2(9) = 4.5*(2.d0+3.d0*xsi1**2+xsi2*(-10.d0+9.d0*xsi2)+xsi1*(-5.d0+12.d0*xsi2))
        dfis2(10) = 0.5d0*(-11.d0+36.d0*xsi2-9.d0*(xsi1*(-4.d0+3.d0*xsi1)+6.d0*xsi1*xsi2+3.d0*xsi2**2))
    
    end subroutine dfi2
    
    !---------------------------------------------------------------------------------------------------------!
    
end module
    