
    module leitura
    use VAR
    implicit none
    contains
    
    !-----------------------------------------------------acadm----------------------------------------------------!
    
    subroutine acadm()
    
        !Arquivo de saída do acadmesh
        open(10,file="acadmesh.txt", status = "old")
    
        read(10,*) texto !NODE X Y
    
        io=0
        do while (io==0)
            read(10,*,iostat=io) nnos
        end do
    
        read(10,*) texto !ELEM. NODES E v d THICK.
    
        io=0
        do while (io==0)
            read(10,*,iostat=io) nel
        end do
    
        read(10,*) texto !FORCE NODE FX FY
    
        io=0
        do while (io==0)
            read(10,*,iostat=io) nfnod
        end do
    
        read(10,*) texto !DISP. NODE DOF VALUE
    
        io=0
        do while (io==0)
            read(10,*,iostat=io) nbc
        end do
    
        rewind(10)
    
        call alocar()
    
        read(10,*) texto !NODE X Y
    
        do i=1,nnos
            read(10,*) auxi, X(auxi*2-1), X(auxi*2)
        end do
       
        read(10,*) texto !ELEM. NODES E v d THICK.
        
        do i=1,nel
            read(10,*) auxi, inc(auxi,1), inc(auxi,2), inc(auxi,3), inc(auxi,4), inc(auxi,5), inc(auxi,6), inc(auxi,7), &
                inc(auxi,8), inc(auxi,9), inc(auxi,10), prop(auxi,1), prop(auxi,2), prop(auxi,3), prop(auxi,4)
        end do
        
        read(10,*) texto !FORCE NODE FX FY
        
        fnodal=0.
        do i=1,nfnod
            read(10,*) texto, auxi, fnodal(auxi*2-1), fnodal(auxi*2)
        end do
                
        read(10,*) texto !DISP. NODE DOF VALUE
        
        rest = 0
        dnodal = 0.
        do i=1,nbc
            read(10,*) texto, auxi, texto, auxr
            if (texto=='X') then
                rest(auxi*2-1) = 1
                dnodal(auxi*2-1) = auxr
            else if (texto=='Y') then
                rest(auxi*2) = 1
                dnodal(auxi*2) = auxr
            else
                rest(auxi*2-1) = 1
                rest(auxi*2) = 1
                dnodal(auxi*2-1) = auxr
                dnodal(auxi*2) = auxr
            end if
        end do

        close(10)    
    end subroutine acadm
    
    !-----------------------------------------------------------------------------------------------------------!
    
    end module leitura