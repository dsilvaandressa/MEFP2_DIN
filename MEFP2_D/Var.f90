
    module VAR
    implicit none
    
    !leitura
    integer :: nnos, io, nel, nfnod, nbc, i, auxi, j
    character(len=3) :: texto
    real :: auxr
    real, allocatable, dimension(:) :: X, fnodal, dnodal, Y
    real, allocatable, dimension(:,:) :: prop  !prop = E v d THICK
    integer, allocatable, dimension(:) :: rest !restrição do nó (1 restrito, 0 livre)
    integer, allocatable, dimension(:,:) :: inc
    
    !funções de forma e derivadas
    real, allocatable, dimension(:) :: fis, dfis1, dfis2
    real :: xsi1, xsi2, xsi3
    
    !elemento
    real, allocatable, dimension(:,:) :: coord_el
    
    !pontos de integração
    real, allocatable, dimension(:,:) :: pintegra
    
    !procedimentos
    real, dimension(2,2) :: S, Id, maux1, E
    integer :: el, k
    real, dimension(10) :: X1, X2, Y1, Y2
    real :: G, v, alfa, beta
    real, dimension(2,2,20) :: DA1
    real, dimension(2,2) ::  A0, A1
    
    !corpo do programa - processamento
    real, allocatable, dimension(:) :: fint, dbm, ip, u
    real, allocatable, dimension(:,:) :: HESS
    integer :: pct, pc, it, pi, p, gl, glg, gl2, glg2, info

    real :: err, tol, hessianinha
    character(len=3) :: ep
    real, dimension(2,2) :: maux2, D2E, DS
    real, dimension(2,2,20) :: DE
    real, dimension(20) :: ccar
    real, dimension(12,20,20) :: TM
    integer, dimension(20) :: ind
    
    !pós-processamento
    real, dimension(12) :: SX,SY,SXY
    real, dimension(12,10) :: MFIs_pi, MFI
    real, dimension(10,3) :: S_nodal
    real, dimension(10,10) :: maux3
    real, dimension(10) :: ip2
    real, allocatable, dimension(:,:) :: tensoes
    real, dimension(2,2) :: A
    
    !dinâmica
    real, allocatable, dimension(:) :: ace, vel, Q, R
    real :: dt, t, gama_newmark, beta_newmark, b1, b2
    real, dimension(20) :: fnodal_a
    real, dimension(20,20) :: MASS
    real, allocatable, dimension(:,:) :: MASS_G, C_G
    
    contains
    subroutine alocar()
        allocate(X(2*nnos), Y(2*nnos))
        allocate(inc(nel,10))
        allocate(prop(nel,4))
        allocate(fnodal(2*nnos), dnodal(2*nnos))
        allocate(rest(2*nnos))
        allocate(fint(2*nnos), dbm(2*nnos), ip(2*nnos))
        allocate(HESS(2*nnos,2*nnos))
        allocate(u(2*nnos))
        allocate(tensoes(nnos,4))
        allocate(ace(2*nnos), vel(2*nnos))
        allocate(Q(2*nnos), R(2*nnos))
        allocate(MASS_G(2*nnos,2*nnos), C_G(2*nnos,2*nnos))

    end subroutine alocar
    
    end module 
    