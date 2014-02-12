      !
      ! Definition of the interpolation grid 
      !
  
      Module DefGrid

!     INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)  ! If not using
                                                         ! with Elmer   
      USE Types    ! If using with Elmer 
      
      Real, Parameter :: kmin=0.002_dp ! valeur de ki mimum
      Integer, Parameter :: Ndiv=30    ! Ndiv+2 Number of points along ik1
      Integer, Parameter :: Ntot=813   ! Total number of points
      Integer, Parameter :: NetaI=4878 ! 6*4884 length of EtaI
      Integer, Parameter, Dimension(32) :: Nk2 = & 
                (/ -1,  46,  93, 139, 183, 226, 267, 307, 345, 382,&
                 417, 451, 483, 514, 543, 571, 597, 622, 645, 667,& 
                 687, 706, 723, 739, 753, 766, 777, 787, 795, 802,&
                 807, 811/) 
      
      End Module DefGrid
      

