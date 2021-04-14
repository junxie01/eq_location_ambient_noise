SUBROUTINE init_random_seed(case1)
!c ... resets (or not) the pseudo random sequence for each run
!      REAL*8 X 
INTEGER case1
INTEGER time(8) 
INTEGER seed_size
INTEGER, ALLOCATABLE :: newseed(:)
         
! ...if the user that calls with case1 = 0 then random seed is called without 
!    arguments ...should give a random default for every run (good for debuging) 
IF(case1.eq.0)THEN
        CALL RANDOM_SEED()                                                    
ELSEIF(case1.eq.1)THEN
! ... if the user calls with case1 = 1 then a new random seed 
!    is used in every new run (good for actual runs)
        CALL RANDOM_SEED(SIZE=seed_size) ! Determine seed array size (which depends on the compiler).                                                                                 
        ALLOCATE(newseed(seed_size))
! Generate a new seed as a "random" event of date information.       
        CALL DATE_AND_TIME(values=time)     ! Get the current time 
        newseed = time(4) * (360000*time(5) + 6000*time(6) + 100*time(7) + time(8))  
        CALL RANDOM_SEED(PUT=newseed) 
ENDIF         
END SUBROUTINE
