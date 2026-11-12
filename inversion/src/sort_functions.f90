SUBROUTINE sort_shell(arr)
  INTEGER, DIMENSION(:), INTENT(INOUT) :: arr   
  INTEGER :: i, j
  REAL*8 :: v
  DO i = 2, SIZE(arr)
     v = arr(i)
     j = i - 1
     DO WHILE (j >= 1)
        IF (arr(j) <= v) EXIT
        arr(j+1) = arr(j)
        j = j - 1
     END DO
     arr(j+1) = v
  END DO
END SUBROUTINE sort_shell

SUBROUTINE sort_shell_real(arr)
  REAL*8, DIMENSION(:), INTENT(INOUT) :: arr   
  INTEGER :: i, j
  REAL*8 :: v
  DO i = 2, SIZE(arr)
     v = arr(i)
     j = i - 1
     DO WHILE (j >= 1)
        IF (arr(j) <= v) EXIT
        arr(j+1) = arr(j)
        j = j - 1
     END DO
     arr(j+1) = v
  END DO
END SUBROUTINE sort_shell_real

