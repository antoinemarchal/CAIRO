!! Header/ender call module
module mod_start
  !! Header/ender call module
  implicit none
  
  private
  
  public :: header, ender
  
contains

  subroutine header()
    implicit none

    write(*,*) "-------------------------------------------------------------------------"
    call timestamp()
    write(*,*) ""
    write(*,*) "          ___           ___                       ___           ___      "
    write(*,*) "         /  /\         /  /\        ___          /  /\         /  /\     "
    write(*,*) "        /  /:/        /  /::\      /  /\        /  /::\       /  /::\    "
    write(*,*) "       /  /:/        /  /:/\:\    /  /:/       /  /:/\:\     /  /:/\:\   "
    write(*,*) "      /  /:/  ___   /  /:/~/::\  /__/::\      /  /:/~/:/    /  /:/  \:\  "
    write(*,*) "     /__/:/  /  /\ /__/:/ /:/\:\ \__\/\:\__  /__/:/ /:/___ /__/:/ \__\:\ "
    write(*,*) "     \  \:\ /  /:/ \  \:\/:/__\/    \  \:\/\ \  \:\/:::::/ \  \:\ /  /:/ "
    write(*,*) "      \  \:\  /:/   \  \::/          \__\::/  \  \::/~~~~   \  \:\  /:/  "
    write(*,*) "       \  \:\/:/     \  \:\          /__/:/    \  \:\        \  \:\/:/   "
    write(*,*) "        \  \::/       \  \:\         \__\/      \  \:\        \  \::/    "
    write(*,*) "         \__\/         \__\/                     \__\/         \__\/     "
    write(*,*) ""
    write(*,*) " Version 1.0.0"
    write(*,*) " CAIRO is released as open source code"
    ! write(*,*) " Check out the documentation: https://antoinemarchal.github.io/ROHSA/"
    write(*,*) ""
    write(*,*) "run: ./CAIRO parameters.txt"
    write(*,*) "-------------------------------------------------------------------------"
  end subroutine header

  
  subroutine ender()
    implicit none
    
    write(*,*) "##################################################################"
    call timestamp()
    write(*,*) " Terminate"
    write(*,*) "                        CAIRO ALGORITHM"                   
    write(*,*) ""
    write(*,*) "##################################################################"
  end subroutine ender

  
  subroutine timestamp ( )

    !*****************************************************************************80
    !
    !! TIMESTAMP prints the current YMDHMS date as a time stamp.
    !
    !  Example:
    !
    !    31 May 2001   9:45:54.872 AM
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 May 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    None
    !
    implicit none

    character ( len = 8 ) ampm
    integer ( kind = 4 ) d
    integer ( kind = 4 ) h
    integer ( kind = 4 ) m
    integer ( kind = 4 ) mm
    character ( len = 9 ), parameter, dimension(12) :: month = (/ &
         'January  ', 'February ', 'March    ', 'April    ', &
         'May      ', 'June     ', 'July     ', 'August   ', &
         'September', 'October  ', 'November ', 'December ' /)
    integer ( kind = 4 ) n
    integer ( kind = 4 ) s
    integer ( kind = 4 ) values(8)
    integer ( kind = 4 ) y

    call date_and_time ( values = values )

    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)

    if ( h < 12 ) then
       ampm = 'AM'
    else if ( h == 12 ) then
       if ( n == 0 .and. s == 0 ) then
          ampm = 'Noon'
       else
          ampm = 'PM'
       end if
    else
       h = h - 12
       if ( h < 12 ) then
          ampm = 'PM'
       else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
             ampm = 'Midnight'
          else
             ampm = 'AM'
          end if
       end if
    end if

    write (*, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
         d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

    return
  end subroutine timestamp
  
end module mod_start
