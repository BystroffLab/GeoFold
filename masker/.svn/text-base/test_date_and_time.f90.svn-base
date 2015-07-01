          program test_time_and_date
              character(8)  :: date
              character(10) :: time
              character(5)  :: zone
              Character(len=3),dimension(12),parameter :: &
              month=(/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'/)
              integer,dimension(8) :: values
              ! using keyword arguments
              call date_and_time(date,time,zone,values)
              call date_and_time(DATE=date,ZONE=zone)
              call date_and_time(TIME=time)
              call date_and_time(VALUES=values)
              print '(a,2x,a,2x,a)', date, time, zone
              print '(8i5))', values
              print '(i2,"-",a3,"-",i4))', values(3),month(values(2)),values(1)
          end program test_time_and_date
