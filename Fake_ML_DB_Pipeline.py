'''
Pipeline to:
1) Add Fakes
2) Run subtraction, Run ML and load DB

We want to do this for a small sample of images to get an idea of times

Make different fakes numbers

+---------+-----------------+
| Version | Number of Fakes |
+---------+-----------------+
|	V1	  |			0		|
|	V2	  |			50		|
|   V3	  |			75		|
|	V4	  |			100		|
|	V5	  |			150		|
|	V6	  |			200		|
+---------+-----------------+
'''

import time

t0=time.time()

from IceCube_Fakes_ML_DB_V1 import *

Run_All_V1()

from IceCube_Fakes_ML_DB_V2 import *

Run_All_V2()
'''
from IceCube_Fakes_ML_DB_V3 import *

Run_All_V3()

from IceCube_Fakes_ML_DB_V4 import *

Run_All_V4()

from IceCube_Fakes_ML_DB_V5 import *

Run_All_V5()

from IceCube_Fakes_ML_DB_V6 import *

Run_All_V6()

print 'V1,2,3,4,5,6 took ', time.time()-t0, 'seconds'
'''

print 'Took', time.time()-t0, 'seconds'