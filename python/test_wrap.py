import os.path
import sys
currentdir = os.path.dirname(os.path.abspath(__file__))
parentdir = os.path.dirname(currentdir)

sys.path.insert(0,parentdir+os.path.sep + "lib")

import dagmc_call_wrap as w

d=w.DataClass('msg')
print(d.message())


