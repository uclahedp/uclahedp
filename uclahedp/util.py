import os
import psutil
import time

def mem():
    process = psutil.Process(os.getpid())
    print(str(round(process.memory_info().rss*1e-6,2)) + ' MB') 
    

def timeTest(t0=None):
    if t0 is None:
        return time.time()
    else:
        deltaT = time.time() - t0
        print('Execution time: '  + str( round(deltaT,3 )  ) + ' s' )
        return deltaT