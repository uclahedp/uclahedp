import os
import psutil

def mem():
    process = psutil.Process(os.getpid())
    print(str(round(process.memory_info().rss*1e-6,2)) + ' MB') 